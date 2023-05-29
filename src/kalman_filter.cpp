#include "kalman_filter.h"
#include <math.h>

// /*constructor
//   @param x_0 initial state
// */
// KalmanFilter::KalmanFilter(Eigen::VectorXf x_0)
//     : _x_c_k_1(x_0) {
//   // process noise covariance
//   _Q << 0.5, 0, 0.01,
//         0, 0.5, 0.01,
//         0.01, 0.01, 0.3;
//   // measurement noise covariance
//   _R << 0.5, 0, 0.1,
//         0, 0.5, 0.1,
//         0.1, 0.1, 0.6;
//   // initial estimate error covariance
//   _P_k_1 = Eigen::Matrix3f::Identity();
// }

KalmanFilter::KalmanFilter(INS& ins) : _ins(ins) {
  // initialize error states
  _dspi_nb << 0, 0, 0;
  _dv_eb_n << 0, 0;
  _dllh << 0, 0;
  _ba << 0, 0, 0;
  _bg << 0, 0, 0;

  // initialize covariance of estimate
  _P_covar.resize(13, 13);
  _P_covar.setZero();
  // roll pitch yaw
  _P_covar(0, 0) = deg2rad(10) * deg2rad(10);
  _P_covar(1, 1) = deg2rad(10) * deg2rad(10);
  _P_covar(2, 2) = deg2rad(10) * deg2rad(10);
  // velocity north and east
  _P_covar(3, 3) = 0.5 * 0.5;
  _P_covar(4, 4) = 0.5 * 0.5;
  // position latitude and longitude
  _P_covar(5, 5) =
      (2.5 * _ins.llh_scale / meridionalRadius(_ins._lat0 / _ins.llh_scale));
  _P_covar(6, 6) =
      (4.0 * _ins.llh_scale / transverseRadius(_ins._lat0 / _ins.llh_scale) /
       cosf(_ins._lat0 / _ins.llh_scale));
  _P_covar(5, 5) = _P_covar(5, 5) * _P_covar(5, 5);
  _P_covar(6, 6) = _P_covar(6, 6) * _P_covar(6, 6);
  // accelerometer bias
  _P_covar(7, 7) = (0.04e-3 * 10) * (0.04e-3 * 10);  // (0.04 mg)^2
  _P_covar(8, 8) = (0.04e-3 * 10) * (0.04e-3 * 10);  // (0.04 mg)^2
  _P_covar(9, 9) = (0.04e-3 * 10) * (0.04e-3 * 10);  // (0.04 mg)^2
  // gyro bias
  _P_covar(10, 10) = (deg2rad(10) / 3600) * (deg2rad(10) / 3600);  // (10 deg/h)^2
  _P_covar(11, 11) = (deg2rad(10) / 3600) * (deg2rad(10) / 3600);  // (10 deg/h)^2
  _P_covar(12, 12) = (deg2rad(10) / 3600) * (deg2rad(10) / 3600);  // (10 deg/h)^2

  // initialize measurement matrix
  _H.resize(4, 13);
  _H << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;

  // initialize system noise
  _Q.resize(13, 13);
  _Q.setZero();
  // gyro noise 0.007 deg/s/sqr(Hz)
  _Q.block(0, 0, 3, 3) = (3 * deg2rad(0.007)) * (3 * deg2rad(0.007)) * _ins._T *
                         Eigen::Matrix3f::Identity();
  // accelorometer noise 120 ug/sqr(Hz)
  _Q.block(3, 3, 2, 2) = (3 * 120e-6 * 10) * (3 * 120e-6 * 10) * _ins._T *
                         Eigen::Matrix2f::Identity();
  // accelerometer dynamic bias 0.04 mg
  _Q.block(7, 7, 3, 3) =
      (3 * 0.04e-3 * 10) * (3 * 0.04e-3 * 10) * Eigen::Matrix3f::Identity();
  // gyro dynamic bias 10 deg/h
  _Q.block(10, 10, 3, 3) = (3 * deg2rad(10) / 3600) * (3 * deg2rad(10) / 3600) *
                           Eigen::Matrix3f::Identity();

  // initialize measurement noise
  _R.resize(4, 4);
  _R << 0.05, 0, 0, 0,
        0, 0.05, 0, 0,
        0, 0, 2.5 * _ins.llh_scale / meridionalRadius(_ins._lat / _ins.llh_scale), 0,
        0, 0, 0, 4.0 * _ins.llh_scale / transverseRadius(_ins._lat / _ins.llh_scale) / cosf(_ins._lat / _ins.llh_scale);
  _R = _R * _R;
}

/*state function of the system*/
void KalmanFilter::stateFcn() {
  _x_p_k(0) = _x_c_k_1(0);
  _x_p_k(1) = _x_c_k_1(1);
  _x_p_k(2) = std::atan2(_x_c_k_1(1), _x_c_k_1(0));
}

/*Jacobian of the state function of the system*/
void KalmanFilter::stateJacFcn() {
  _A_k << 1, 0, 0,
          0, 1, 0,
          -_x_c_k_1(1)/(_x_c_k_1(0)*_x_c_k_1(0) + _x_c_k_1(1)*_x_c_k_1(1)),
          _x_c_k_1(0)/(_x_c_k_1(0)*_x_c_k_1(0) + _x_c_k_1(1)*_x_c_k_1(1)),
          0;
}

/*measurement function of the system*/
void KalmanFilter::measureFcn() {
  _z_p_k = _x_p_k;
}

/*Jacobian of the measurement function of the system*/
void KalmanFilter::measureJacFcn() {
  _H_k = Eigen::Matrix3f::Identity();
}

/*do prediction*/
void KalmanFilter::predict() {
  // calculate _x_p_k
  stateFcn();
  // calculate _A_k
  stateJacFcn();
  // calculate _P_p_k
  _P_p_k = _A_k * _P_k_1 * _A_k.transpose() + _Q;
}

/*do correction*/
void KalmanFilter::correct() {
  // calculate _H_k
  measureJacFcn();
  // calculate _K_k
  _K_k = _P_p_k * _H_k.transpose() *
         (_H_k * _P_p_k * _H_k.transpose() + _R).inverse();
  // calculate _z_p_k
  measureFcn();
  // calculate _x_c_k
  _x_c_k = _x_p_k + _K_k * (_z_k - _z_p_k);
  // calculate _P_k
  _P_k = (Eigen::Matrix3f::Identity() - _K_k * _H_k) * _P_p_k;
}

/*do estimation for the current step
  @param z_k measurement
  @returns estimate of state
*/
Eigen::Vector3f KalmanFilter::estimate(Eigen::Vector3f z_k) {
  _z_k = z_k;
  predict();
  correct();
  _x_c_k_1 = _x_c_k;
  _P_k_1 = _P_k;
  return _x_c_k;
}