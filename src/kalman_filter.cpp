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

/*constructor
  @param ins the INS that needs to be aided by the Kalman Filter*/
KalmanFilter::KalmanFilter(INS& ins) : ins_(ins) {
  // initialize error states
  dspi_nb_ << 0, 0, 0;
  dv_eb_n_ << 0, 0;
  dllh_ << 0, 0;
  ba_ << 0, 0, 0;
  bg_ << 0, 0, 0;

  // initialize covariance of estimate
  P_.resize(13, 13);
  P_.setZero();
  // roll pitch yaw
  P_(0, 0) = deg2rad(10) * deg2rad(10);
  P_(1, 1) = deg2rad(10) * deg2rad(10);
  P_(2, 2) = deg2rad(10) * deg2rad(10);
  // velocity north and east
  P_(3, 3) = 0.5 * 0.5;
  P_(4, 4) = 0.5 * 0.5;
  // position latitude and longitude
  P_(5, 5) =
      (2.5 * ins_.llh_scale_ / meridionalRadius(ins_.lat0_ / ins_.llh_scale_));
  P_(6, 6) =
      (4.0 * ins_.llh_scale_ / transverseRadius(ins_.lat0_ / ins_.llh_scale_) /
       cosf(ins_.lat0_ / ins_.llh_scale_));
  P_(5, 5) = P_(5, 5) * P_(5, 5);
  P_(6, 6) = P_(6, 6) * P_(6, 6);
  // accelerometer bias
  P_(7, 7) = (0.04e-3 * 10) * (0.04e-3 * 10);  // (0.04 mg)^2
  P_(8, 8) = (0.04e-3 * 10) * (0.04e-3 * 10);  // (0.04 mg)^2
  P_(9, 9) = (0.04e-3 * 10) * (0.04e-3 * 10);  // (0.04 mg)^2
  // gyro bias
  P_(10, 10) = (deg2rad(10) / 3600) * (deg2rad(10) / 3600);  // (10 deg/h)^2
  P_(11, 11) = (deg2rad(10) / 3600) * (deg2rad(10) / 3600);  // (10 deg/h)^2
  P_(12, 12) = (deg2rad(10) / 3600) * (deg2rad(10) / 3600);  // (10 deg/h)^2

  // initialize measurement matrix
  H_.resize(4, 13);
  H_ << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;

  // initialize system noise
  Q_.resize(13, 13);
  Q_.setZero();
  // gyro noise 0.007 deg/s/sqr(Hz)
  Q_(Eigen::seq(0, 2), Eigen::seq(0, 2)) = (3 * deg2rad(0.007)) *
                                           (3 * deg2rad(0.007)) * ins_.T_ *
                                           Eigen::Matrix3f::Identity();
  // accelorometer noise 120 ug/sqr(Hz)
  Q_(Eigen::seq(3, 4), Eigen::seq(3, 4)) = (3 * 120e-6 * 10) *
                                           (3 * 120e-6 * 10) * ins_.T_ *
                                           Eigen::Matrix2f::Identity();
  // accelerometer dynamic bias 0.04 mg
  Q_(Eigen::seq(7, 9), Eigen::seq(7, 9)) =
      (3 * 0.04e-3 * 10) * (3 * 0.04e-3 * 10) * Eigen::Matrix3f::Identity();
  // gyro dynamic bias 10 deg/h
  Q_(Eigen::seq(10, 12), Eigen::seq(10, 12)) = (3 * deg2rad(10) / 3600) *
                                               (3 * deg2rad(10) / 3600) *
                                               Eigen::Matrix3f::Identity();

  // initialize measurement noise
  R_.resize(4, 4);
  R_ << 0.05, 0, 0, 0,
        0, 0.05, 0, 0,
        0, 0, 2.5 * ins_.llh_scale_ / meridionalRadius(ins_.lat_ / ins_.llh_scale_), 0,
        0, 0, 0, 4.0 * ins_.llh_scale_ / transverseRadius(ins_.lat_ / ins_.llh_scale_) / cosf(ins_.lat_ / ins_.llh_scale_);
  R_ = R_ * R_;

  // initialize transition matrix
  PHI_.resize(13, 13);
  PHI_.setZero();
  PHI_(Eigen::seq(0, 2), Eigen::seq(0, 2)) = Eigen::Matrix3f::Identity();
  PHI_(Eigen::seq(3, 4), Eigen::seq(3, 4)) = Eigen::Matrix2f::Identity();
  PHI_(Eigen::seq(5, 6), Eigen::seq(5, 6)) = Eigen::Matrix2f::Identity();
  PHI_(Eigen::seq(7, 9), Eigen::seq(7, 9)) = Eigen::Matrix3f::Identity();
  PHI_(Eigen::seq(10, 12), Eigen::seq(10, 12)) = Eigen::Matrix3f::Identity();
}

/*prediction (propagation) phase of the Kalman Filter*/
void KalmanFilter::predict() {
  // update the transition matrix PHI_
  // 0-2 rows, delta phi, delta theta, and delta psi
  auto f_ib_n = ins_.Rnb_ * ins_.f_ib_b_;
  PHI_(Eigen::seq(0, 2), Eigen::seq(10, 12)) = ins_.Rnb_ * ins_.T_;

  // 3-4 rows, delta velocity north and delta velocity east
  PHI_(Eigen::seq(3, 4), Eigen::seq(0, 2)) << 0, f_ib_n(2), -f_ib_n(1),
                                              -f_ib_n(2), 0, f_ib_n(0);
  PHI_(Eigen::seq(3, 4), Eigen::seq(0, 2)) *= ins_.T_;
  PHI_(Eigen::seq(3, 4), Eigen::seq(7, 9)) = ins_.Rnb_(Eigen::seq(0, 1), Eigen::seq(0, 2)) * ins_.T_;

  // 5-6 rows, delta latitude and delta longitude
  PHI_(Eigen::seq(5, 6), Eigen::seq(3, 4)) << 
    ins_.T_ * ins_.llh_scale_ / (meridionalRadius(ins_.lat_ / ins_.llh_scale_) + ins_.h_), 0,
    0, ins_.T_ * ins_.llh_scale_ / ((transverseRadius(ins_.lat_ / ins_.llh_scale_) + ins_.h_) * cosf(ins_.lat_ / ins_.llh_scale_));

  // 7-9 rows, accelerometer x-, y- and z-error
  // 10-12 rows, gyro x-, y-, and z-error
  // These 6 rows are constant so they are already initialized in constructor.
  
  // propagation
  P_ = PHI_ * P_ * PHI_.transpose() + Q_;
  
  // Note: the states are zeroed after each correction, so no need to propagate
  // them as they are zero during propagation anyway.
}

/*correction (measurement) phase of the Kalman Filter*/
void KalmanFilter::correct() {

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