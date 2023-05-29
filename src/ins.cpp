#include "ins.h"

INS::INS(float roll_rad, float pitch_rad, float yaw_rad, float vn, float ve,
         float lat, float lon, float sample_interval) {
  // attitude
  _psi0 = yaw_rad;
  _theta0 = pitch_rad;
  _phi0 = roll_rad;

  _Rnb0 = R3(0.5 * M_PI) * R1(M_PI);
  _Rb0b = R3(_psi0) * R2(_theta0) * R1(_phi0);
  _Rnb = _Rnb0 * _Rb0b;

  // velocity
  _v_eb_n << vn, ve;

  // position
  _lat0 = lat;  // initial position
  _lon0 = lon;
  _lat_incre_total = 0;  // position increment
  _lon_incre_total = 0;
  _lat = lat;  // current position
  _lon = lon;
  _lat0_int = int(_lat0);  // initial position integer part
  _lon0_int = int(_lon0);
  _lat0_frac = _lat0 - _lat0_int;  // initial position fractional part
  _lon0_frac = _lon0 - _lon0_int;

  // sensor errors
  _ba << 0, 0, 0;
  _bg << 0, 0, 0;

  // sample rate
  _T = sample_interval;
}

void INS::step(const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc){
  // corrected angular velocity
  Eigen::Vector3f w_ib_b = gyro;
  w_ib_b += _bg;

  // corrected acceleration
  Eigen::Vector3f f_ib_b = acc;
  f_ib_b += _ba;

  // attitude
  Eigen::Matrix3f Rnb_ = _Rnb;
  Eigen::Matrix3f omega_ib_b;
  omega_ib_b << 0, -w_ib_b(2), w_ib_b(1),
                w_ib_b(2), 0, -w_ib_b(0),
                -w_ib_b(1), w_ib_b(0), 0;
  _Rnb = Rnb_ * (Eigen::Matrix3f::Identity() + omega_ib_b * _T);

  // velocity
  auto v_eb_n_ = _v_eb_n;
  auto v_incre = Rnb_ * f_ib_b * _T;
  _v_eb_n(0) += v_incre(0);
  _v_eb_n(1) += v_incre(1);

  // position
  auto lat_ = _lat;
  auto lat_incre = (v_eb_n_(0) / (meridionalRadius(lat_ / llh_scale) + _h) +
                    _v_eb_n(0) / (meridionalRadius(lat_ / llh_scale) + _h)) *
                   0.5 * _T * llh_scale;
  _lat_incre_total += lat_incre;
  _lat = _lat0 + _lat_incre_total;

  auto lon_incre = (v_eb_n_(1) / ((transverseRadius(lat_ / llh_scale) + _h) *
                                  cosf(lat_ / llh_scale)) +
                    _v_eb_n(1) / ((transverseRadius(_lat / llh_scale) + _h) *
                                  cosf(_lat / llh_scale))) *
                   0.5 * _T * llh_scale;
  _lon_incre_total += lon_incre;
  _lon = _lon0 + _lon_incre_total;
}



