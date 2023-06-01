#include "ins.h"

INS::INS(float roll_rad, float pitch_rad, float yaw_rad, float vn, float ve,
         float lat, float lon, float sample_interval) {
  // attitude
  psi0_ = yaw_rad;
  theta0_ = pitch_rad;
  phi0_ = roll_rad;

  Rnb0_ = R3(0.5 * M_PI) * R1(M_PI);
  Rb0b_ = R3(psi0_) * R2(theta0_) * R1(phi0_);
  Rnb_ = Rnb0_ * Rb0b_;

  // velocity
  v_eb_n_ << vn, ve;

  // position
  lat0_ = lat;  // initial position
  lon0_ = lon;
  lat_incre_total_ = 0;  // position increment
  lon_incre_total_ = 0;
  lat_ = lat;  // current position
  lon_ = lon;
  lat0_int_ = int(lat0_);  // initial position integer part
  lon0_int_ = int(lon0_);
  lat0_frac_ = lat0_ - lat0_int_;  // initial position fractional part
  lon0_frac_ = lon0_ - lon0_int_;

  // sensor errors
  ba_ << 0, 0, 0;
  bg_ << 0, 0, 0;

  // sample rate
  T_ = sample_interval;
}

void INS::step(const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc){
  // corrected angular velocity
  w_ib_b_ = gyro;
  w_ib_b_ += bg_;

  // corrected acceleration
  f_ib_b_ = acc;
  f_ib_b_ += ba_;

  // attitude
  Eigen::Matrix3f Rnb_ = Rnb_;
  Eigen::Matrix3f omega_ib_b;
  omega_ib_b << 0, -w_ib_b_(2), w_ib_b_(1),
                w_ib_b_(2), 0, -w_ib_b_(0),
                -w_ib_b_(1), w_ib_b_(0), 0;
  Rnb_ = Rnb_ * (Eigen::Matrix3f::Identity() + omega_ib_b * T_);

  // velocity
  auto v_eb_n_1 = v_eb_n_;
  auto v_incre = Rnb_ * f_ib_b_ * T_;
  v_eb_n_(0) += v_incre(0);
  v_eb_n_(1) += v_incre(1);

  // position
  auto lat_1 = lat_;
  auto lat_incre = (v_eb_n_1(0) / (meridionalRadius(lat_1 / llh_scale_) + h_) +
                    v_eb_n_(0) / (meridionalRadius(lat_1 / llh_scale_) + h_)) *
                   0.5 * T_ * llh_scale_;
  lat_incre_total_ += lat_incre;
  lat_ = lat0_ + lat_incre_total_;

  auto lon_incre = (v_eb_n_1(1) / ((transverseRadius(lat_1 / llh_scale_) + h_) *
                                  cosf(lat_1 / llh_scale_)) +
                    v_eb_n_(1) / ((transverseRadius(lat_ / llh_scale_) + h_) *
                                  cosf(lat_ / llh_scale_))) *
                   0.5 * T_ * llh_scale_;
  lon_incre_total_ += lon_incre;
  lon_ = lon0_ + lon_incre_total_;
}



