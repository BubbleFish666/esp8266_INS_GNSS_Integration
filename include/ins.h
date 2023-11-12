#include "ArduinoEigen/ArduinoEigenDense.h"
#include "utils.h"
#include <math.h>

class INS {
private:
  // initial orientation
  float psi0_ = 0;
  float theta0_ = 0;
  float phi0_ = 0;

  Eigen::Matrix3f Rnb0_;
  Eigen::Matrix3f Rb0b_;
  Eigen::Matrix3f Rnb_;
  
  // velocity
  Eigen::Vector2f v_eb_n_;
  
  // position
  float lat0_;
  float lon0_;
  float lat_incre_total_ = 0;
  float lon_incre_total_ = 0;
  float lat_;
  float lon_;
  int lat0_int_;
  int lon0_int_;
  float lat0_frac_;
  float lon0_frac_;

  float h_ = 0;

  // accelerometer error and gyro error
  Eigen::Vector3f ba_;
  Eigen::Vector3f bg_;

  // sample time interval
  float T_;

  // accelerometer and gyro readings
  Eigen::Vector3f f_ib_b_;
  Eigen::Vector3f w_ib_b_;

  // scale the longitude and latitude to milli-radian
  float llh_scale_ = 1000;

  // the Kalman Filter used to aid the INS
  friend class KalmanFilter;

public:
  INS(float roll_rad, float pitch_rad, float yaw_rad, float vn, float ve,
     float lat, float lon, float sample_interval);
  
  void step(const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc);
};

