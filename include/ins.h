#include "ArduinoEigen/ArduinoEigenDense.h"
#include "utils.h"
#include <math.h>

class INS {
private:
  float psi0_ = 0;
  float theta0_ = 0;
  float phi0_ = 0;

  Eigen::Matrix3f Rnb0_;
  Eigen::Matrix3f Rb0b_;
  Eigen::Matrix3f Rnb_;
  
  Eigen::Vector2f v_eb_n_;
  
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

  Eigen::Vector3f ba_;
  Eigen::Vector3f bg_;

  float T_;

  float llh_scale_ = 1000;

  friend class KalmanFilter;

public:
  INS(float roll_rad, float pitch_rad, float yaw_rad, float vn, float ve,
     float lat, float lon, float sample_interval);
  
  void step(const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc);
};

