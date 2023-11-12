#include "ArduinoEigen/ArduinoEigenDense.h"
#include "utils.h"
#include <math.h>

class INS {
private:
  float _psi0 = 0;
  float _theta0 = 0;
  float _phi0 = 0;

  Eigen::Matrix3f _Rnb0;
  Eigen::Matrix3f _Rb0b;
  Eigen::Matrix3f _Rnb;
  
  Eigen::Vector2f _v_eb_n;
  
  float _lat0;
  float _lon0;
  float _lat_incre_total = 0;
  float _lon_incre_total = 0;
  float _lat;
  float _lon;
  int _lat0_int;
  int _lon0_int;
  float _lat0_frac;
  float _lon0_frac;

  float _h = 0;

  Eigen::Vector3f _ba;
  Eigen::Vector3f _bg;

  float _T;

  float llh_scale = 1000;

  friend class KalmanFilter;

public:
  INS(float roll_rad, float pitch_rad, float yaw_rad, float vn, float ve,
     float lat, float lon, float sample_interval);
  
  void step(const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc);
};

