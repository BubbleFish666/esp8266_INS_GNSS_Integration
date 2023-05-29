#include "ArduinoEigen/ArduinoEigenDense.h"
#include <math.h>

float deg2rad(float deg) {
  return deg / 180 * M_PI;
}

float rad2deg(float rad) {
  return rad / M_PI * 180;
}

Eigen::Matrix3f R3(float psi_rad) {
  Eigen::Matrix3f rotation;
  rotation << cosf(psi_rad), -sinf(psi_rad), 0,
              sinf(psi_rad), cosf(psi_rad), 0,
              0, 0, 1;
  return rotation;
}

Eigen::Matrix3f R2(float theta_rad) {
  Eigen::Matrix3f rotation;
  rotation << cosf(theta_rad), 0, sinf(theta_rad),
              0, 1, 0,
              -sinf(theta_rad), 0, cosf(theta_rad);
  return rotation;
}

Eigen::Matrix3f R1(float phi_rad) {
  Eigen::Matrix3f rotation;
  rotation << 1, 0, 0,
              0, cosf(phi_rad), -sinf(phi_rad),
              0, sinf(phi_rad), cosf(phi_rad);
  return rotation;
}

float meridionalRadius(float lat) {
  float b = 6356752.3131;
  float a = 6378137;
  float e2 = (a * a - b * b) / (a * a);
  float m_radius = b * b / (a * powf(1 - e2 * sinf(lat) * sinf(lat), 1.5));
  return m_radius;
}

float transverseRadius(float lat) {
  float b = 6356752.3131;
  float a = 6378137;
  float e2 = (a * a - b * b) / (a * a);
  float t_radius = a / sqrt(1 - e2 * sinf(lat) * sinf(lat));
  return t_radius;
}
