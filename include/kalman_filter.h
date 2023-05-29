#include "ArduinoEigen/ArduinoEigenDense.h"
#include "utils.h"
#include "ins.h"

class KalmanFilter {
private:
  // INS to be corrected
  INS& _ins;
  
  // attitude errors
  Eigen::Vector3f _dspi_nb;

  // velocity errors
  Eigen::Vector2f _dv_eb_n;

  // position errors
  Eigen::Vector2f _dllh;

  // accelerometer and gyro errors
  Eigen::Vector3f _ba;
  Eigen::Vector3f _bg;

  // covariance of estimate
  Eigen::MatrixXf _P_covar;
  
  // measurement matrix
  Eigen::MatrixXf _H;

  // process (system) noise Q
  Eigen::MatrixXf _Q;

  // measurement noise R
  Eigen::MatrixXf _R;

  // prediction
  // transition matrix
  Eigen::MatrixXf _PHI;
  // // state estimate of last step
  // Eigen::VectorXf _x_c_k_1;
  // // estimate error covariance of last step
  // Eigen::MatrixXf _P_k_1;

  // correction
  // Kalman gain
  Eigen::MatrixXf _K;
  // // predicted state
  // Eigen::VectorXf _x_p_k;
  // // predicted measurement
  // Eigen::VectorXf _z_p_k;
  // // predicted estimate error covariance
  // Eigen::MatrixXf _P_p_k;
  // // measurement matrix
  // Eigen::MatrixXf _H_k;
  // // measurement
  // Eigen::VectorXf _z_k;

  // // estimate (output of correction)
  // // state estimate
  // Eigen::VectorXf _x_c_k;
  // // estimate error covariance
  // Eigen::MatrixXf _P_k;

  void stateFcn();
  void stateJacFcn();
  void measureFcn();
  void measureJacFcn();

public:
  // KalmanFilter(Eigen::VectorXf x_0);
  KalmanFilter(INS& ins);
  void predict();
  void correct();
  Eigen::VectorXf estimate(Eigen::VectorXf z_k);
};

