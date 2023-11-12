#include "ArduinoEigen/ArduinoEigenDense.h"
#include "utils.h"
#include "ins.h"

class KalmanFilter {
private:
  // INS to be corrected
  INS& ins_;
  
  // attitude errors
  Eigen::Vector3f dspi_nb_;

  // velocity errors
  Eigen::Vector2f dv_eb_n_;

  // position errors
  Eigen::Vector2f dllh_;

  // accelerometer and gyro errors
  Eigen::Vector3f ba_;
  Eigen::Vector3f bg_;

  // covariance of estimate
  Eigen::MatrixXf P_;
  
  // measurement matrix
  Eigen::MatrixXf H_;

  // process (system) noise Q
  Eigen::MatrixXf Q_;

  // measurement noise R
  Eigen::MatrixXf R_;

  // prediction
  // transition matrix
  Eigen::MatrixXf PHI_;
  // // state estimate of last step
  // Eigen::VectorXf _x_c_k_1;
  // // estimate error covariance of last step
  // Eigen::MatrixXf _P_k_1;

  // correction
  // Kalman gain
  Eigen::MatrixXf K_;
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

