#include <iostream>
#include "Eigen/Dense"
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() :
    is_initialized_(false),
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_(true),
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_(true),
    // TODO: Find out what we are doing with time_us_
    // Process noise standard deviation longitudinal acceleration in m/s^2
    // TODO: Review initialization, should probably be a constructor param
    std_a_(6),
    // Process noise standard deviation yaw acceleration in rad/s^2
    // TODO: Review initialization, should probably be a constructor param
    std_yawdd_(kPI / 16),
    // State dimension
    // TODO: Review initialization, since it depends on the model
    n_x_(5),  // number of degrees of freedom for state space
    // Augmented state dimension
    // TODO: Review initialization, since it depends on the model
    n_aug_(n_x_ + 2),
    // Number of sigma points
    n_sig_(2 * n_aug_ + 1),
    // Sigma point spreading parameter
    lambda_(3 - n_aug_),
    // initial state vector
    x_(VectorXd::Zero(n_x_)),
    // initial covariance matrix
    P_(MatrixXd::Identity(n_x_, n_x_)),
    // predicted sigma points matrix
    // Generic initialization since values are reset every loop
    Xsig_pred_(MatrixXd(n_x_, n_sig_)),
    // Weights of sigma points
    // Generic initialization since values are reset every loop
    weights_(VectorXd(n_sig_)),

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
    // Laser measurement noise standard deviation position1 in m
    std_laspx_(0.15),
    // Laser measurement noise standard deviation position2 in m
    std_laspy_(0.15),
    // Radar measurement noise standard deviation radius in m
    std_radr_(0.3),
    // Radar measurement noise standard deviation angle in rad
    std_radphi_(0.03),
    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_(0.3) {
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // TODO: Set initial values of x and P here
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  double dt = (measurement_pack.timestamp_ - time_us_)/1000000.;
  time_us_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      double ro_dot = measurement_pack.raw_measurements_(2);
      double phi_cos = cos(phi);
      double phi_sin = sin(phi);
      // The lectures say we shouldn't use ro_dot to initialize velocity.
      x_ << ro * phi_cos, ro * phi_sin, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_(0),
            measurement_pack.raw_measurements_(1),
            0,
            0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /*
  ekf_.F_(0, 2) = ekf_.F_(1, 3) = dt;

  const float noise_ax = 9.F;
  const float noise_ay = 9.F;

  double dt_4_by_4 = pow(dt, 4.0)/4.0;
  double dt_3_by_2 = pow(dt, 3.0)/2.0;
  double dt_2 = pow(dt, 2);
  ekf_.Q_ << dt_4_by_4*noise_ax, 0, dt_3_by_2*noise_ax, 0,
            0, dt_4_by_4*noise_ay, 0, dt_3_by_2*noise_ay,
            dt_3_by_2*noise_ax, 0, dt_2*noise_ax, 0,
            0, dt_3_by_2*noise_ay, 0, dt_2*noise_ay;
  */

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
      use_radar_) {
    // Radar updates
    /*
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    */
    UpdateRadar(measurement_pack);
  } else if (use_laser_) {
    // Laser updates
    /*
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    */
    UpdateLidar(measurement_pack);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // 1. Generate Sigma points for current state (previous observation)

  // augmented mean vector
  int aug_n = n_aug_ - n_x_;
  VectorXd x_aug = VectorXd(7);
  x_aug.head(n_x_) = x_;

  VectorXd nu = VectorXd::Zero(aug_n);  // Noise is zero mean
  x_aug.tail(aug_n) = nu;

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.setZero();
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;

  MatrixXd Q(aug_n, aug_n);
  Q << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
  P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;

  // calculate covariance square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;
  double factor = sqrt(lambda_ + n_aug_);
  //for (int i = 0; i < n_aug; i++)
  // set sigma points as columns of matrix Xsig
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = (factor*A).colwise() + x_aug;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) =
      (-1 * (factor*A)).colwise() + x_aug;

  // 2. Predict Sigma points for next state (current observation)

  // create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);
  Xsig_pred_.setZero();
  // predict sigma points
  for (int i = 0; i < n_sig_; ++i) {
    VectorXd current = Xsig_aug.col(i);
    VectorXd x = current.head(n_x_);
    VectorXd delta(n_x_);
    VectorXd noise(n_x_);
    // TODO: Get rid of magic numbers
    double px = current(0);
    double py = current(1);
    double v = current(2);
    double yaw = current(3);
    double yawd = current(4);
    double nu_a = current(5);
    double nu_ydd = current(6);
    // avoid division by zero
    if (yawd == 0) {  
      delta <<  v*cos(yaw)*delta_t,
                v*sin(yaw)*delta_t,
                0,
                yawd*delta_t,
                0;
    } else {
      double v_by_yd = v/yawd;
      delta <<  v_by_yd*(sin(yaw + yawd*delta_t)-sin(yaw)),
                v_by_yd*(-cos(yaw + yawd*delta_t)+cos(yaw)),
                0,
                yawd*delta_t,
                0;
    }
    double delta_t_2_by_2 = pow(delta_t, 2)/2;
    noise <<  delta_t_2_by_2 * cos(yaw) * nu_a,
              delta_t_2_by_2 * sin(yaw) * nu_a,
              delta_t * nu_a,
              delta_t_2_by_2 * nu_ydd,
              delta_t * nu_ydd;
    // write predicted sigma points into right column
    Xsig_pred_.col(i) = x + delta + noise;
  }

  // 3. Predict mean and covariance

  // vector for weights
  //VectorXd weights = VectorXd(n_sig_);
  weights_.setZero();
  weights_.setConstant(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // predicted mean state
  x_ = (Xsig_pred_ * weights_).rowwise().sum();

  // predicted state covariance matrix
  // TODO: Do this with Eigen methods instead of looping (assuming it is faster)
  //  Something like P = ((Xsig_pred.colwise() - x)*(Xsig_pred.colwise() - x).transpose() * weights); //.rowwise().sum();
  P_.setZero();
  for (int i = 0; i < n_sig_; ++i) {
    MatrixXd x_diff = Xsig_pred_.col(i) - x_;
    double theta = x_diff(3);
    x_diff(3) = theta - (k2PI * floor((theta + kPI) / k2PI));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // 1. Predict measurement
  // set measurement dimension, lidar can measure px and py
  int n_z = 2;

  // matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.setZero();
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.setZero();
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  // XXX: Or Identity times a std_err vector?
  R(0,0) = pow(std_laspx_, 2);
  R(1,1) = pow(std_laspy_, 2);

  // transform sigma points into measurement space
  //Zsig.row(0) = (Xsig_pred.row(0).array().pow(2) + Xsig_pred.row(1).array().pow(2)).matrix().cwiseSqrt();
  //Zsig.row(1) = (Xsig_pred.row(1).array() / Xsig_pred.row(0).array()).matrix().unaryExpr(std::ptr_fun(atan));
  for (int i = 0; i < n_sig_; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    Zsig(0, i) = px;
    Zsig(1, i) = py;
    // calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S
  for (int i = 0; i < n_sig_; ++i) {
    MatrixXd Y = Zsig.col(i) - z_pred;
    S += (weights_(i) * Y * Y.transpose()); // + R;
  }
  S += R;

  // 2. Update state
  // matrix for cross correlation of sigma points Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  //Tc =  (Xsig_pred.colwise() - x).rowwise() * weights * (Zsig.colwise() - z_pred).transpose();
  for (int i = 0; i < n_sig_; ++i) {
    Tc +=  weights_(i) * (Xsig_pred_.col(i) - x_) *
        (Zsig.col(i) - z_pred).transpose();
  }
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  x_ = x_ + K * (meas_package.raw_measurements_ - z_pred);
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // 1. Predict measurement
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.setZero();
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.setZero();
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  // XXX: Or Identity times a std_err vector?
  R(0,0) = pow(std_radr_, 2);
  R(1,1) = pow(std_radphi_, 2);
  R(2,2) = pow(std_radrd_, 2);

  // transform sigma points into measurement space
  //Zsig.row(0) = (Xsig_pred.row(0).array().pow(2) + Xsig_pred.row(1).array().pow(2)).matrix().cwiseSqrt();
  //Zsig.row(1) = (Xsig_pred.row(1).array() / Xsig_pred.row(0).array()).matrix().unaryExpr(std::ptr_fun(atan));
  for (int i = 0; i < n_sig_; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double rho = sqrt(pow(px, 2) + pow(py, 2));
    if (abs(rho) > 0.001) {
      Zsig(0, i) = rho;
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
      // calculate mean predicted measurement
      z_pred += weights_(i) * Zsig.col(i);
    } else {
      std::cout << "ERROR: Division by zero" << std::endl;
    }
  }

  // calculate innovation covariance matrix S
  for (int i = 0; i < n_sig_; ++i) {
    MatrixXd Y = Zsig.col(i) - z_pred;
    Y(1) = atan2(sin(Y(1)), cos(Y(1)));  // Normalize angle
    S += (weights_(i) * Y * Y.transpose()); // + R;
  }
  S += R;

  // 2. Update state
  // matrix for cross correlation of sigma points Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  //Tc =  (Xsig_pred.colwise() - x).rowwise() * weights * (Zsig.colwise() - z_pred).transpose();
  for (int i = 0; i < n_sig_; ++i) {
    Tc +=  weights_(i) * (Xsig_pred_.col(i) - x_) *
        (Zsig.col(i) - z_pred).transpose();
  }
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  x_ = x_ + K * (meas_package.raw_measurements_ - z_pred);
  P_ = P_ - K * S * K.transpose();

  // TODO: Calculate NIS
}
