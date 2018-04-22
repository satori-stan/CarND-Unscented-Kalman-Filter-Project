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
    // We don't initialize time_us_ because it is done with the initial measurement
    // Process noise standard deviation longitudinal acceleration in m/s^2
    // TODO: Review initialization, should probably be a constructor param
    std_a_(3), // for a max acceleration of 6 m/s2
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
    // Initial state vector
    x_(VectorXd::Zero(n_x_)),
    // initial covariance matrix
    P_(MatrixXd::Identity(n_x_, n_x_)),
    // Predicted sigma points matrix
    // Generic initialization since values are reset every loop
    Xsig_pred_(MatrixXd(n_x_, n_sig_)),
    // Weights of sigma points
    // Generic initialization since values are reset every loop
    weights_(VectorXd(n_sig_)),

    // DO NOT MODIFY measurement noise values below these are provided by the
    // sensor manufacturer.
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

  // vector for weights
  weights_.setConstant(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
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
      x_ << ro * phi_cos, ro * phi_sin, 0, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_(0),
            measurement_pack.raw_measurements_(1),
            0,
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

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
      use_radar_) {
    // Radar updates
    UpdateRadar(measurement_pack);

  } else if (use_laser_) {
    // Laser updates
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
  // 1. Generate Sigma points for current state (previous observation)

  // augmented mean vector
  int aug_n = n_aug_ - n_x_;
  VectorXd x_aug = VectorXd(7);
  x_aug.head(n_x_) = x_;

  VectorXd nu = VectorXd::Zero(aug_n);  // Noise is zero mean
  x_aug.tail(aug_n) = nu;

  // augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
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

  // set sigma points as columns of matrix Xsig
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = (factor*A).colwise() + x_aug;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) =
      (-1 * (factor*A)).colwise() + x_aug;

  // 2. Predict Sigma points for next state (current observation)

  // create matrix with predicted sigma points as columns
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

  // predicted mean state
  x_ = (Xsig_pred_ * weights_).rowwise().sum();

  // predicted state covariance matrix
  // TODO: Do this with Eigen methods instead of looping (assuming it is faster)
  //  Something like P = ((Xsig_pred.colwise() - x)*(Xsig_pred.colwise() - x).transpose() * weights); //.rowwise().sum();
  P_.setZero();
  for (int i = 0; i < n_sig_; ++i) {
    MatrixXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
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

  // 2. Update state
  // calculate innovation covariance matrix S
  // matrix for cross correlation of sigma points Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  for (int i = 0; i < n_sig_; ++i) {
    MatrixXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    MatrixXd z_diff_t = z_diff.transpose();
    S += weights_(i) * z_diff * z_diff_t;
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * z_diff_t;
  }
  S += R;

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // update state mean and covariance matrix
  VectorXd y = (meas_package.raw_measurements_ - z_pred);
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  std::cout << "Lidar NIS: " << y.transpose() * S.inverse() * y << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
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

  // 2. Update state
  // calculate innovation covariance matrix S
  // matrix for cross correlation of sigma points Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  for (int i = 0; i < n_sig_; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    x_diff(3) = NormalizeAngle(x_diff(3));
    z_diff(1) = NormalizeAngle(z_diff(1));

    MatrixXd z_diff_t = z_diff.transpose();
    S += weights_(i) * z_diff * z_diff_t;
    Tc +=  weights_(i) * x_diff * z_diff_t;
  }
  S += R;

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // update state mean and covariance matrix
  VectorXd y = (meas_package.raw_measurements_ - z_pred);
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  std::cout << "Radar NIS: " << y.transpose() * S.inverse() * y << std::endl;
}

/**
 * Calculates an angle normalized to -PI and PI
 * @param theta The angle in radians
 * @returns The angle normalized between -PI and PI
 */
double UKF::NormalizeAngle(double theta) {
  // From a hint by fernandodamasio 
  // https://discussions.udacity.com/t/already-used-atan2-to-calculate-phi-in-hx-do-i-still-need-to-normalize-the-phi-in-y/242332/7
  return atan2(sin(theta), cos(theta));
  // From a StackOverflow answer https://stackoverflow.com/a/24234924
  // return theta - (k2PI * floor((theta + kPI) / k2PI));
}