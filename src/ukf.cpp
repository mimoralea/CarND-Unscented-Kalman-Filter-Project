#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl;

#define ZERO 1e-9

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
        // if this is false, laser measurements will be ignored (except during init)
        use_laser_ = true;

        // if this is false, radar measurements will be ignored (except during init)
        use_radar_ = true;

        //
        is_initialized_ = false;

        // Process noise standard deviation longitudinal acceleration in m/s^2
        std_a_ = 0.5;

        // Process noise standard deviation yaw acceleration in rad/s^2
        std_yawdd_ = 0.25;

        // Laser measurement noise standard deviation position1 in m
        std_laspx_ = 0.225;

        // Laser measurement noise standard deviation position2 in m
        std_laspy_ = 0.15;

        // Radar measurement noise standard deviation radius in m
        std_radr_ = 0.15;

        // Radar measurement noise standard deviation angle in rad
        std_radphi_ = 0.0015;

        // Radar measurement noise standard deviation radius change in m/s
        std_radrd_ = 0.08;

        //set state dimension
        n_x_ = 5;

        //set augmented dimension
        n_aug_ = n_x_ + 2;

        //set augmented dimension
        n_sig_ = 2 * n_aug_ + 1;

        //
        n_z_radar_ = 3;

        //
        n_z_lidar_ = 2;

        //define spreading parameter
        lambda_ = 3 - n_aug_;

        // initial state vector
        x_ = VectorXd::Zero(n_x_);

        // initial covariance matrix
        P_ = MatrixXd::Zero(n_x_, n_x_);

        Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_);

        weights_ = VectorXd::Zero(n_sig_);

        R_lidar_ = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);

        R_radar_ = MatrixXd::Zero(n_z_radar_, n_z_radar_);
}

UKF::~UKF() {}


void UKF::Init(const MeasurementPackage m_pack) {

        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        float px = 0, py = 0, vx = 0, vy = 0, v = 0, yaw = 0, yaw_rate = 0;
        if (m_pack.sensor_type_ == MeasurementPackage::RADAR) {

                // extract the RADAR measurements and convert from
                // Polar to Cartesian coordinates
                float range = m_pack.raw_measurements_[0];
                float bearing = m_pack.raw_measurements_[1];
                float range_rate = m_pack.raw_measurements_[2];

                // calculate position and velocity
                px = range * cos(bearing);
                py = range * sin(bearing);
                vx = range_rate * cos(bearing);
                vy = range_rate * sin(bearing);
                v = sqrt(vx*vx + vy*vy);
                yaw = vx == 0 ? yaw : vy / vx;

        } else if (m_pack.sensor_type_ == MeasurementPackage::LASER) {

                // if it is laser, just grab the raw x, y coordinates
                px = m_pack.raw_measurements_[0];
                py = m_pack.raw_measurements_[1];
        }

        // set the state and state covariance matrices
        x_ << px , py , v, yaw, yaw_rate;
        P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1000, 0, 0,
              0, 0, 0, 1000, 0,
              0, 0, 0, 0, 1;

        // set weights
        weights_.fill(0.5 / (n_aug_ + lambda_));
        weights_(0) = lambda_ / (lambda_ + n_aug_);

        // radar covariance matrix
        R_radar_ << std_radr_ * std_radr_ , 0, 0,
                    0, std_radphi_ * std_radphi_, 0,
                    0, 0, std_radrd_ * std_radrd_;

        // lidar covariance matrix
        R_lidar_ << std_laspx_ * std_laspx_ , 0,
                    0, std_laspy_ * std_laspy_;

        // ensure we mark the timestamp
        time_us_ = m_pack.timestamp_;
        is_initialized_ = true;
}

/**
 * Mark timestamp and calculate elapsed
 * seconds since last measurement
 */
float UKF::ProcessTimestamp(const long timestamp) {
        float dt = (timestamp - time_us_) / 1000000.0;
        time_us_ = timestamp;
        return dt;
}

/**
 * @param {MeasurementPackage} m_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage m_pack) {


        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        if(!is_initialized_) {
                UKF::Init(m_pack);
        }

        /*****************************************************************************
         *  Prediction
         ****************************************************************************/
        Prediction(ProcessTimestamp(m_pack.timestamp_));

        /*****************************************************************************
         *  Update
         ****************************************************************************/
        if (use_radar_ && m_pack.sensor_type_ == MeasurementPackage::RADAR) {
                UpdateRadar(m_pack);
        }
        else if (use_laser_ && m_pack.sensor_type_ == MeasurementPackage::LASER) {
                UpdateLidar(m_pack);
        }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
        // create augmented sigma points and prediction
        SigmaPointPrediction(GenerateSigmaPoints(), dt);

        // calculate predicted state and covariance
        PredictMeanAndCovariance();
}

MatrixXd UKF::GenerateSigmaPoints() {

        //create augmented mean vector
        VectorXd x_aug = VectorXd::Zero(n_aug_);

        //create augmented state covariance
        MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

        //create process noise covariance
        MatrixXd Q = MatrixXd::Zero(n_aug_ - n_x_, n_aug_ - n_x_);
        Q << std_a_ * std_a_, 0,
                0, std_yawdd_ * std_yawdd_;

        //create sigma point matrix
        MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

        //create augmented mean state
        x_aug << x_, 0, 0;

        //create augmented covariance matrix
        P_aug.topLeftCorner(n_x_, n_x_) = P_;
        P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q;

        //create square root matrix
        MatrixXd A_aug = P_aug.llt().matrixL();

        //create augmented sigma points
        MatrixXd term_aug = sqrt(lambda_ + n_aug_) * A_aug;
        Xsig_aug.col(0) = x_aug;
        Xsig_aug.block(0, 1, n_aug_, n_aug_) = term_aug.colwise() + x_aug;
        Xsig_aug.block(0, 1 + n_aug_, n_aug_, n_aug_) = (-1 * term_aug).colwise() + x_aug;
        return Xsig_aug;
}

void UKF::SigmaPointPrediction(const MatrixXd &Xsig_aug, float dt) {

        // neither 0 nor 1 are actually needed
        VectorXd vs = Xsig_aug.row(2);
        VectorXd yaws = Xsig_aug.row(3);
        VectorXd yawds = Xsig_aug.row(4);
        VectorXd nu_as = Xsig_aug.row(5);
        VectorXd nu_yawdds = Xsig_aug.row(6);

        ArrayXd new_yaws = (yaws + yawds * dt).array();
        ArrayXd ratios = vs.array() / yawds.array();
        MatrixXd dvalues = MatrixXd::Zero(n_x_, n_sig_);
        dvalues.row(0) = (yawds.array() > ZERO).select(
                ratios * (new_yaws.sin() - yaws.array().sin()),
                dt * vs.array() * yaws.array().cos());
        dvalues.row(1) = (yawds.array() > ZERO).select(
                ratios * (yaws.array().cos() - new_yaws.cos()),
                dt * vs.array() * yaws.array().sin());
        dvalues.row(3) = yawds * dt;

        MatrixXd noises = MatrixXd::Zero(n_x_, n_sig_);
        noises.row(0) = dt * dt * 0.5 * nu_as.array() * yaws.array().cos();
        noises.row(1) = dt * dt * 0.5 * nu_as.array() * yaws.array().sin();
        noises.row(2) = dt * nu_as;
        noises.row(3) = dt * dt * 0.5 * nu_yawdds;
        noises.row(4) = dt * nu_yawdds;

        Xsig_pred_ = Xsig_aug.topRows(n_x_) + dvalues + noises;
}

void UKF::PredictMeanAndCovariance() {

        //predicted state mean
        x_ = (Xsig_pred_ * weights_.asDiagonal()).rowwise().sum();

        //predicted state covariance matrix
        MatrixXd diff = Xsig_pred_.colwise() - x_;
        diff.row(3) = diff.row(3).array().unaryExpr(&normalize_angle);
        P_ = (diff * weights_.asDiagonal()) * diff.transpose();
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} m_pack
 */
void UKF::UpdateLidar(MeasurementPackage m_pack) {

        // PredictLidarMeasurement

        //create matrix for sigma points in measurement space
        MatrixXd Zsig = MatrixXd::Zero(n_z_lidar_, n_sig_);

        //mean predicted measurement
        VectorXd z_pred = VectorXd::Zero(n_z_lidar_);

        //measurement covariance matrix S
        MatrixXd S = MatrixXd::Zero(n_z_radar_, n_z_radar_);

        Zsig.row(0) = Xsig_pred_.row(0);
        Zsig.row(1) = Xsig_pred_.row(1);

        //predicted state mean
        z_pred = (Zsig * weights_.asDiagonal()).rowwise().sum();

        //measurement covariance matrix S
        MatrixXd diff = Zsig.colwise() - z_pred;
        diff.row(1) = diff.row(1).array().unaryExpr(&normalize_angle);
        S = (diff * weights_.asDiagonal()) * diff.transpose();
        S = S + R_lidar_;

        // UpdateState

        //calculate cross correlation matrix
        MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_lidar_);
        MatrixXd x_diffs = Xsig_pred_.colwise() - x_;
        MatrixXd z_diffs = Zsig.colwise() - z_pred;
        x_diffs.row(3) = x_diffs.row(3).array().unaryExpr(&normalize_angle);
        z_diffs.row(1) = z_diffs.row(1).array().unaryExpr(&normalize_angle);
        Tc = (x_diffs * weights_.asDiagonal()) * z_diffs.transpose();

        //Kalman gain K;
        MatrixXd K = Tc * S.inverse();

        //residual
        VectorXd z = VectorXd(n_z_lidar_);
        z << m_pack.raw_measurements_[0],
                m_pack.raw_measurements_[1];
        VectorXd z_diff = z - z_pred;

        //angle normalization
        z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

        //update state mean and covariance matrix
        x_ = x_ + K * z_diff;
        P_ = P_ - K * S * K.transpose();
        NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} m_pack
 */
void UKF::UpdateRadar(MeasurementPackage m_pack) {

        // PredictRadarMeasurement

        //create matrix for sigma points in measurement space
        MatrixXd Zsig = MatrixXd::Zero(n_z_radar_, n_sig_);

        //mean predicted measurement
        VectorXd z_pred = VectorXd::Zero(n_z_radar_);

        //measurement covariance matrix S
        MatrixXd S = MatrixXd::Zero(n_z_radar_, n_z_radar_);

        ArrayXd pxs = Xsig_pred_.row(0);
        ArrayXd pxs_sq = pxs * pxs;

        ArrayXd pys = Xsig_pred_.row(1);
        ArrayXd pys_sq = pys * pys;

        ArrayXd vs = Xsig_pred_.row(2);
        ArrayXd yaws = Xsig_pred_.row(3);
        ArrayXd vs1 = yaws.cos() * vs;
        ArrayXd vs2 = yaws.sin() * vs;

        ArrayXd rs = (pxs_sq + pys_sq).sqrt();
        rs = (rs > ZERO).select(rs, ZERO);
        Zsig.row(0) = rs;
        Zsig.row(1) = pys.binaryExpr(pxs, &atan2_m);
        Zsig.row(2) = (pxs * vs1 + pys * vs2) / rs;

        //predicted state mean
        z_pred = (Zsig * weights_.asDiagonal()).rowwise().sum();

        //measurement covariance matrix S
        MatrixXd diff = Zsig.colwise() - z_pred;
        diff.row(1) = diff.row(1).array().unaryExpr(&normalize_angle);
        S = (diff * weights_.asDiagonal()) * diff.transpose();
        S = S + R_radar_;

        // UpdateState

        //calculate cross correlation matrix
        MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_radar_);
        MatrixXd x_diffs = Xsig_pred_.colwise() - x_;
        MatrixXd z_diffs = Zsig.colwise() - z_pred;
        x_diffs.row(3) = x_diffs.row(3).array().unaryExpr(&normalize_angle);
        z_diffs.row(1) = z_diffs.row(1).array().unaryExpr(&normalize_angle);
        Tc = (x_diffs * weights_.asDiagonal()) * z_diffs.transpose();

        //Kalman gain K;
        MatrixXd K = Tc * S.inverse();

        //residual
        VectorXd z = VectorXd(n_z_radar_);
        z << m_pack.raw_measurements_[0],
             m_pack.raw_measurements_[1],
             m_pack.raw_measurements_[2];
        VectorXd z_diff = z - z_pred;

        //angle normalization
        z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

        //update state mean and covariance matrix
        x_ = x_ + K * z_diff;
        P_ = P_ - K * S * K.transpose();
        NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
