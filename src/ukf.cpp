#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
        // if this is false, laser measurements will be ignored (except during init)
        use_laser_ = true;

        // if this is false, radar measurements will be ignored (except during init)
        use_radar_ = true;

        // initial state vector
        x_ = VectorXd(5);

        // initial covariance matrix
        P_ = MatrixXd(5, 5);

        // Process noise standard deviation longitudinal acceleration in m/s^2
        std_a_ = 30;

        // Process noise standard deviation yaw acceleration in rad/s^2
        std_yawdd_ = 30;

        // Laser measurement noise standard deviation position1 in m
        std_laspx_ = 0.15;

        // Laser measurement noise standard deviation position2 in m
        std_laspy_ = 0.15;

        // Radar measurement noise standard deviation radius in m
        std_radr_ = 0.3;

        // Radar measurement noise standard deviation angle in rad
        std_radphi_ = 0.03;

        // Radar measurement noise standard deviation radius change in m/s
        std_radrd_ = 0.3;

        /**
           TODO:

           Complete the initialization. See ukf.h for other member properties.

           Hint: one or more values initialized above might be wildly off...
        */
}

UKF::~UKF() {}


void UKF::Init(const MeasurementPackage measurement_pack) {

        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        if (is_initialized_) {
                // only initialize once
                return;
        }

        float px, py, vx = 0, vy = 0, v = 0, yaw = 0, yaw_rate = 0;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

                // extract the RADAR measurements and convert from
                // Polar to Cartesian coordinates
                float range = measurement_pack.raw_measurements_[0];
                float bearing = measurement_pack.raw_measurements_[1];
                float range_rate = measurement_pack.raw_measurements_[2];

                // calculate position and velocity
                px = range * cos(bearing);
                py = range * sin(bearing);
                vx = range_rate * cos(bearing);
                vy = range_rate * sin(bearing);

                v = sqrt(vx*vx + vy*vy);

                if(vx != 0){
                        yaw = vy / vx;
                }

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

                // if it is laser, just grab the raw x, y coordinates
                px = measurement_pack.raw_measurements_[0];
                py = measurement_pack.raw_measurements_[1];
        }

        // set the state and state covariance matrices
        x_ << px , py , v, yaw, yaw_rate;
        P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
               -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
               0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
               -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
               -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

        // ensure we mark the timestamp
        time_us_ = measurement_pack.timestamp_;
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
 * @param {MeasurementPackage} measurement_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {


        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        UKF::Init(measurement_pack);

        /*****************************************************************************
         *  Prediction
         ****************************************************************************/
        float dt = ProcessTimestamp(measurement_pack.timestamp_);
        Prediction(dt);

        /*****************************************************************************
         *  Update
         ****************************************************************************/
        VectorXd z_pred;
        switch (measurement_pack.sensor_type_) {
        case MeasurementPackage::RADAR:
                UpdateRadar(measurement_pack);
                break;
        case MeasurementPackage::LASER:
                UpdateLidar(measurement_pack);
                break;
        default:
                break;
        }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
        /**
           TODO:

           ProcessStateTransMatrix(dt);
           ProcessCovarianceMatrix(dt);

           Complete this function! Estimate the object's location. Modify the state
           vector, x_. Predict sigma points, the state, and the state covariance matrix.
        */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} measurement_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_package) {
        /**
           TODO:

           Complete this function! Use lidar data to update the belief about the object's
           position. Modify the state vector, x_, and covariance, P_.

           You'll also need to calculate the lidar NIS.
        */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} measurement_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_package) {
        /**
           TODO:

           Complete this function! Use radar data to update the belief about the object's
           position. Modify the state vector, x_, and covariance, P_.

           You'll also need to calculate the radar NIS.
        */
}
