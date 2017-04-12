#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

public:
        // the current NIS for radar
        double NIS_radar_;

        // the current NIS for laser
        double NIS_laser_;

        // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
        VectorXd x_;

        /**
         * Constructor
         */
        UKF();

        /**
         * Destructor
         */
        virtual ~UKF();

        /**
         * ProcessMeasurement
         * @param meas_package The latest measurement data of either radar or laser
         */
        void ProcessMeasurement(MeasurementPackage meas_package);

        /**
         * Prediction Predicts sigma points, the state, and the state covariance
         * matrix
         * @param delta_t Time between k and k+1 in s
         */
        void Prediction(double delta_t);

        /**
         * Updates the state and the state covariance matrix using a laser measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateLidar(MeasurementPackage meas_package);

        /**
         * Updates the state and the state covariance matrix using a radar measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateRadar(MeasurementPackage meas_package);


private:
        // initially set to false, set to true in first call of ProcessMeasurement
        bool is_initialized_;

        // if this is false, laser measurements will be ignored (except for init)
        bool use_laser_;

        // if this is false, radar measurements will be ignored (except for init)
        bool use_radar_;

        // state covariance matrix
        MatrixXd P_;

        // predicted sigma points matrix
        MatrixXd Xsig_pred_;

        // R matrix for lidar measurements
        MatrixXd R_lidar_;

        // R matrix for radar measurements
        MatrixXd R_radar_;

        // time when the state is true, in us
        long long time_us_;

        // Process noise standard deviation longitudinal acceleration in m/s^2
        double std_a_;

        // Process noise standard deviation yaw acceleration in rad/s^2
        double std_yawdd_;

        // Weights of sigma points
        VectorXd weights_;

        // State dimension
        int n_z_radar_;

        // State dimension
        int n_z_lidar_;

        // State dimension
        int n_x_;

        // Augmented state dimension
        int n_aug_;

        // Augmented state dimension
        int n_sig_;

        // Sigma point spreading parameter
        double lambda_;

        /**
         * Initialization of the unscented kalman filter
         *
         * @param m_pack
         */
        void Init(const MeasurementPackage m_pack);

        /**
         * Generate sigma points method simplifying kalman logic
         */
        MatrixXd GenerateSigmaPoints();

        /**
         * Sigma point prediction using the augmented matrix and the delta time
         *
         * @param Xsig_aug
         * @param dt
         */
        void SigmaPointPrediction(const MatrixXd &Xsig_aug, float dt);

        /**
         * Predict mean and covariance separating logic
         */
        void PredictMeanAndCovariance();

        /**
         * Predict radar measurement separating logic
         */
        void PredictRadarMeasurement();

        /**
         *
         * @param timestamp The timestamp as passed by the measurement package
         */
        float ProcessTimestamp(const long timestamp);
};

#endif /* UKF_H */
