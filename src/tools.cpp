#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}
Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
        VectorXd rmse(4);
        rmse << 0,0,0,0;

        // make sure we have the same number of measurements
        if(estimations.size() != ground_truth.size()
           || estimations.size() == 0) {
                cout << "Invalid estimation or ground_truth data" << endl;
                return rmse;
        }

        // calculate the squares differences
        for(unsigned int i=0; i < estimations.size(); ++i) {
                VectorXd residual = estimations[i] - ground_truth[i];
                residual = residual.array() * residual.array();
                rmse += residual;
        }

        // average and take square root
        rmse = rmse / estimations.size();
        rmse = rmse.array().sqrt();
        return rmse;
}

float normalize_angle(float angle) {
        return atan2(sin(angle), cos(angle));
}

float atan2_m(float y, float x) {
        return atan2(y, x);
}
