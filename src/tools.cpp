#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0) {
    cout << "ERROR: The estimations vector is empty" << endl;
  } else {
    if (estimations.size() != ground_truth.size()) {
      cout << "ERROR: The number of estimations and ground truths is different"
          << endl;
    } else {
      uint32_t size = estimations.size();
      for (uint32_t i = 0; i < size; ++i) {
        /*
        rmse += static_cast<VectorXd>(
            (estimations[i] - ground_truth[i]).array().pow(2));
        */
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        rmse += residual;
      }
      rmse /= size;
      rmse = rmse.array().sqrt();
    }
  }
  return rmse;
}