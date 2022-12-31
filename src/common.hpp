#ifndef CFD_COMMON_HPP
#define CFD_COMMON_HPP

#include <Eigen/Core>

#include "cfd/problem_parameters.hpp"

namespace cfd {

ProblemParameters make_params() noexcept;

Eigen::VectorXd make_x(const ProblemParameters& params) noexcept;

Eigen::VectorXd make_sine_wave(const Eigen::VectorXd& x) noexcept;

Eigen::VectorXd make_pulse_wave(const Eigen::VectorXd& x) noexcept;

}  // namespace cfd

#endif  // CFD_COMMON_HPP