#include "common.hpp"

#include <cmath>

namespace cfd {

namespace {

constexpr double x_left() noexcept { return -1.0; }

constexpr double x_right() noexcept { return 1.0; }

}  // namespace

ProblemParameters make_params() noexcept {
  const auto xl = x_left();
  const auto xr = x_right();
  const int n_domain_cells = 100;
  const int n_boundary_cells = 2;
  const auto dx = (xr - xl) / static_cast<double>(n_domain_cells);
  const auto dt = 0.2 * dx;
  const int n_timesteps = 500;
  const double velocity = 1.0;

  return {n_timesteps, n_domain_cells, n_boundary_cells, dt, dx, velocity};
}

Eigen::VectorXd make_x(const ProblemParameters& params) noexcept {
  using Eigen::VectorXd;
  const auto nd = params.n_domain_cells;
  const VectorXd x = VectorXd::LinSpaced(nd + 1, x_left(), x_right());
  return 0.5 * (x.head(nd) + x.tail(nd));
}

Eigen::VectorXd make_sine_wave(const Eigen::VectorXd& x) noexcept {
  return ((2.0 * M_PI) * x.array()).sin().matrix();
}

Eigen::VectorXd make_pulse_wave(const Eigen::VectorXd& x) noexcept {
  using Eigen::VectorXd;
  VectorXd u = VectorXd::Zero(x.size());
  u(Eigen::seq(x.size() / 2 - 10, x.size() / 2 + 10)).array() = 1.0;
  return u;
}

}  // namespace cfd
