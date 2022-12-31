#ifndef CFD_SPACIAL_RECONSTRUCTION_SCHEMES_HPP
#define CFD_SPACIAL_RECONSTRUCTION_SCHEMES_HPP

#include <Eigen/Core>

#include "cfd/problem_parameters.hpp"

namespace cfd {

class FirstOrderSpacialReconstructor {
 public:
  FirstOrderSpacialReconstructor(const ProblemParameters& params)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells} {
    assert(params.n_boundary_cells >= 1 &&
           "First order spacial reconstruction requires (# of boundary cells "
           ">= 1).");
  }

  FirstOrderSpacialReconstructor(int n_boundary_cells, int n_domain_cells)
      : n_boundary_cells_{n_boundary_cells}, n_domain_cells_{n_domain_cells} {
    assert(n_boundary_cells >= 1 &&
           "First order spacial reconstruction requires (# of boundary cells "
           ">= 1).");
  }

  template <typename Derived>
  Eigen::VectorXd calc_left(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    return u(Eigen::seqN(n_boundary_cells_ - 1, n_domain_cells_ + 1));
  }

  template <typename Derived>
  Eigen::VectorXd calc_right(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    return u(Eigen::seqN(n_boundary_cells_, n_domain_cells_ + 1));
  }

 private:
  int n_boundary_cells_;
  int n_domain_cells_;
};

class BeamWarmingSpacialReconstructor {
 public:
  BeamWarmingSpacialReconstructor(const ProblemParameters& params)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells},
        dt_{params.dt},
        dx_{params.dx},
        velocity_{params.velocity} {
    assert(params.n_boundary_cells >= 2 &&
           "Beam-Warming method requires (# of boundary cells >= 2).");
  }

  template <typename Derived>
  Eigen::VectorXd calc_left(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    const Eigen::VectorXd delta =
        0.5 * (u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1) -
               u.segment(n_boundary_cells_ - 2, n_domain_cells_ + 1));
    return u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1) +
           (1 - velocity_ * dt_ / dx_) * delta;
  }

  template <typename Derived>
  Eigen::VectorXd calc_right(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    const Eigen::VectorXd delta =
        0.5 * (u.segment(n_boundary_cells_ + 1, n_domain_cells_ + 1) -
               u.segment(n_boundary_cells_, n_domain_cells_ + 1));
    return u.segment(n_boundary_cells_, n_domain_cells_ + 1) -
           (1 + velocity_ * dt_ / dx_) * delta;
  }

 private:
  int n_boundary_cells_;
  int n_domain_cells_;
  double dt_;
  double dx_;
  double velocity_;
};

class FrommSpacialReconstructor {
 public:
  FrommSpacialReconstructor(const ProblemParameters& params)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells},
        dt_{params.dt},
        dx_{params.dx},
        velocity_{params.velocity} {
    assert(params.n_boundary_cells >= 2 &&
           "Fromm method requires (# of boundary cells >= 2).");
  }

  template <typename Derived>
  Eigen::VectorXd calc_left(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    const Eigen::VectorXd delta =
        0.25 * (u.segment(n_boundary_cells_, n_domain_cells_ + 1) -
                u.segment(n_boundary_cells_ - 2, n_domain_cells_ + 1));
    return u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1) +
           (1 - velocity_ * dt_ / dx_) * delta;
  }

  template <typename Derived>
  Eigen::VectorXd calc_right(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    const Eigen::VectorXd delta =
        0.25 * (u.segment(n_boundary_cells_ + 1, n_domain_cells_ + 1) -
                u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1));
    return u.segment(n_boundary_cells_, n_domain_cells_ + 1) -
           (1 + velocity_ * dt_ / dx_) * delta;
  }

 private:
  int n_boundary_cells_;
  int n_domain_cells_;
  double dt_;
  double dx_;
  double velocity_;
};

class LaxWendroffSpacialReconstructor {
 public:
  LaxWendroffSpacialReconstructor(const ProblemParameters& params)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells},
        dt_{params.dt},
        dx_{params.dx},
        velocity_{params.velocity} {
    assert(params.n_boundary_cells >= 1 &&
           "Beam-Warming method requires (# of boundary cells >= 1).");
  }

  template <typename Derived>
  Eigen::VectorXd calc_left(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    const Eigen::VectorXd delta =
        0.5 * (u.segment(n_boundary_cells_, n_domain_cells_ + 1) -
               u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1));
    return u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1) +
           (1 - velocity_ * dt_ / dx_) * delta;
  }

  template <typename Derived>
  Eigen::VectorXd calc_right(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    const Eigen::VectorXd delta =
        0.5 * (u.segment(n_boundary_cells_, n_domain_cells_ + 1) -
               u.segment(n_boundary_cells_ - 1, n_domain_cells_ + 1));
    return u.segment(n_boundary_cells_, n_domain_cells_ + 1) -
           (1 + velocity_ * dt_ / dx_) * delta;
  }

 private:
  int n_boundary_cells_;
  int n_domain_cells_;
  double dt_;
  double dx_;
  double velocity_;
};

template <typename SlopeLimiter>
class TvdSpacialReconstructor {
 public:
  TvdSpacialReconstructor(const ProblemParameters& params)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells},
        dt_{params.dt},
        dx_{params.dx},
        velocity_{params.velocity} {
    assert(params.n_boundary_cells >= 2 &&
           "TVD method requires (# of boundary cells >= 2).");
  }

  template <typename Derived>
  Eigen::VectorXd calc_left(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    const auto nb = n_boundary_cells_;
    const auto nd = n_domain_cells_;
    assert(u.size() == (nb * 2 + nd));
    using Eigen::lastN;
    using Eigen::seqN;
    using Eigen::VectorXd;
    const VectorXd du = u(seqN(nb - 1, nd + 2)) - u(seqN(nb - 2, nd + 2));
    const VectorXd r =
        du(seqN(0, nd + 1))
            .cwiseQuotient(du(lastN(nd + 1)) +
                           du(lastN(nd + 1)).unaryExpr([](double x) {
                             return x >= 0 ? 1e-5 : -1e-5;
                           }));
    const VectorXd phi = SlopeLimiter::eval(r);
    const VectorXd delta =
        0.5 * (u(seqN(nb, nd + 1)) - u(seqN(nb - 1, nd + 1)));
    return u(seqN(nb - 1, nd + 1)) +
           (1 - velocity_ * dt_ / dx_) * phi.cwiseProduct(delta);
  }

  template <typename Derived>
  Eigen::VectorXd calc_right(
      const Eigen::MatrixBase<Derived>& u) const noexcept {
    const auto nb = n_boundary_cells_;
    const auto nd = n_domain_cells_;
    assert(u.size() == (nb * 2 + nd));
    using Eigen::lastN;
    using Eigen::seqN;
    using Eigen::VectorXd;
    const VectorXd du = u(seqN(nb, nd + 2)) - u(seqN(nb - 1, nd + 2));
    const VectorXd r =
        du(lastN(nd + 1))
            .cwiseQuotient(du(seqN(0, nd + 1)) +
                           du(seqN(0, nd + 1)).unaryExpr([](double x) {
                             return x >= 0 ? 1e-5 : -1e-5;
                           }));
    const VectorXd phi = SlopeLimiter::eval(r);
    const VectorXd delta =
        0.5 * (u(seqN(nb, nd + 1)) - u(seqN(nb - 1, nd + 1)));
    return u(seqN(nb, nd + 1)) -
           (1 + velocity_ * dt_ / dx_) * phi.cwiseProduct(delta);
  }

 private:
  int n_boundary_cells_;
  int n_domain_cells_;
  double dt_;
  double dx_;
  double velocity_;
};

}  // namespace cfd

#endif  // CFD_SPACIAL_RECONSTRUCTION_SCHEMES_HPP