#ifndef CFD_TIME_INTEGRATION_SCHEMES_HPP
#define CFD_TIME_INTEGRATION_SCHEMES_HPP

#include <Eigen/Core>

#include "cfd/problem_parameters.hpp"

namespace cfd {

class ExplicitEulerScheme {
 public:
  ExplicitEulerScheme(const ProblemParameters& params)
      : dx_{params.dx},
        dt_{params.dt},
        n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells} {}

  ExplicitEulerScheme(double dx, double dt, int n_boundary_cells,
                      int n_domain_cells)
      : dx_{dx},
        dt_{dt},
        n_boundary_cells_{n_boundary_cells},
        n_domain_cells_{n_domain_cells} {}

  template <typename Derived1, typename Derived2>
  void update(Eigen::MatrixBase<Derived1>& u,
              const Eigen::MatrixBase<Derived2>& f) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    assert(f.size() == (n_domain_cells_ + 1));
    u.segment(n_boundary_cells_, n_domain_cells_) -=
        (dt_ / dx_) * (f.tail(n_domain_cells_) - f.head(n_domain_cells_));
  }

 private:
  double dx_;
  double dt_;
  int n_boundary_cells_;
  int n_domain_cells_;
};

}  // namespace cfd

#endif  // CFD_TIME_INTEGRATION_SCHEMES_HPP