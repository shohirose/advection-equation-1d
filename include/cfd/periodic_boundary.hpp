#ifndef CFD_PERIODIC_BOUNDARY_HPP
#define CFD_PERIODIC_BOUNDARY_HPP

#include <assert.h>

#include <Eigen/Core>

#include "cfd/problem_parameters.hpp"

namespace cfd {

class PeriodicBoundary {
 public:
  PeriodicBoundary(int n_boundary_cells, int n_domain_cells)
      : n_boundary_cells_{n_boundary_cells}, n_domain_cells_{n_domain_cells} {}

  PeriodicBoundary(const ProblemParameters& params)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells} {}

  template <typename Derived>
  void apply(Eigen::MatrixBase<Derived>& u) const noexcept {
    assert(u.size() == (n_boundary_cells_ * 2 + n_domain_cells_));
    using Eigen::seqN;
    // Left boundary
    u.head(n_boundary_cells_) = u(seqN(n_domain_cells_, n_boundary_cells_));
    // Right boundary
    u.tail(n_boundary_cells_) = u(seqN(n_boundary_cells_, n_boundary_cells_));
  }

 private:
  int n_boundary_cells_;
  int n_domain_cells_;
};

}  // namespace cfd

#endif  // CFD_PERIODIC_BOUNDARY_HPP