#ifndef CFD_TIME_INTEGRATION_SCHEMES_HPP
#define CFD_TIME_INTEGRATION_SCHEMES_HPP

#include <Eigen/Core>

#include "cfd/problem_parameters.hpp"

namespace cfd {

/**
 * @brief Time integration with explicit Euler scheme
 *
 * The explicit Euler scheme for the linear scalar advection equation can be
 * expressed by
 * @f[
 * u_j^{n+1} = u_j^n - \frac{\Delta t}{\Delta x} \cdot
 *    (\hat{f}_{j+1/2} - \hat{f}_{j-1/2})
 * @f]
 * where @f$ \hat{f} @f$ is numerical flux at cell boundaries.
 */
class ExplicitEulerScheme {
 public:
  /**
   * @brief Construct a new Explicit Euler Scheme object
   *
   * @param params Problem parameters
   */
  ExplicitEulerScheme(const ProblemParameters& params)
      : dx_{params.dx},
        dt_{params.dt},
        n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells} {}

  /**
   * @brief Construct a new Explicit Euler Scheme object
   *
   * @param dx Grid length
   * @param dt Time step length
   * @param n_boundary_cells Number of boundary cells
   * @param n_domain_cells Number of domain cells
   */
  ExplicitEulerScheme(double dx, double dt, int n_boundary_cells,
                      int n_domain_cells)
      : dx_{dx},
        dt_{dt},
        n_boundary_cells_{n_boundary_cells},
        n_domain_cells_{n_domain_cells} {}

  /**
   * @brief Update @f$ u @f$
   *
   * @tparam Derived1
   * @tparam Derived2
   * @param u Variable to solve
   * @param f Numerical flux
   */
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