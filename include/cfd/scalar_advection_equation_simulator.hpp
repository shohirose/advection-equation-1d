#ifndef CFD_SCALAR_ADVECTION_EQUATION_SIMULATOR_HPP
#define CFD_SCALAR_ADVECTION_EQUATION_SIMULATOR_HPP

#include <Eigen/Core>

#include "cfd/periodic_boundary.hpp"
#include "cfd/problem_parameters.hpp"

namespace cfd {

template <typename RiemannSolver, typename SpacialReconstructor,
          typename TimeIntegrator>
class ScalarAdvectionEquationSimulator {
 public:
  /**
   * @brief Construct a new Scalar Advection Equation Simulator object
   *
   * @param params Problem parameters
   * @param solver Riemann solver
   * @param reconstructor Spacial reconstructor
   * @param integrator Time integrator
   */
  ScalarAdvectionEquationSimulator(const ProblemParameters& params,
                                   const RiemannSolver& solver,
                                   const SpacialReconstructor& reconstructor,
                                   const TimeIntegrator& integrator)
      : n_boundary_cells_{params.n_boundary_cells},
        n_domain_cells_{params.n_domain_cells},
        n_timesteps_{params.n_timesteps},
        solver_{solver},
        reconstructor_{reconstructor},
        integrator_{integrator},
        boundary_{params} {}

  /**
   * @brief Run simulator
   *
   * @tparam Derived
   * @param u0 Initial condition
   * @return Eigen::VectorXd Values at the end of time steps.
   */
  template <typename Derived>
  Eigen::VectorXd run(const Eigen::MatrixBase<Derived>& u0) const noexcept {
    using Eigen::seqN;
    using Eigen::VectorXd;

    VectorXd u(this->n_total_cells());
    u(seqN(n_boundary_cells_, n_domain_cells_)) = u0;
    boundary_.apply(u);

    for (int i = 1; i <= n_timesteps_; ++i) {
      const VectorXd ul = reconstructor_.calc_left(u);
      const VectorXd ur = reconstructor_.calc_right(u);
      const VectorXd f = solver_.calc_flux(ul, ur);
      integrator_.update(u, f);
      boundary_.apply(u);
    }

    return u(seqN(n_boundary_cells_, n_domain_cells_));
  }

 private:
  int n_total_cells() const noexcept {
    return n_boundary_cells_ * 2 + n_domain_cells_;
  }

  int n_boundary_cells_;
  int n_domain_cells_;
  int n_timesteps_;
  RiemannSolver solver_;
  SpacialReconstructor reconstructor_;
  TimeIntegrator integrator_;
  PeriodicBoundary boundary_;
};

template <typename RiemannSolver, typename SpacialReconstructor,
          typename TimeIntegrator>
ScalarAdvectionEquationSimulator<RiemannSolver, SpacialReconstructor,
                                 TimeIntegrator>
make_simulator(const ProblemParameters& params, const RiemannSolver& solver,
               const SpacialReconstructor& reconstructor,
               const TimeIntegrator& integrator) {
  return {params, solver, reconstructor, integrator};
}

}  // namespace cfd

#endif  // CFD_SCALAR_ADVECTION_EQUATION_SIMULATOR_HPP