#ifndef CFD_PROBLEM_PARAMETERS_HPP
#define CFD_PROBLEM_PARAMETERS_HPP

namespace cfd {

/**
 * @brief Parameters to solve the scalar advection equation using the finite
 * difference method.
 *
 * 1-D cells are ordered as:
 *
 * +----+----+----+-- ... --+----+----+----+
 * |    |    |    |   ...   |    |    |    |
 * +----+----+----+-- ... --+----+----+----+
 * <-- nb --><------- nd -------><-- nb -->
 *
 * where nb is number of boundary cells, and nd is the number of domain cells.
 */
struct ProblemParameters {
  int n_timesteps;       ///> Number of time steps
  int n_domain_cells;    ///> Number of domain cells
  int n_boundary_cells;  ///> Number of boundary cells to add one side of the
                         ///> domain
  double dt;             ///> Time step length
  double dx;             ///> Cell length
  double velocity;       ///> Velocity
  double eps;            ///> Entropy fix parameter for Harten-Riemann solver

  /**
   * @brief Returns the number of total cells (domain + boundary cells)
   *
   * @return int
   */
  int n_total_cells() const noexcept {
    return n_boundary_cells * 2 + n_domain_cells;
  }
};

}  // namespace cfd

#endif  // CFD_PROBLEM_PARAMETERS_HPP