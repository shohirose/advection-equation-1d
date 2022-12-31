#ifndef CFD_RIEMANN_SOLVERS_HPP
#define CFD_RIEMANN_SOLVERS_HPP

#include <Eigen/Core>

namespace cfd {

class RoeRiemannSolver {
 public:
  RoeRiemannSolver(double velocity) : velocity_{velocity} {}

  template <typename Derived1, typename Derived2>
  Eigen::VectorXd calc_flux(
      const Eigen::MatrixBase<Derived1>& ul,
      const Eigen::MatrixBase<Derived2>& ur) const noexcept {
    return 0.5 * (velocity_ * (ul + ur) - std::fabs(velocity_) * (ur - ul));
  }

 private:
  double velocity_;
};

class LocalLaxFriedrichsRiemannSolver {
 public:
  LocalLaxFriedrichsRiemannSolver(double velocity) : velocity_{velocity} {}

  template <typename Derived1, typename Derived2>
  Eigen::VectorXd calc_flux(
      const Eigen::MatrixBase<Derived1>& ul,
      const Eigen::MatrixBase<Derived2>& ur) const noexcept {
    const Eigen::VectorXd a = ul.cwiseAbs().cwiseMax(ur.cwiseAbs());
    return 0.5 * (velocity_ * (ul + ur) - a.cwiseProduct(ur - ul));
  }

 private:
  double velocity_;
};

class HartenRiemannSolver {
 public:
  /**
   * @brief Construct a new Harten Riemann Solver object
   *
   * @param velocity Velocity
   * @param eps Parameter for Entropy fix. (0 <= eps <= 0.5)
   */
  HartenRiemannSolver(double velocity, double eps = 0.25)
      : velocity_{velocity}, eps_{eps} {
    assert(eps >= 0 && eps <= 0.5);
  }

  template <typename Derived1, typename Derived2>
  Eigen::VectorXd calc_flux(
      const Eigen::MatrixBase<Derived1>& ul,
      const Eigen::MatrixBase<Derived2>& ur) const noexcept {
    using Eigen::VectorXd;
    const VectorXd a = 0.5 * (ul + ur).cwiseAbs();
    const VectorXd nu = a.unaryExpr([eps = eps_](double x) {
      if (x < 2 * eps) {
        return 0.25 * x * x / eps + eps;
      } else {
        return x;
      }
    });
    return 0.5 * (velocity_ * (ul + ur) - nu.cwiseProduct(ur - ul));
  }

 private:
  double velocity_;
  double eps_;
};

}  // namespace cfd

#endif  // CFD_RIEMANN_SOLVERS_HPP