#ifndef CFD_SLOPE_LIMITERS_HPP
#define CFD_SLOPE_LIMITERS_HPP

#include <Eigen/Core>

namespace cfd {

/**
 * @brief Minmod limiter
 *
 * @f[
 * \Phi (r) = minmod(1, r) = max(0, min(1, r))
 * @f]
 */
struct MinmodLimiter {
  template <typename Derived>
  static Eigen::VectorXd eval(const Eigen::MatrixBase<Derived>& r) noexcept {
    return r.cwiseMin(1).cwiseMax(0);
  }
};

/**
 * @brief Superbee limiter
 *
 * @f[
 * \Phi (r) = max(0, min(1, 2r), min(2, r))
 * @f]
 */
struct SuperbeeLimiter {
  template <typename Derived>
  static Eigen::VectorXd eval(const Eigen::MatrixBase<Derived>& r) noexcept {
    return (2 * r).cwiseMin(1).cwiseMax(r.cwiseMin(2)).cwiseMax(0);
  }
};

/**
 * @brief van Leer limiter
 *
 * @f[
 * \Phi (r) = \frac{r + |r|}{1 + |r|}
 * @f]
 */
struct VanLeerLimiter {
  template <typename Derived>
  static Eigen::VectorXd eval(const Eigen::MatrixBase<Derived>& r) noexcept {
    return ((r.array() + r.array().abs()) / (1 + r.array().abs())).matrix();
  }
};

/**
 * @brief van Albada limiter
 *
 * @f[
 * \Phi (r) = \frac{r + r^2}{1 + r^2}
 * @f]
 */
struct VanAlbadaLimiter {
  template <typename Derived>
  static Eigen::VectorXd eval(const Eigen::MatrixBase<Derived>& r) noexcept {
    const Eigen::VectorXd r2 = r.array().square().matrix();
    return ((r.array() + r2.array()) / (1 + r2.array())).matrix();
  }
};

}  // namespace cfd

#endif  // CFD_SLOPE_LIMITERS_HPP