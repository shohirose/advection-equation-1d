#include <gtest/gtest.h>

#include "cfd/cfd.hpp"

using Eigen::VectorXd;

TEST(CfdTest, PeriodicBoundaryTest) {
  VectorXd u = VectorXd::LinSpaced(11, 0.0, 10.0);
  const cfd::PeriodicBoundary boundary{2, 7};

  boundary.apply(u);

  EXPECT_DOUBLE_EQ(u(0), 7.0);
  EXPECT_DOUBLE_EQ(u(1), 8.0);
  EXPECT_DOUBLE_EQ(u(9), 2.0);
  EXPECT_DOUBLE_EQ(u(10), 3.0);
}

TEST(CfdTest, RoeRiemannSolverTest) {
  cfd::RoeRiemannSolver solver{1.0};
  VectorXd ul(3);
  ul << 0.0, 1.0, 2.0;
  VectorXd ur(3);
  ur << 1.0, 2.0, 3.0;

  const VectorXd f = solver.calc_flux(ul, ur);
  ASSERT_EQ(f.size(), 3);

  EXPECT_DOUBLE_EQ(f(0), 0.0);
  EXPECT_DOUBLE_EQ(f(1), 1.0);
  EXPECT_DOUBLE_EQ(f(2), 2.0);
}

TEST(CfdTest, FirstOrderSpacialReconstuctorTest) {
  cfd::FirstOrderSpacialReconstructor reconstructor{2, 3};
  VectorXd u = VectorXd::LinSpaced(7, 0.0, 6.0);
  const VectorXd ul = reconstructor.calc_left(u);
  const VectorXd ur = reconstructor.calc_right(u);

  ASSERT_EQ(ul.size(), 4);
  ASSERT_EQ(ur.size(), 4);
  EXPECT_DOUBLE_EQ(ul(0), u(1));
  EXPECT_DOUBLE_EQ(ul(1), u(2));
  EXPECT_DOUBLE_EQ(ul(2), u(3));
  EXPECT_DOUBLE_EQ(ul(3), u(4));
  EXPECT_DOUBLE_EQ(ur(0), u(2));
  EXPECT_DOUBLE_EQ(ur(1), u(3));
  EXPECT_DOUBLE_EQ(ur(2), u(4));
  EXPECT_DOUBLE_EQ(ur(3), u(5));
}

TEST(CfdTest, ExplicitEulerSchemeTest) {
  cfd::ExplicitEulerScheme euler{0.1, 0.1, 2, 3};
  VectorXd u = VectorXd::Zero(7);
  VectorXd f = VectorXd::LinSpaced(4, 0.0, 3.0);
  euler.update(u, f);

  ASSERT_EQ(u.size(), 7);
  EXPECT_DOUBLE_EQ(u(0), 0.0);
  EXPECT_DOUBLE_EQ(u(1), 0.0);
  EXPECT_DOUBLE_EQ(u(2), -1.0);
  EXPECT_DOUBLE_EQ(u(3), -1.0);
  EXPECT_DOUBLE_EQ(u(4), -1.0);
  EXPECT_DOUBLE_EQ(u(5), 0.0);
  EXPECT_DOUBLE_EQ(u(6), 0.0);
}