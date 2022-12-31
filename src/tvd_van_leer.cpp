#include <Eigen/Core>

#include "cfd/cfd.hpp"
#include "common.hpp"
#include "text_file_writer.hpp"

int main(int argc, char** argv) {
  using Eigen::VectorXd;
  namespace fs = std::filesystem;

  const auto params = cfd::make_params();
  const VectorXd x = cfd::make_x(params);
  const auto simulator = cfd::make_simulator(
      params, cfd::RoeRiemannSolver{params.velocity},
      cfd::TvdSpacialReconstructor<cfd::VanLeerLimieter>{params},
      cfd::ExplicitEulerScheme{params});

  // Sine wave
  {
    const VectorXd u0 = cfd::make_sine_wave(x);
    const VectorXd uN = simulator.run(u0);
    const auto writer =
        cfd::TextFileWriter{fs::path("result/tvd_van_leer/sine")};
    writer.write(x, "x.txt");
    writer.write(u0, "u0.txt");
    writer.write(uN, "u500.txt");
  }

  // Pulse wave
  {
    const VectorXd u0 = cfd::make_pulse_wave(x);
    const VectorXd uN = simulator.run(u0);
    const auto writer =
        cfd::TextFileWriter{fs::path("result/tvd_van_leer/pulse")};
    writer.write(x, "x.txt");
    writer.write(u0, "u0.txt");
    writer.write(uN, "u500.txt");
  }
}