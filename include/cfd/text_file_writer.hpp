#ifndef CFD_TEXT_FILE_WRITER_HPP
#define CFD_TEXT_FILE_WRITER_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace cfd {

class TextFileWriter {
 public:
  /**
   * @brief Construct a new Text File Writer object
   *
   * @param directory Directory to output files
   */
  TextFileWriter(const std::filesystem::path& directory)
      : directory_{directory} {}

  /**
   * @brief Write data to a file.
   *
   * @param x Data
   * @param filename File name
   */
  template <typename Derived>
  void write(const Eigen::MatrixBase<Derived>& x,
             const std::string& filename) const noexcept {
    namespace fs = std::filesystem;
    if (!fs::exists(directory_)) {
      std::error_code ec;
      fs::create_directories(directory_, ec);
      if (ec) {
        fmt::print(stderr, "Failed to create a directory: {}\n",
                   directory_.string());
        fmt::print(stderr, "Error code: {}\n", ec.message());
        std::exit(EXIT_FAILURE);
      }
    }
    std::ofstream file(directory_ / fs::path(filename));
    file << x << std::endl;
  }

 private:
  std::filesystem::path directory_;  ///> Directory to output files
};

}  // namespace cfd

#endif  // CFD_TEXT_FILE_WRITER_HPP