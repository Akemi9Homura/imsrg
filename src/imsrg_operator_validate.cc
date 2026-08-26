#include "ModelSpace.hh"
#include "Operator.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct Options
{
  int emax = -1;
  double hw = 0.0;
  double tolerance = 1e-10;
  std::string reference;
  std::vector<std::string> files;
};

[[noreturn]] void UsageError(const std::string &message)
{
  throw std::invalid_argument(
      message +
      "\nusage: imsrg_operator_validate --emax E --hw HW --reference NAME "
      "[--tolerance X] OMEGA [OMEGA ...]");
}

Options ParseOptions(int argc, char **argv)
{
  Options options;
  for (int i = 1; i < argc; ++i)
  {
    const std::string key = argv[i];
    if (key == "--emax" || key == "--hw" || key == "--reference" ||
        key == "--tolerance")
    {
      if (i + 1 >= argc)
        UsageError("missing value after " + key);
      const std::string value = argv[++i];
      if (key == "--emax")
        options.emax = std::stoi(value);
      else if (key == "--hw")
        options.hw = std::stod(value);
      else if (key == "--reference")
        options.reference = value;
      else
        options.tolerance = std::stod(value);
    }
    else if (key.compare(0, 2, "--") == 0)
    {
      UsageError("unknown option " + key);
    }
    else
    {
      options.files.push_back(key);
    }
  }
  if (options.emax < 0 || !(options.hw > 0.0) ||
      !(options.tolerance > 0.0) || options.reference.empty() ||
      options.files.empty())
    UsageError("emax, hw, reference, a positive tolerance, and Omega files are required");
  return options;
}

double MaxAbs(const arma::mat &matrix)
{
  return matrix.empty() ? 0.0 : arma::abs(matrix).max();
}
} // namespace

int main(int argc, char **argv)
{
  try
  {
    const Options options = ParseOptions(argc, argv);
    ModelSpace modelspace(options.emax, options.reference, options.reference);
    modelspace.SetHbarOmega(options.hw);

    double maximum_zero_body = 0.0;
    double maximum_one_body_violation = 0.0;
    double maximum_two_body_violation = 0.0;
    bool metadata_passed = true;
    bool finite_passed = true;

    for (const std::string &filename : options.files)
    {
      Operator omega(modelspace);
      std::ifstream stream(filename, std::ios::binary);
      if (!stream)
        throw std::runtime_error("cannot open Omega file: " + filename);
      omega.ReadBinary(stream);
      if (!stream || stream.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("truncated or trailing Omega payload: " + filename);

      metadata_passed = metadata_passed && omega.GetJRank() == 0 &&
                        omega.GetTRank() == 0 && omega.GetParity() == 0 &&
                        omega.GetParticleRank() <= 2 &&
                        omega.IsAntiHermitian() &&
                        omega.TwoBody.IsAntiHermitian();
      finite_passed = finite_passed && std::isfinite(omega.ZeroBody) &&
                      omega.OneBody.is_finite();
      maximum_zero_body = std::max(maximum_zero_body, std::abs(omega.ZeroBody));
      maximum_one_body_violation = std::max(
          maximum_one_body_violation,
          MaxAbs(omega.OneBody + omega.OneBody.t()));

      for (std::size_t channel = 0;
           channel < modelspace.GetNumberTwoBodyChannels(); ++channel)
      {
        const arma::mat &matrix = omega.TwoBody.GetMatrix(channel, channel);
        const std::size_t expected =
            modelspace.GetTwoBodyChannel(static_cast<int>(channel)).GetNumberKets();
        if (matrix.n_rows != expected || matrix.n_cols != expected)
          throw std::runtime_error("Omega two-body channel dimension mismatch: " + filename);
        finite_passed = finite_passed && matrix.is_finite();
        maximum_two_body_violation = std::max(
            maximum_two_body_violation, MaxAbs(matrix + matrix.t()));
      }
    }

    const double maximum_violation = std::max(
        maximum_zero_body,
        std::max(maximum_one_body_violation, maximum_two_body_violation));
    const bool passed = metadata_passed && finite_passed &&
                        maximum_violation <= options.tolerance;
    std::cout << std::boolalpha << std::setprecision(17)
              << "{\"schema\":\"imsrg_scalar_omega_validation_v1\","
              << "\"files\":" << options.files.size() << ","
              << "\"metadata_passed\":" << metadata_passed << ","
              << "\"finite_passed\":" << finite_passed << ","
              << "\"zero_body_max_abs\":" << maximum_zero_body << ","
              << "\"one_body_antihermiticity_max_abs\":"
              << maximum_one_body_violation << ","
              << "\"two_body_antihermiticity_max_abs\":"
              << maximum_two_body_violation << ","
              << "\"tolerance\":" << options.tolerance << ","
              << "\"passed\":" << passed << "}\n";
    return passed ? 0 : 2;
  }
  catch (const std::exception &error)
  {
    std::cerr << "imsrg_operator_validate: " << error.what() << "\n";
    return 1;
  }
}
