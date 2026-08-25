#include "Hamiltonian.hpp"
#include "jcoupled64_io.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
using nucleus::Hamiltonian;
using mrimsrg::CoupledHamiltonian;

namespace
{
struct Options
{
    fs::path interaction;
    fs::path output;
    int mass_number = 0;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_minipack_to_j64 --interaction FILE --output FILE --A MASS");
}

Options parse_options(int argc, char **argv)
{
    Options options;
    for (int i = 1; i < argc; i += 2)
    {
        if (i + 1 >= argc)
            usage_error("missing value for " + std::string(argv[i]));
        const std::string key = argv[i];
        const std::string value = argv[i + 1];
        if (key == "--interaction")
            options.interaction = value;
        else if (key == "--output")
            options.output = value;
        else if (key == "--A")
            options.mass_number = std::stoi(value);
        else
            usage_error("unknown option " + key);
    }
    if (options.interaction.empty() || options.output.empty() || options.mass_number <= 0)
        usage_error("--interaction, --output, and positive --A are required");
    if (!fs::is_regular_file(options.interaction))
        usage_error("interaction is not a regular file: " + options.interaction.string());
    return options;
}

CoupledHamiltonian copy_j_scheme(const Hamiltonian &hamiltonian)
{
    const auto &basis = hamiltonian.get_jbasis();
    CoupledHamiltonian result;
    result.zero_body = hamiltonian.get_zero_body();
    const int norbits = basis.orbit_number();
    result.one_body = Eigen::MatrixXd::Zero(norbits, norbits);
    for (int a = 0; a < norbits; ++a)
        for (int b = 0; b < norbits; ++b)
            result.one_body(a, b) = hamiltonian.get_j1b(a, b);
    result.two_body.reserve(basis.j2b_channel_number());
    for (std::size_t channel = 0; channel < basis.j2b_channel_number(); ++channel)
        result.two_body.push_back(hamiltonian.get_j2b_vmat_of_channel(channel));
    return result;
}

double maximum_difference(const CoupledHamiltonian &left,
                          const CoupledHamiltonian &right)
{
    double maximum = std::abs(left.zero_body - right.zero_body);
    if (left.one_body.rows() != right.one_body.rows() ||
        left.one_body.cols() != right.one_body.cols() ||
        left.two_body.size() != right.two_body.size())
        throw std::runtime_error("jcoupled64 roundtrip changed matrix dimensions");
    if (left.one_body.size() > 0)
        maximum = std::max(
            maximum, (left.one_body - right.one_body).cwiseAbs().maxCoeff());
    for (std::size_t channel = 0; channel < left.two_body.size(); ++channel)
    {
        if (left.two_body[channel].rows() != right.two_body[channel].rows() ||
            left.two_body[channel].cols() != right.two_body[channel].cols())
            throw std::runtime_error("jcoupled64 roundtrip changed channel dimensions");
        if (left.two_body[channel].size() > 0)
            maximum = std::max(
                maximum,
                (left.two_body[channel] - right.two_body[channel])
                    .cwiseAbs()
                    .maxCoeff());
    }
    return maximum;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        nucleus::log_init(nucleus::LogLevel::Info);
        const Options options = parse_options(argc, argv);
        Hamiltonian hamiltonian;
        // Reuse shell-model-obs' validated A-dependent intrinsic kinetic
        // energy. No m-scheme Hamiltonian is initialized in this bridge.
        hamiltonian.read_minipack(
            options.interaction.string(), options.mass_number, 0.0);
        const CoupledHamiltonian coupled = copy_j_scheme(hamiltonian);
        mrimsrg::write_jcoupled64(
            options.output, hamiltonian.get_jbasis(), coupled);
        const CoupledHamiltonian roundtrip = mrimsrg::read_jcoupled64(
            options.output, hamiltonian.get_jbasis());
        const double error = maximum_difference(coupled, roundtrip);
        if (error != 0.0)
            throw std::runtime_error(
                "lossless minipack-to-jcoupled64 roundtrip is not bitwise exact");
        std::cout << std::setprecision(17)
                  << "interaction=" << fs::absolute(options.interaction) << "\n"
                  << "output=" << fs::absolute(options.output) << "\n"
                  << "A=" << options.mass_number << "\n"
                  << "hw=" << hamiltonian.get_jbasis().get_hbar_omega() << "\n"
                  << "emax=" << hamiltonian.get_jbasis().max_2n_l() << "\n"
                  << "j_orbits=" << hamiltonian.get_jbasis().orbit_number() << "\n"
                  << "j_channels=" << hamiltonian.get_jbasis().j2b_channel_number() << "\n"
                  << "roundtrip_max_abs_mev=" << error << "\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_minipack_to_j64: " << error.what() << "\n";
        return 1;
    }
}
