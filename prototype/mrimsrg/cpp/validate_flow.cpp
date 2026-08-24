#include "Hamiltonian.hpp"
#include "fci.hpp"
#include "fci_util.hpp"
#include "util.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using nucleus::CI_truncation;
using nucleus::CIspace;
using nucleus::Hamiltonian;
using nucleus::simpleFCI;

namespace
{
struct Options
{
    fs::path interaction;
    fs::path flow_output;
    int Z = -1;
    int N = -1;
    int nmax = 0;
    int max_iter = 500;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_validate --interaction FILE --flow-output DIR --Z Z --N N "
        "[--nmax N] [--max-iter N]");
}

Options parse_options(int argc, char **argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string key = argv[i];
        if (i + 1 >= argc)
            usage_error("missing value after " + key);
        const std::string value = argv[++i];
        if (key == "--interaction")
            options.interaction = value;
        else if (key == "--flow-output")
            options.flow_output = value;
        else if (key == "--Z")
            options.Z = std::stoi(value);
        else if (key == "--N")
            options.N = std::stoi(value);
        else if (key == "--nmax")
            options.nmax = std::stoi(value);
        else if (key == "--max-iter")
            options.max_iter = std::stoi(value);
        else
            usage_error("unknown option " + key);
    }
    if (options.interaction.empty() || options.flow_output.empty() || options.Z < 0 || options.N < 0)
        usage_error("--interaction, --flow-output, --Z and --N are required");
    if (options.nmax < 0 || options.max_iter <= 0)
        usage_error("invalid Nmax or iteration limit");
    return options;
}

template <typename T> void read_exact(std::ifstream &stream, T *destination, std::size_t count, const char *field)
{
    stream.read(reinterpret_cast<char *>(destination), static_cast<std::streamsize>(count * sizeof(T)));
    if (!stream)
        throw std::runtime_error(std::string("truncated bridge payload while reading ") + field);
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        nucleus::log_init(nucleus::LogLevel::Info);
        const Options options = parse_options(argc, argv);
        const int A = options.Z + options.N;
        Hamiltonian hamiltonian;
        hamiltonian.read_minipack(options.interaction.string(), A, 0.0);
        hamiltonian.init_mscheme();

        const fs::path payload_path = options.flow_output / "vacuum_mscheme.bin";
        std::ifstream payload(payload_path, std::ios::binary);
        if (!payload)
            throw std::runtime_error("cannot open bridge payload: " + payload_path.string());
        std::array<char, 16> magic{};
        read_exact(payload, magic.data(), magic.size(), "magic");
        const std::array<char, 16> expected = {'m', 'r', 'i', 'm', 's', 'r', 'g', '_',
                                               'm', '_', 'v', '1', 0,   0,   0,   0};
        if (magic != expected)
            throw std::runtime_error("unsupported MR-IMSRG bridge payload");
        std::uint64_t norb_file = 0;
        double zero_body = 0.0;
        read_exact(payload, &norb_file, 1, "norb");
        read_exact(payload, &zero_body, 1, "zero body");
        const std::size_t norb = hamiltonian.get_mbasis().m_orbit_number();
        if (norb_file != norb)
            throw std::runtime_error("bridge payload orbit count does not match interaction basis");
        Eigen::MatrixXd one_body(norb, norb);
        for (std::size_t p = 0; p < norb; ++p)
            for (std::size_t q = 0; q < norb; ++q)
                read_exact(payload, &one_body(p, q), 1, "one body");
        std::vector<double> two_body(norb * norb * norb * norb);
        read_exact(payload, two_body.data(), two_body.size(), "two body");
        if (payload.peek() != std::ifstream::traits_type::eof())
            throw std::runtime_error("bridge payload contains trailing data");
        hamiltonian.replace_mscheme(zero_body, one_body, two_body);

        CI_truncation truncation;
        truncation.hw_trunc.hw = hamiltonian.get_jbasis().min_2n_l(options.Z, options.N) + options.nmax;
        truncation.parity_trunc.parity = 0;
        CIspace space(&hamiltonian, options.Z, options.N, truncation);
        simpleFCI solver(space, 1);
        if (solver.get_configs().size() == 1)
            solver.brute_force();
        else
            solver.lanczos_otf(options.max_iter, 1e-11, 1e-9, nucleus::Lanczos::default_seed, false);
        solver.compute_J_otf();
        const auto state_j2 = solver.get_states_J2();
        std::cout << std::setprecision(12) << "readback A=" << A << " Nmax=" << options.nmax
                  << " dimension=" << solver.get_configs().size()
                  << " E0=" << solver.get_eigenvalue(0) << " J2=" << state_j2.at(0) << "\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_validate: " << error.what() << "\n";
        return 1;
    }
}
