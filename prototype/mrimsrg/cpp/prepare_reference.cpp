#include "Hamiltonian.hpp"
#include "fci.hpp"
#include "fci_util.hpp"
#include "fixed_interaction.hpp"
#include "util.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;
using nucleus::CI_truncation;
using nucleus::CIspace;
using nucleus::Hamiltonian;
using nucleus::simpleFCI;

#ifndef SHELL_MODEL_OBS_ROOT
#define SHELL_MODEL_OBS_ROOT "unknown"
#endif
#ifndef SHELL_MODEL_OBS_REVISION
#define SHELL_MODEL_OBS_REVISION "unknown"
#endif

namespace
{

struct Options
{
    fs::path interaction;
    fs::path output;
    int Z = -1;
    int N = -1;
    int nrefmax = 0;
    int max_iter = 400;
    double energy_tolerance = 1e-11;
    double residual_tolerance = 1e-9;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_prepare --interaction FILE --output DIR --Z Z --N N "
        "[--nrefmax N] [--max-iter N]");
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
        else if (key == "--output")
            options.output = value;
        else if (key == "--Z")
            options.Z = std::stoi(value);
        else if (key == "--N")
            options.N = std::stoi(value);
        else if (key == "--nrefmax")
            options.nrefmax = std::stoi(value);
        else if (key == "--max-iter")
            options.max_iter = std::stoi(value);
        else
            usage_error("unknown option " + key);
    }
    if (options.interaction.empty() || options.output.empty() || options.Z < 0 || options.N < 0)
        usage_error("--interaction, --output, --Z and --N are required");
    if (options.Z + options.N <= 0 || options.nrefmax < 0 || options.max_iter <= 0)
        usage_error("invalid particle number, Nrefmax or iteration limit");
    if (!fs::is_regular_file(options.interaction))
        usage_error("interaction file does not exist: " + options.interaction.string());
    return options;
}

template <typename T> const char *npy_descr();
template <> const char *npy_descr<double>() { return "<f8"; }
template <> const char *npy_descr<std::int32_t>() { return "<i4"; }
template <> const char *npy_descr<std::uint8_t>() { return "|u1"; }

std::string shape_string(const std::vector<std::size_t> &shape)
{
    std::ostringstream stream;
    stream << "(";
    for (std::size_t i = 0; i < shape.size(); ++i)
    {
        if (i != 0)
            stream << ", ";
        stream << shape[i];
    }
    if (shape.size() == 1)
        stream << ",";
    stream << ")";
    return stream.str();
}

template <typename T>
void write_npy(const fs::path &path, const std::vector<T> &values, const std::vector<std::size_t> &shape)
{
    static_assert(std::is_trivially_copyable_v<T>);
    std::size_t expected = 1;
    for (const std::size_t extent : shape)
        expected *= extent;
    if (expected != values.size())
        throw std::runtime_error("array size does not match shape for " + path.string());

    std::string header = "{'descr': '" + std::string(npy_descr<T>()) +
                         "', 'fortran_order': False, 'shape': " + shape_string(shape) + ", }";
    constexpr std::size_t preamble_size = 10;
    const std::size_t padding = (16 - (preamble_size + header.size() + 1) % 16) % 16;
    header.append(padding, ' ');
    header.push_back('\n');
    if (header.size() > 65535)
        throw std::runtime_error("NumPy v1 header is too large");

    std::ofstream output(path, std::ios::binary);
    if (!output)
        throw std::runtime_error("cannot create " + path.string());
    const char magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y', '\x01', '\x00'};
    output.write(magic, sizeof(magic));
    const std::uint16_t header_length = static_cast<std::uint16_t>(header.size());
    const char length_bytes[2] = {static_cast<char>(header_length & 0xff),
                                  static_cast<char>((header_length >> 8) & 0xff)};
    output.write(length_bytes, sizeof(length_bytes));
    output.write(header.data(), static_cast<std::streamsize>(header.size()));
    output.write(reinterpret_cast<const char *>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(T)));
    if (!output)
        throw std::runtime_error("failed while writing " + path.string());
}

std::string json_escape(const std::string &value)
{
    std::ostringstream stream;
    for (const char character : value)
    {
        if (character == '\\' || character == '"')
            stream << '\\';
        stream << character;
    }
    return stream.str();
}

void require_new_output_directory(const fs::path &path)
{
    if (fs::exists(path))
        throw std::runtime_error("refusing to overwrite existing output directory: " + path.string());
    if (!fs::create_directories(path))
        throw std::runtime_error("cannot create output directory: " + path.string());
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        nucleus::log_init(nucleus::LogLevel::Info);
        const Options options = parse_options(argc, argv);
        const int A = options.Z + options.N;
        const std::string interaction_sha256 = mrimsrg::require_fixed_interaction(options.interaction);

        Hamiltonian hamiltonian;
        hamiltonian.read_minipack(options.interaction.string(), A, 0.0);
        const int emax = hamiltonian.get_jbasis().max_2n_l();
        if (std::abs(hamiltonian.get_jbasis().get_hbar_omega() - 20.0) > 1e-12 || emax != 2)
            throw std::runtime_error("the rapid prototype requires hw=20 MeV and emax=2");
        hamiltonian.init_mscheme();

        const int hw_min = hamiltonian.get_jbasis().min_2n_l(options.Z, options.N);
        CI_truncation truncation;
        truncation.hw_trunc.hw = hw_min + options.nrefmax;
        truncation.parity_trunc.parity = 0;
        CIspace space(&hamiltonian, options.Z, options.N, truncation);
        simpleFCI solver(space, 1);
        if (solver.get_configs().size() == 1)
            solver.brute_force();
        else
            solver.lanczos_otf(options.max_iter, options.energy_tolerance,
                               options.residual_tolerance, nucleus::Lanczos::default_seed, false);
        solver.compute_J_otf();
        const std::vector<int> state_j2 = solver.get_states_J2();
        if (state_j2.size() != 1 || state_j2.front() != 0)
            throw std::runtime_error("lowest positive-parity NCSM state is not J=0");

        const auto &mbasis = hamiltonian.get_mbasis();
        const std::size_t norb = mbasis.m_orbit_number();
        const auto &configs = solver.get_configs();
        const auto coefficients = solver.get_eigenvector(0);

        std::vector<std::int32_t> orbits;
        orbits.reserve(norb * 6);
        for (std::size_t p = 0; p < norb; ++p)
        {
            const auto &orbit = mbasis.morbit(static_cast<int>(p));
            orbits.insert(orbits.end(), {orbit.index, orbit.n, orbit.l, orbit.j, orbit.m, orbit.tz});
        }

        std::vector<double> one_body(norb * norb);
        for (std::size_t p = 0; p < norb; ++p)
            for (std::size_t q = 0; q < norb; ++q)
                one_body[p * norb + q] = hamiltonian.get_m1b(static_cast<int>(p), static_cast<int>(q));

        std::vector<double> two_body(norb * norb * norb * norb);
        for (std::size_t p = 0; p < norb; ++p)
            for (std::size_t q = 0; q < norb; ++q)
                for (std::size_t r = 0; r < norb; ++r)
                    for (std::size_t s = 0; s < norb; ++s)
                    {
                        const std::size_t index = ((p * norb + q) * norb + r) * norb + s;
                        two_body[index] = hamiltonian.get_mvmat(static_cast<int>(p), static_cast<int>(q),
                                                               static_cast<int>(r), static_cast<int>(s));
                    }

        std::vector<std::uint8_t> determinants(configs.size() * norb);
        for (std::size_t i = 0; i < configs.size(); ++i)
            for (std::size_t p = 0; p < norb; ++p)
                determinants[i * norb + p] = configs[i].test(p) ? 1 : 0;
        std::vector<double> coefficient_values(coefficients.data(), coefficients.data() + coefficients.size());

        require_new_output_directory(options.output);
        write_npy(options.output / "orbits.npy", orbits, {norb, 6});
        write_npy(options.output / "one_body.npy", one_body, {norb, norb});
        write_npy(options.output / "two_body.npy", two_body, {norb, norb, norb, norb});
        write_npy(options.output / "determinants.npy", determinants, {configs.size(), norb});
        write_npy(options.output / "coefficients.npy", coefficient_values, {configs.size()});

        std::ofstream metadata(options.output / "metadata.json");
        metadata << std::setprecision(17)
                 << "{\n"
                 << "  \"schema\": \"mrimsrg_reference_v1\",\n"
                 << "  \"interaction\": \"" << json_escape(fs::absolute(options.interaction).string()) << "\",\n"
                 << "  \"interaction_sha256\": \"" << interaction_sha256 << "\",\n"
                 << "  \"shell_model_obs_root\": \"" << json_escape(SHELL_MODEL_OBS_ROOT) << "\",\n"
                 << "  \"shell_model_obs_revision\": \"" << SHELL_MODEL_OBS_REVISION << "\",\n"
                 << "  \"one_body_convention\": \"t[p,q] a^dagger_p a_q\",\n"
                 << "  \"two_body_convention\": \"(1/4) V[p,q,r,s] a^dagger_p a^dagger_q a_s a_r\",\n"
                 << "  \"A\": " << A << ",\n"
                 << "  \"Z\": " << options.Z << ",\n"
                 << "  \"N\": " << options.N << ",\n"
                 << "  \"hw\": " << hamiltonian.get_jbasis().get_hbar_omega() << ",\n"
                 << "  \"emax\": " << emax << ",\n"
                 << "  \"e2max\": 4,\n"
                 << "  \"Nrefmax\": " << options.nrefmax << ",\n"
                 << "  \"J2\": " << state_j2.front() << ",\n"
                 << "  \"parity\": 1,\n"
                 << "  \"zero_body\": " << hamiltonian.get_zero_body() << ",\n"
                 << "  \"reference_energy\": " << solver.get_eigenvalue(0) << ",\n"
                 << "  \"reference_dimension\": " << configs.size() << ",\n"
                 << "  \"lanczos_iterations\": " << solver.get_last_lanczos_iterations() << "\n"
                 << "}\n";
        if (!metadata)
            throw std::runtime_error("failed while writing metadata.json");

        std::cout << std::setprecision(12)
                  << "prepared A=" << A << " Z=" << options.Z << " N=" << options.N
                  << " Nrefmax=" << options.nrefmax << " dimension=" << configs.size()
                  << " Eref=" << solver.get_eigenvalue(0) << " in " << options.output << "\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_prepare: " << error.what() << "\n";
        return 1;
    }
}
