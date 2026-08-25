#include "Hamiltonian.hpp"
#include "fci.hpp"
#include "fci_util.hpp"
#include "fixed_interaction.hpp"
#include "jcoupled64_io.hpp"
#include "util.hpp"
#include "vacuum_mscheme_io.hpp"

#include <filesystem>
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
    fs::path no2bpack;
    fs::path jcoupled64;
    int Z = -1;
    int N = -1;
    int nmax = 0;
    int max_iter = 500;
    int states = 3;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_validate (--interaction FILE --flow-output DIR | --no2bpack FILE | "
        "--interaction FILE --jcoupled64 FILE) "
        "--Z Z --N N "
        "[--nmax N] [--max-iter N] [--states N]");
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
        else if (key == "--no2bpack")
            options.no2bpack = value;
        else if (key == "--jcoupled64")
            options.jcoupled64 = value;
        else if (key == "--Z")
            options.Z = std::stoi(value);
        else if (key == "--N")
            options.N = std::stoi(value);
        else if (key == "--nmax")
            options.nmax = std::stoi(value);
        else if (key == "--max-iter")
            options.max_iter = std::stoi(value);
        else if (key == "--states")
            options.states = std::stoi(value);
        else
            usage_error("unknown option " + key);
    }
    const bool dense_input = !options.flow_output.empty();
    const bool packed_input = !options.no2bpack.empty();
    const bool jcoupled_input = !options.jcoupled64.empty();
    if (static_cast<int>(dense_input) + static_cast<int>(packed_input) +
            static_cast<int>(jcoupled_input) !=
        1)
        usage_error("select exactly one input: --flow-output, --no2bpack, or --jcoupled64");
    if (dense_input && options.interaction.empty())
        usage_error("dense input requires both --interaction and --flow-output");
    if (jcoupled_input && options.interaction.empty())
        usage_error("jcoupled64 input requires --interaction to define the active basis");
    if (packed_input && !options.interaction.empty())
        usage_error("standard no2bpack input does not use --interaction");
    if (options.Z < 0 || options.N < 0)
        usage_error("--Z and --N are required");
    if (options.nmax < 0 || options.max_iter <= 0 || options.states <= 0)
        usage_error("invalid Nmax, iteration limit, or state count");
    return options;
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
        if (!options.no2bpack.empty())
        {
            hamiltonian.read_no2bpack(options.no2bpack.string(), A, 0.0);
            hamiltonian.init_mscheme();
        }
        else
        {
            mrimsrg::require_fixed_interaction(options.interaction);
            hamiltonian.read_minipack(options.interaction.string(), A, 0.0);
            hamiltonian.init_mscheme();
            const auto dense = options.jcoupled64.empty()
                                   ? mrimsrg::read_vacuum_mscheme(
                                         mrimsrg::resolve_vacuum_payload(options.flow_output))
                                   : mrimsrg::reconstruct_vacuum_mscheme(
                                         mrimsrg::read_jcoupled64(
                                             options.jcoupled64, hamiltonian.get_jbasis()),
                                         hamiltonian.get_jbasis(), hamiltonian.get_mbasis());
            const std::size_t norb = hamiltonian.get_mbasis().m_orbit_number();
            if (dense.norb != norb)
                throw std::runtime_error("bridge payload orbit count does not match interaction basis");
            hamiltonian.replace_mscheme(dense.zero_body, dense.one_body, dense.two_body);
        }

        CI_truncation truncation;
        truncation.hw_trunc.hw = hamiltonian.get_jbasis().min_2n_l(options.Z, options.N) + options.nmax;
        truncation.parity_trunc.parity = 0;
        CIspace space(&hamiltonian, options.Z, options.N, truncation);
        simpleFCI solver(space, options.states);
        if (solver.get_configs().size() == 1)
            solver.brute_force();
        else
            solver.lanczos_otf(options.max_iter, 1e-11, 1e-9, nucleus::Lanczos::default_seed, false);
        solver.compute_J_otf();
        const auto state_j2 = solver.get_states_J2();
        std::cout << std::setprecision(17) << "readback A=" << A << " Nmax=" << options.nmax
                  << " dimension=" << solver.get_configs().size()
                  << " states=" << state_j2.size() << " E0=" << solver.get_eigenvalue(0)
                  << " twoJ=" << state_j2.at(0) << "\n";
        for (std::size_t state = 0; state < state_j2.size(); ++state)
        {
            const double energy = solver.get_eigenvalue(static_cast<int>(state));
            std::cout << "state=" << state << " E=" << energy
                      << " Ex=" << energy - solver.get_eigenvalue(0)
                      << " twoJ=" << state_j2[state] << "\n";
        }
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_validate: " << error.what() << "\n";
        return 1;
    }
}
