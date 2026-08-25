#include "Hamiltonian.hpp"
#include "HarmonicOscillator.hpp"
#include "fixed_interaction.hpp"
#include "jcoupled64_io.hpp"
#include "util.hpp"
#include "vacuum_mscheme_io.hpp"

#include <eigen3/Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using nucleus::Hamiltonian;
using nucleus::HarmonicOscillator;
using nucleus::JBasisSet;
using nucleus::MBasisSet;
using nucleus::Orbit;
using mrimsrg::CoupledHamiltonian;

namespace
{
struct Options
{
    fs::path interaction;
    fs::path flow_output;
    fs::path jcoupled64;
    fs::path output;
    fs::path diagnostic_jcoupled64;
    int Z = -1;
    int N = -1;
    double scalar_tolerance = 1e-9;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_export_no2bpack --interaction FILE "
        "(--flow-output DIR_OR_PAYLOAD | --jcoupled64 FILE) "
        "--output FILE --Z Z --N N [--diagnostic-jcoupled64 FILE] [--scalar-tolerance X]");
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
        else if (key == "--jcoupled64")
            options.jcoupled64 = value;
        else if (key == "--output")
            options.output = value;
        else if (key == "--diagnostic-jcoupled64")
            options.diagnostic_jcoupled64 = value;
        else if (key == "--Z")
            options.Z = std::stoi(value);
        else if (key == "--N")
            options.N = std::stoi(value);
        else if (key == "--scalar-tolerance")
            options.scalar_tolerance = std::stod(value);
        else
            usage_error("unknown option " + key);
    }
    if (options.interaction.empty() || options.output.empty())
        usage_error("--interaction and --output are required");
    if (options.flow_output.empty() == options.jcoupled64.empty())
        usage_error("select exactly one input: --flow-output or --jcoupled64");
    if (options.Z < 0 || options.N < 0)
        usage_error("--Z and --N are required");
    if (!(options.scalar_tolerance > 0.0) || !std::isfinite(options.scalar_tolerance))
        usage_error("--scalar-tolerance must be finite and positive");
    if (!options.diagnostic_jcoupled64.empty() &&
        fs::absolute(options.output).lexically_normal() ==
            fs::absolute(options.diagnostic_jcoupled64).lexically_normal())
        usage_error("--output and --diagnostic-jcoupled64 must name different files");
    return options;
}

std::vector<std::vector<int>> magnetic_substates(const JBasisSet &jbasis,
                                                 const MBasisSet &mbasis)
{
    std::vector<std::vector<int>> result(jbasis.orbit_number());
    for (int p = 0; p < static_cast<int>(mbasis.m_orbit_number()); ++p)
    {
        const int orbit = mbasis.morbit(p).index;
        if (orbit < 0 || orbit >= static_cast<int>(result.size()))
            throw std::runtime_error("m-scheme orbit has an invalid parent j-orbit index");
        result[orbit].push_back(p);
    }
    for (int a = 0; a < static_cast<int>(jbasis.orbit_number()); ++a)
    {
        const std::size_t expected = static_cast<std::size_t>(jbasis.orbit(a).j + 1);
        if (result[a].size() != expected)
            throw std::runtime_error("incomplete magnetic multiplet for j-orbit " + std::to_string(a));
        std::sort(result[a].begin(), result[a].end(), [&](int p, int q) {
            return mbasis.morbit(p).m < mbasis.morbit(q).m;
        });
    }
    return result;
}

double project_two_body(const mrimsrg::VacuumMScheme &dense,
                        const JBasisSet &jbasis,
                        const MBasisSet &mbasis,
                        const std::vector<std::vector<int>> &substates,
                        int a, int b, int c, int d, int J)
{
    const Orbit &oa = jbasis.orbit(a);
    const Orbit &ob = jbasis.orbit(b);
    const Orbit &oc = jbasis.orbit(c);
    const Orbit &od = jbasis.orbit(d);
    double sum = 0.0;
    for (int M2 = -2 * J; M2 <= 2 * J; M2 += 2)
    {
        for (int p : substates[a])
        {
            const Orbit &op = mbasis.morbit(p);
            for (int q : substates[b])
            {
                const Orbit &oq = mbasis.morbit(q);
                if (op.m + oq.m != M2)
                    continue;
                const double cg_bra = util::CGfast(oa.j, ob.j, 2 * J, op.m, oq.m);
                if (cg_bra == 0.0)
                    continue;
                for (int r : substates[c])
                {
                    const Orbit &orr = mbasis.morbit(r);
                    for (int s : substates[d])
                    {
                        const Orbit &os = mbasis.morbit(s);
                        if (orr.m + os.m != M2)
                            continue;
                        const double cg_ket = util::CGfast(oc.j, od.j, 2 * J, orr.m, os.m);
                        if (cg_ket != 0.0)
                            sum += cg_bra * cg_ket * dense.two(p, q, r, s);
                    }
                }
            }
        }
    }
    const double pair_normalization = std::sqrt(a == b ? 2.0 : 1.0) *
                                      std::sqrt(c == d ? 2.0 : 1.0);
    return sum / (pair_normalization * (2 * J + 1));
}

CoupledHamiltonian couple_to_j_scheme(const mrimsrg::VacuumMScheme &dense,
                                      const JBasisSet &jbasis,
                                      const MBasisSet &mbasis)
{
    if (dense.norb != mbasis.m_orbit_number())
        throw std::runtime_error("dense Hamiltonian orbit count does not match the active m-scheme basis");
    util::cg_table_init(jbasis.max_j());
    const auto substates = magnetic_substates(jbasis, mbasis);

    CoupledHamiltonian result;
    result.zero_body = dense.zero_body;
    result.one_body = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(jbasis.orbit_number()),
        static_cast<Eigen::Index>(jbasis.orbit_number()));
    for (int a = 0; a < static_cast<int>(jbasis.orbit_number()); ++a)
    {
        for (int b = 0; b < static_cast<int>(jbasis.orbit_number()); ++b)
        {
            if (!nucleus::check_symmetry_1b_jptz(jbasis.orbit(a), jbasis.orbit(b)))
                continue;
            double sum = 0.0;
            int count = 0;
            for (int p : substates[a])
                for (int q : substates[b])
                    if (mbasis.morbit(p).m == mbasis.morbit(q).m)
                    {
                        sum += dense.one_body(p, q);
                        ++count;
                    }
            if (count == 0)
                throw std::runtime_error("one-body multiplets have no common magnetic substates");
            result.one_body(a, b) = sum / count;
        }
    }
    result.one_body = 0.5 * (result.one_body + result.one_body.transpose()).eval();

    result.two_body.resize(jbasis.j2b_channel_number());
    for (int channel_index = 0;
         channel_index < static_cast<int>(jbasis.j2b_channel_number());
         ++channel_index)
    {
        const auto &channel = jbasis.j2b_channel(channel_index);
        const int size = channel.size();
        result.two_body[channel_index] = Eigen::MatrixXd::Zero(size, size);
        for (int bra = 0; bra < size; ++bra)
        {
            const auto [a, b] = channel[bra];
            for (int ket = bra; ket < size; ++ket)
            {
                const auto [c, d] = channel[ket];
                const double forward = project_two_body(
                    dense, jbasis, mbasis, substates, a, b, c, d, channel.J);
                const double reverse = project_two_body(
                    dense, jbasis, mbasis, substates, c, d, a, b, channel.J);
                const double value = 0.5 * (forward + reverse);
                result.two_body[channel_index](bra, ket) = value;
                result.two_body[channel_index](ket, bra) = value;
            }
        }
    }

    const int norb = static_cast<int>(mbasis.m_orbit_number());
    for (int p = 0; p < norb; ++p)
    {
        const Orbit &op = mbasis.morbit(p);
        for (int q = 0; q < norb; ++q)
        {
            const Orbit &oq = mbasis.morbit(q);
            const double reconstructed = nucleus::check_symmetry_1b_jmptz(op, oq)
                                             ? result.one_body(op.index, oq.index)
                                             : 0.0;
            result.one_body_projection_error = std::max(
                result.one_body_projection_error,
                std::abs(reconstructed - dense.one_body(p, q)));
        }
    }
    for (int p = 0; p < norb; ++p)
        for (int q = 0; q < norb; ++q)
            for (int r = 0; r < norb; ++r)
                for (int s = 0; s < norb; ++s)
                    result.two_body_projection_error = std::max(
                        result.two_body_projection_error,
                        std::abs(mrimsrg::reconstruct_two_body(result.two_body, jbasis, mbasis,
                                                               p, q, r, s) -
                                 dense.two(p, q, r, s)));
    return result;
}

template <typename T> void write_value(std::ofstream &stream, const T &value)
{
    stream.write(reinterpret_cast<const char *>(&value), sizeof(value));
    if (!stream)
        throw std::runtime_error("failed while writing no2bpack payload");
}

void require_float_representable(double value, const char *field)
{
    const double limit = static_cast<double>(std::numeric_limits<float>::max());
    if (!std::isfinite(value) || std::abs(value) > limit)
        throw std::runtime_error(std::string(field) + " cannot be represented by no2bpack float payload");
}

void write_no2bpack(const fs::path &path,
                    const JBasisSet &jbasis,
                    const CoupledHamiltonian &hamiltonian)
{
    if (fs::exists(path))
        throw std::runtime_error("refusing to overwrite existing output: " + path.string());
    if (!path.parent_path().empty())
        fs::create_directories(path.parent_path());

    const int orbit_number = static_cast<int>(jbasis.orbit_number());
    const int obme_number = orbit_number * (orbit_number + 1) / 2;
    std::int64_t tbme_count64 = 0;
    for (int channel_index = 0;
         channel_index < static_cast<int>(jbasis.j2b_channel_number());
         ++channel_index)
    {
        const std::int64_t size = jbasis.j2b_channel(channel_index).size();
        tbme_count64 += size * (size + 1) / 2;
    }
    if (tbme_count64 > std::numeric_limits<int>::max())
        throw std::overflow_error("too many TBMEs for no2bpack header");
    const int tbme_number = static_cast<int>(tbme_count64);

    std::ofstream output(path, std::ios::binary);
    if (!output)
        throw std::runtime_error("cannot create no2bpack output: " + path.string());
    const double hw = jbasis.get_hbar_omega();
    const int emax = jbasis.max_2n_l();
    write_value(output, hw);
    write_value(output, emax);
    write_value(output, orbit_number);
    write_value(output, obme_number);
    write_value(output, tbme_number);
    for (int a = 0; a < orbit_number; ++a)
    {
        const Orbit &orbit = jbasis.orbit(a);
        write_value(output, orbit.n);
        write_value(output, orbit.l);
        write_value(output, orbit.j);
        write_value(output, orbit.tz);
    }
    if (!std::isfinite(hamiltonian.zero_body))
        throw std::runtime_error("zero-body term is non-finite");
    write_value(output, hamiltonian.zero_body);
    for (int a = 0; a < orbit_number; ++a)
        for (int b = a; b < orbit_number; ++b)
        {
            const double value = hamiltonian.one_body(a, b);
            require_float_representable(value, "one-body matrix element");
            write_value(output, static_cast<float>(value));
        }

    HarmonicOscillator ho(hw);
    const float zero_hcm = 0.0F;
    for (int channel_index = 0;
         channel_index < static_cast<int>(jbasis.j2b_channel_number());
         ++channel_index)
    {
        const auto &channel = jbasis.j2b_channel(channel_index);
        const int channel_size = static_cast<int>(channel.size());
        for (int bra = 0; bra < channel_size; ++bra)
        {
            const auto [a, b] = channel[bra];
            for (int ket = bra; ket < channel_size; ++ket)
            {
                const auto [c, d] = channel[ket];
                const double value = hamiltonian.two_body[channel_index](bra, ket);
                require_float_representable(value, "two-body matrix element");
                write_value(output, a);
                write_value(output, b);
                write_value(output, c);
                write_value(output, d);
                write_value(output, channel.J);
                write_value(output, static_cast<float>(value));
                if (!ho.k1k2_is_zero(jbasis.orbit(a), jbasis.orbit(b),
                                     jbasis.orbit(c), jbasis.orbit(d)))
                    write_value(output, zero_hcm);
            }
        }
    }
    output.close();
    if (!output)
        throw std::runtime_error("failed while closing no2bpack output: " + path.string());
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        nucleus::log_init(nucleus::LogLevel::Info);
        const Options options = parse_options(argc, argv);
        // This interaction supplies only the verified orbit/channel table;
        // every matrix element is replaced by the lossless J64 payload.
        mrimsrg::require_validated_interaction_basis(options.interaction);
        Hamiltonian basis_source;
        basis_source.read_minipack(options.interaction.string(), options.Z + options.N, 0.0);
        basis_source.init_mbasis();
        CoupledHamiltonian coupled;
        if (!options.jcoupled64.empty())
        {
            coupled = mrimsrg::read_jcoupled64(
                options.jcoupled64, basis_source.get_jbasis());
            std::cout << "read jcoupled64=" << fs::absolute(options.jcoupled64) << "\n";
        }
        else
        {
            const auto dense = mrimsrg::read_vacuum_mscheme(
                mrimsrg::resolve_vacuum_payload(options.flow_output));
            coupled = couple_to_j_scheme(
                dense, basis_source.get_jbasis(), basis_source.get_mbasis());

            std::cout << std::setprecision(12)
                      << "m-to-J projection: max_one_body_error="
                      << coupled.one_body_projection_error
                      << " max_two_body_error=" << coupled.two_body_projection_error
                      << " tolerance=" << options.scalar_tolerance << "\n";
            if (coupled.one_body_projection_error > options.scalar_tolerance ||
                coupled.two_body_projection_error > options.scalar_tolerance)
                throw std::runtime_error(
                    "m-scheme Hamiltonian is not representable as a scalar J-coupled operator within tolerance");
        }

        if (!options.diagnostic_jcoupled64.empty())
        {
            mrimsrg::write_jcoupled64(
                options.diagnostic_jcoupled64, basis_source.get_jbasis(), coupled);
            std::cout << "wrote diagnostic_jcoupled64="
                      << fs::absolute(options.diagnostic_jcoupled64) << "\n";
        }
        write_no2bpack(options.output, basis_source.get_jbasis(), coupled);
        std::cout << "wrote no2bpack=" << fs::absolute(options.output)
                  << " zero_body=" << coupled.zero_body
                  << " j_orbits=" << basis_source.get_jbasis().orbit_number()
                  << " obmes=" << basis_source.get_jbasis().orbit_number() *
                                         (basis_source.get_jbasis().orbit_number() + 1) / 2
                  << " j_channels=" << basis_source.get_jbasis().j2b_channel_number()
                  << "\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_export_no2bpack: " << error.what() << "\n";
        return 1;
    }
}
