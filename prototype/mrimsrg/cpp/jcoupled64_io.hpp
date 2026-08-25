#pragma once

#include "Hamiltonian.hpp"
#include "util.hpp"
#include "vacuum_mscheme_io.hpp"

#include <eigen3/Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mrimsrg
{
namespace fs = std::filesystem;

// Double-precision, J-coupled diagnostic representation. This is deliberately
// distinct from production no2bpack, whose OBME/TBME payload is float32.
struct CoupledHamiltonian
{
    double zero_body = 0.0;
    Eigen::MatrixXd one_body;
    std::vector<Eigen::MatrixXd> two_body;
    double one_body_projection_error = 0.0;
    double two_body_projection_error = 0.0;
};

inline double reconstruct_two_body(const std::vector<Eigen::MatrixXd> &coupled,
                                   const nucleus::JBasisSet &jbasis,
                                   const nucleus::MBasisSet &mbasis,
                                   int p, int q, int r, int s)
{
    if (p == q || r == s)
        return 0.0;
    const nucleus::Orbit &op = mbasis.morbit(p);
    const nucleus::Orbit &oq = mbasis.morbit(q);
    const nucleus::Orbit &orr = mbasis.morbit(r);
    const nucleus::Orbit &os = mbasis.morbit(s);
    if (!nucleus::check_symmetry_2b_mptz(op, oq, orr, os))
        return 0.0;

    const int ja = op.j;
    const int jb = oq.j;
    const int jc = orr.j;
    const int jd = os.j;
    const int jmin = std::max(std::abs(ja - jb), std::abs(jc - jd)) / 2;
    const int jmax = std::min(ja + jb, jc + jd) / 2;
    int a = op.index;
    int b = oq.index;
    int c = orr.index;
    int d = os.index;
    const bool swap_ab = a > b;
    const bool swap_cd = c > d;
    if (swap_ab)
        std::swap(a, b);
    if (swap_cd)
        std::swap(c, d);

    double value = 0.0;
    for (int J = jmin; J <= jmax; ++J)
    {
        if (a == b && nucleus::isodd(J))
            continue;
        if (c == d && nucleus::isodd(J))
            continue;
        double phase = 1.0;
        if (a == b)
            phase *= nucleus::sqrt_2;
        if (c == d)
            phase *= nucleus::sqrt_2;
        if (swap_ab)
            phase *= nucleus::iphase((ja + jb) / 2 - J + 1);
        if (swap_cd)
            phase *= nucleus::iphase((jc + jd) / 2 - J + 1);
        const auto [channel, bra, ket] = jbasis.jvmat_index(a, b, c, d, J);
        value += phase * coupled[channel](bra, ket) *
                 util::CGfast(ja, jb, 2 * J, op.m, oq.m) *
                 util::CGfast(jc, jd, 2 * J, orr.m, os.m);
    }
    return value;
}

inline VacuumMScheme reconstruct_vacuum_mscheme(const CoupledHamiltonian &coupled,
                                                const nucleus::JBasisSet &jbasis,
                                                const nucleus::MBasisSet &mbasis)
{
    util::cg_table_init(jbasis.max_j());
    VacuumMScheme dense;
    dense.norb = mbasis.m_orbit_number();
    dense.zero_body = coupled.zero_body;
    dense.one_body = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(dense.norb), static_cast<Eigen::Index>(dense.norb));
    for (int p = 0; p < static_cast<int>(dense.norb); ++p)
        for (int q = 0; q < static_cast<int>(dense.norb); ++q)
        {
            const nucleus::Orbit &op = mbasis.morbit(p);
            const nucleus::Orbit &oq = mbasis.morbit(q);
            if (nucleus::check_symmetry_1b_jmptz(op, oq))
                dense.one_body(p, q) = coupled.one_body(op.index, oq.index);
        }

    const std::size_t n = static_cast<std::size_t>(dense.norb);
    if (n != 0 && n > std::numeric_limits<std::size_t>::max() / n)
        throw std::overflow_error("jcoupled64 orbit count overflows tensor dimensions");
    const std::size_t n2 = n * n;
    if (n2 != 0 && n2 > std::numeric_limits<std::size_t>::max() / n2)
        throw std::overflow_error("jcoupled64 tensor size overflows address space");
    dense.two_body.resize(n2 * n2);
    for (int p = 0; p < static_cast<int>(dense.norb); ++p)
        for (int q = 0; q < static_cast<int>(dense.norb); ++q)
            for (int r = 0; r < static_cast<int>(dense.norb); ++r)
                for (int s = 0; s < static_cast<int>(dense.norb); ++s)
                    dense.two_body[((static_cast<std::size_t>(p) * n + q) * n + r) * n + s] =
                        reconstruct_two_body(coupled.two_body, jbasis, mbasis, p, q, r, s);
    return dense;
}

template <typename T>
inline void write_j64_value(std::ofstream &stream, const T &value, const char *field)
{
    stream.write(reinterpret_cast<const char *>(&value), sizeof(value));
    if (!stream)
        throw std::runtime_error(std::string("failed while writing jcoupled64 ") + field);
}

template <typename T>
inline void read_j64_value(std::ifstream &stream, T &value, const char *field)
{
    stream.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!stream)
        throw std::runtime_error(std::string("truncated jcoupled64 payload while reading ") + field);
}

inline std::uint64_t jcoupled_tbme_count(const nucleus::JBasisSet &jbasis)
{
    std::uint64_t count = 0;
    for (int channel = 0; channel < static_cast<int>(jbasis.j2b_channel_number()); ++channel)
    {
        const std::uint64_t size = jbasis.j2b_channel(channel).size();
        count += size * (size + 1) / 2;
    }
    return count;
}

inline void write_jcoupled64(const fs::path &path,
                             const nucleus::JBasisSet &jbasis,
                             const CoupledHamiltonian &hamiltonian)
{
    if (fs::exists(path))
        throw std::runtime_error("refusing to overwrite existing output: " + path.string());
    if (!path.parent_path().empty())
        fs::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary);
    if (!output)
        throw std::runtime_error("cannot create jcoupled64 output: " + path.string());

    const std::array<char, 16> magic = {'m', 'r', 'i', 'm', 's', 'r', 'g', '_',
                                        'j', '6', '4', '_', 'v', '1', 0,   0};
    output.write(magic.data(), static_cast<std::streamsize>(magic.size()));
    const double hw = jbasis.get_hbar_omega();
    const std::int32_t emax = jbasis.max_2n_l();
    const std::uint64_t orbit_count = jbasis.orbit_number();
    const std::uint64_t obme_count = orbit_count * (orbit_count + 1) / 2;
    const std::uint64_t tbme_count = jcoupled_tbme_count(jbasis);
    write_j64_value(output, hw, "hw");
    write_j64_value(output, emax, "emax");
    write_j64_value(output, orbit_count, "orbit count");
    write_j64_value(output, obme_count, "OBME count");
    write_j64_value(output, tbme_count, "TBME count");
    for (int a = 0; a < static_cast<int>(orbit_count); ++a)
    {
        const nucleus::Orbit &orbit = jbasis.orbit(a);
        for (int quantum_number : {orbit.n, orbit.l, orbit.j, orbit.tz})
        {
            const std::int32_t value = quantum_number;
            write_j64_value(output, value, "orbit table");
        }
    }
    write_j64_value(output, hamiltonian.zero_body, "zero body");
    for (int a = 0; a < static_cast<int>(orbit_count); ++a)
        for (int b = a; b < static_cast<int>(orbit_count); ++b)
            write_j64_value(output, hamiltonian.one_body(a, b), "OBME");
    for (int channel_index = 0;
         channel_index < static_cast<int>(jbasis.j2b_channel_number()); ++channel_index)
    {
        const auto &channel = jbasis.j2b_channel(channel_index);
        for (int bra = 0; bra < static_cast<int>(channel.size()); ++bra)
        {
            const auto [a, b] = channel[bra];
            for (int ket = bra; ket < static_cast<int>(channel.size()); ++ket)
            {
                const auto [c, d] = channel[ket];
                for (int index : {a, b, c, d, channel.J})
                {
                    const std::int32_t value = index;
                    write_j64_value(output, value, "TBME index");
                }
                write_j64_value(output, hamiltonian.two_body[channel_index](bra, ket), "TBME");
            }
        }
    }
    output.close();
    if (!output)
        throw std::runtime_error("failed while closing jcoupled64 output: " + path.string());
}

inline CoupledHamiltonian read_jcoupled64(const fs::path &path,
                                         const nucleus::JBasisSet &jbasis)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
        throw std::runtime_error("cannot open jcoupled64 input: " + path.string());
    std::array<char, 16> magic{};
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    const std::array<char, 16> expected = {'m', 'r', 'i', 'm', 's', 'r', 'g', '_',
                                           'j', '6', '4', '_', 'v', '1', 0,   0};
    if (!input || magic != expected)
        throw std::runtime_error("unsupported jcoupled64 payload");

    double hw = 0.0;
    std::int32_t emax = -1;
    std::uint64_t orbit_count = 0;
    std::uint64_t obme_count = 0;
    std::uint64_t tbme_count = 0;
    read_j64_value(input, hw, "hw");
    read_j64_value(input, emax, "emax");
    read_j64_value(input, orbit_count, "orbit count");
    read_j64_value(input, obme_count, "OBME count");
    read_j64_value(input, tbme_count, "TBME count");
    if (orbit_count != jbasis.orbit_number() || emax != jbasis.max_2n_l() ||
        std::abs(hw - jbasis.get_hbar_omega()) > 1e-12 ||
        obme_count != orbit_count * (orbit_count + 1) / 2 ||
        tbme_count != jcoupled_tbme_count(jbasis))
        throw std::runtime_error("jcoupled64 header does not match the active interaction basis");
    for (int a = 0; a < static_cast<int>(orbit_count); ++a)
    {
        const nucleus::Orbit &orbit = jbasis.orbit(a);
        const std::array<int, 4> expected_qn = {orbit.n, orbit.l, orbit.j, orbit.tz};
        for (int expected_value : expected_qn)
        {
            std::int32_t value = 0;
            read_j64_value(input, value, "orbit table");
            if (value != expected_value)
                throw std::runtime_error("jcoupled64 orbit table does not match the active basis");
        }
    }

    CoupledHamiltonian result;
    read_j64_value(input, result.zero_body, "zero body");
    result.one_body = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(orbit_count), static_cast<Eigen::Index>(orbit_count));
    for (int a = 0; a < static_cast<int>(orbit_count); ++a)
        for (int b = a; b < static_cast<int>(orbit_count); ++b)
        {
            double value = 0.0;
            read_j64_value(input, value, "OBME");
            result.one_body(a, b) = value;
            result.one_body(b, a) = value;
        }
    result.two_body.resize(jbasis.j2b_channel_number());
    for (int channel_index = 0;
         channel_index < static_cast<int>(jbasis.j2b_channel_number()); ++channel_index)
    {
        const auto &channel = jbasis.j2b_channel(channel_index);
        result.two_body[channel_index] = Eigen::MatrixXd::Zero(channel.size(), channel.size());
        for (int bra = 0; bra < static_cast<int>(channel.size()); ++bra)
        {
            const auto [a, b] = channel[bra];
            for (int ket = bra; ket < static_cast<int>(channel.size()); ++ket)
            {
                const auto [c, d] = channel[ket];
                const std::array<int, 5> expected_indices = {a, b, c, d, channel.J};
                for (int expected_value : expected_indices)
                {
                    std::int32_t value = 0;
                    read_j64_value(input, value, "TBME index");
                    if (value != expected_value)
                        throw std::runtime_error("jcoupled64 TBME ordering does not match the active basis");
                }
                double value = 0.0;
                read_j64_value(input, value, "TBME");
                result.two_body[channel_index](bra, ket) = value;
                result.two_body[channel_index](ket, bra) = value;
            }
        }
    }
    if (input.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("jcoupled64 payload contains trailing data");
    if (!std::isfinite(result.zero_body) || !result.one_body.allFinite())
        throw std::runtime_error("jcoupled64 payload contains non-finite zero- or one-body values");
    for (const Eigen::MatrixXd &channel : result.two_body)
        if (!channel.allFinite())
            throw std::runtime_error("jcoupled64 payload contains non-finite two-body values");
    return result;
}
} // namespace mrimsrg
