#pragma once

#include <eigen3/Eigen/Core>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mrimsrg
{
namespace fs = std::filesystem;

struct VacuumMScheme
{
    std::uint64_t norb = 0;
    double zero_body = 0.0;
    Eigen::MatrixXd one_body;
    std::vector<double> two_body;

    double two(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const
    {
        return two_body[((p * norb + q) * norb + r) * norb + s];
    }
};

template <typename T>
void read_exact(std::ifstream &stream, T *destination, std::size_t count, const char *field)
{
    stream.read(reinterpret_cast<char *>(destination),
                static_cast<std::streamsize>(count * sizeof(T)));
    if (!stream)
        throw std::runtime_error(std::string("truncated bridge payload while reading ") + field);
}

inline VacuumMScheme read_vacuum_mscheme(const fs::path &path)
{
    std::ifstream payload(path, std::ios::binary);
    if (!payload)
        throw std::runtime_error("cannot open bridge payload: " + path.string());

    std::array<char, 16> magic{};
    read_exact(payload, magic.data(), magic.size(), "magic");
    const std::array<char, 16> expected = {'m', 'r', 'i', 'm', 's', 'r', 'g', '_',
                                           'm', '_', 'v', '1', 0,   0,   0,   0};
    if (magic != expected)
        throw std::runtime_error("unsupported MR-IMSRG bridge payload");

    VacuumMScheme result;
    read_exact(payload, &result.norb, 1, "norb");
    read_exact(payload, &result.zero_body, 1, "zero body");
    result.one_body.resize(static_cast<Eigen::Index>(result.norb),
                           static_cast<Eigen::Index>(result.norb));
    for (std::size_t p = 0; p < result.norb; ++p)
        for (std::size_t q = 0; q < result.norb; ++q)
            read_exact(payload, &result.one_body(static_cast<Eigen::Index>(p),
                                                 static_cast<Eigen::Index>(q)),
                       1, "one body");

    const std::size_t n2 = static_cast<std::size_t>(result.norb) * result.norb;
    if (result.norb != 0 && n2 / result.norb != result.norb)
        throw std::overflow_error("bridge payload orbit count overflows tensor dimensions");
    const std::size_t n4 = n2 * n2;
    if (n2 != 0 && n4 / n2 != n2)
        throw std::overflow_error("bridge payload tensor size overflows address space");
    result.two_body.resize(n4);
    read_exact(payload, result.two_body.data(), result.two_body.size(), "two body");
    if (payload.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("bridge payload contains trailing data");

    if (!std::isfinite(result.zero_body) || !result.one_body.allFinite())
        throw std::runtime_error("bridge payload contains non-finite zero- or one-body values");
    for (double value : result.two_body)
        if (!std::isfinite(value))
            throw std::runtime_error("bridge payload contains non-finite two-body values");
    return result;
}

inline fs::path resolve_vacuum_payload(const fs::path &flow_output)
{
    if (fs::is_directory(flow_output))
        return flow_output / "vacuum_mscheme.bin";
    return flow_output;
}
} // namespace mrimsrg
