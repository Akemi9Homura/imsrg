#pragma once

#include <array>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace mrimsrg
{
inline constexpr const char *fixed_interaction_sha256 =
    "76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84";
inline constexpr const char *validated_emax4_interaction_sha256 =
    "d3dff5faa2a58d8c234914170caffa3649d0cd3805a1344818f4f4d3c37fd19e";

inline std::string shell_quote(const std::string &value)
{
    std::string quoted = "'";
    for (const char character : value)
    {
        if (character == '\'')
            quoted += "'\\''";
        else
            quoted += character;
    }
    quoted += "'";
    return quoted;
}

inline std::string sha256_file(const std::filesystem::path &path)
{
    const std::string command = "sha256sum -- " + shell_quote(path.string());
    std::FILE *pipe = popen(command.c_str(), "r");
    if (pipe == nullptr)
        throw std::runtime_error("cannot start sha256sum for " + path.string());

    std::array<char, 256> output{};
    const bool read_ok = std::fgets(output.data(), static_cast<int>(output.size()), pipe) != nullptr;
    const int status = pclose(pipe);
    if (!read_ok || status != 0)
        throw std::runtime_error("sha256sum failed for " + path.string());

    const std::string line(output.data());
    if (line.size() < 65 || !std::isspace(static_cast<unsigned char>(line[64])))
        throw std::runtime_error("unexpected sha256sum output for " + path.string());
    const std::string digest = line.substr(0, 64);
    for (const unsigned char character : digest)
        if (!std::isxdigit(character))
            throw std::runtime_error("invalid SHA-256 digest for " + path.string());
    return digest;
}

inline std::string require_fixed_interaction(const std::filesystem::path &path)
{
    const std::string digest = sha256_file(path);
    if (digest != fixed_interaction_sha256)
        throw std::runtime_error("unexpected interaction SHA-256: " + digest);
    return digest;
}

// The NCSM validator uses the interaction only to construct the verified HO
// orbit/channel tables before replacing every 0B/1B/2B matrix element with
// the lossless jcoupled64 payload.  Keep the production reference/export
// paths on require_fixed_interaction(), while permitting this independently
// frozen emax=4 basis source in the read-only spectral check.
inline std::string require_validated_interaction_basis(
    const std::filesystem::path &path)
{
    const std::string digest = sha256_file(path);
    if (digest != fixed_interaction_sha256 &&
        digest != validated_emax4_interaction_sha256)
        throw std::runtime_error("unexpected interaction SHA-256: " + digest);
    return digest;
}
} // namespace mrimsrg
