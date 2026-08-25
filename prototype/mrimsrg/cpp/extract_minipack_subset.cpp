#include "HarmonicOscillator.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
using nucleus::HarmonicOscillator;
using nucleus::JBasisSet;
using nucleus::Orbit;
using nucleus::TwoBodyStates_J;

namespace
{
constexpr std::array<char, 8> minipack_magic{'m', 'i', 'n', 'i', 'p', 'a', 'c', 'k'};

struct Options
{
    fs::path parent;
    fs::path output;
    int child_emax = -1;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_extract_minipack_subset --parent FILE --output FILE --emax N");
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
        if (key == "--parent")
            options.parent = value;
        else if (key == "--output")
            options.output = value;
        else if (key == "--emax")
            options.child_emax = std::stoi(value);
        else
            usage_error("unknown option " + key);
    }
    if (options.parent.empty() || options.output.empty() || options.child_emax < 0)
        usage_error("--parent, --output, and nonnegative --emax are required");
    if (!fs::is_regular_file(options.parent))
        usage_error("parent must be a regular file");
    if (fs::exists(options.output))
        usage_error("refusing to overwrite output " + options.output.string());
    return options;
}

template <typename T>
T read_exact(std::ifstream &stream, const std::string &field)
{
    T value{};
    stream.read(reinterpret_cast<char *>(&value), sizeof(T));
    if (!stream)
        throw std::runtime_error("truncated parent while reading " + field);
    return value;
}

template <typename T>
void write_exact(std::ofstream &stream, const T &value, const std::string &field)
{
    stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
    if (!stream)
        throw std::runtime_error("failed while writing " + field);
}

std::uint64_t upper_triangle_count(const JBasisSet &basis)
{
    std::uint64_t count = 0;
    for (std::size_t channel = 0; channel < basis.j2b_channel_number(); ++channel)
    {
        const std::uint64_t dimension = basis.j2b_channel(channel).size();
        count += dimension * (dimension + 1) / 2;
    }
    return count;
}

bool retained(const Orbit &orbit, int child_emax)
{
    return 2 * orbit.n + orbit.l <= child_emax;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        std::ifstream input(options.parent, std::ios::binary);
        if (!input)
            throw std::runtime_error("cannot open parent " + options.parent.string());

        std::array<char, 8> magic{};
        input.read(magic.data(), magic.size());
        if (!input || magic != minipack_magic)
            throw std::runtime_error("parent has an invalid minipack header");
        const double hw = read_exact<double>(input, "oscillator frequency");
        const int parent_emax = read_exact<int>(input, "parent emax");
        const int declared_parent_records = read_exact<int>(input, "parent record count");
        if (options.child_emax >= parent_emax)
            throw std::runtime_error("child emax must be smaller than parent emax");

        HarmonicOscillator oscillator(hw);
        const JBasisSet parent_basis = oscillator.make_ho_basis(parent_emax, "Oslo");
        const JBasisSet child_basis = oscillator.make_ho_basis(options.child_emax, "Oslo");
        const std::uint64_t expected_parent_records = upper_triangle_count(parent_basis);
        const std::uint64_t expected_child_records = upper_triangle_count(child_basis);
        if (expected_parent_records != static_cast<std::uint64_t>(declared_parent_records))
            throw std::runtime_error("parent record count does not match its HO basis");
        if (expected_child_records > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
            throw std::overflow_error("child record count exceeds minipack header range");

        std::ofstream output(options.output, std::ios::binary | std::ios::out);
        if (!output)
            throw std::runtime_error("cannot create output " + options.output.string());
        output.write(minipack_magic.data(), minipack_magic.size());
        write_exact(output, hw, "oscillator frequency");
        write_exact(output, options.child_emax, "child emax");
        const int child_record_header = static_cast<int>(expected_child_records);
        write_exact(output, child_record_header, "child record count");

        std::uint64_t parent_records = 0;
        std::uint64_t child_records = 0;
        for (int tz : {-2, 0, 2})
            for (int parity : {0, 1})
                for (int angular_momentum = 0;
                     angular_momentum <= parent_basis.max_j(); ++angular_momentum)
                {
                    const int channel_index = TwoBodyStates_J::channel_index(
                        parity, tz, angular_momentum);
                    const auto &channel = parent_basis.j2b_channel(channel_index);
                    for (std::size_t bra_index = 0; bra_index < channel.size(); ++bra_index)
                    {
                        const auto [a, b] = channel[bra_index];
                        const Orbit &oa = parent_basis.orbit(a);
                        const Orbit &ob = parent_basis.orbit(b);
                        for (std::size_t ket_index = bra_index;
                             ket_index < channel.size(); ++ket_index)
                        {
                            const auto [c, d] = channel[ket_index];
                            const Orbit &oc = parent_basis.orbit(c);
                            const Orbit &od = parent_basis.orbit(d);
                            const float interaction = read_exact<float>(input, "interaction");
                            const bool has_kinetic_payload =
                                !oscillator.k1k2_is_zero(oa, ob, oc, od);
                            float hcom = 0.0F;
                            float pipj = 0.0F;
                            if (has_kinetic_payload)
                            {
                                hcom = read_exact<float>(input, "Hcm payload");
                                pipj = read_exact<float>(input, "p1.p2 payload");
                            }
                            ++parent_records;
                            if (!(retained(oa, options.child_emax) &&
                                  retained(ob, options.child_emax) &&
                                  retained(oc, options.child_emax) &&
                                  retained(od, options.child_emax)))
                                continue;
                            write_exact(output, interaction, "interaction");
                            if (has_kinetic_payload)
                            {
                                write_exact(output, hcom, "Hcm payload");
                                write_exact(output, pipj, "p1.p2 payload");
                            }
                            ++child_records;
                        }
                    }
                }

        if (parent_records != expected_parent_records)
            throw std::runtime_error("consumed parent record count is inconsistent");
        if (child_records != expected_child_records)
            throw std::runtime_error("selected child record count is inconsistent");
        if (input.peek() != std::ifstream::traits_type::eof())
            throw std::runtime_error("parent minipack contains trailing data");
        output.close();
        if (!output)
            throw std::runtime_error("failed while closing output");

        std::cout << "parent=" << fs::absolute(options.parent) << "\n"
                  << "output=" << fs::absolute(options.output) << "\n"
                  << "hw=" << hw << "\n"
                  << "parent_emax=" << parent_emax << "\n"
                  << "child_emax=" << options.child_emax << "\n"
                  << "parent_records=" << parent_records << "\n"
                  << "child_records=" << child_records << "\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_extract_minipack_subset: " << error.what() << "\n";
        return 1;
    }
}
