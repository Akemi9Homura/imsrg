#include "Hamiltonian.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
using nucleus::Hamiltonian;

namespace
{
struct Options
{
    fs::path parent;
    fs::path child;
    int mass_number = 0;
};

[[noreturn]] void usage_error(const std::string &message)
{
    throw std::invalid_argument(
        message +
        "\nusage: mrimsrg_check_minipack_subset --parent FILE --child FILE --A MASS");
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
        else if (key == "--child")
            options.child = value;
        else if (key == "--A")
            options.mass_number = std::stoi(value);
        else
            usage_error("unknown option " + key);
    }
    if (options.parent.empty() || options.child.empty() || options.mass_number <= 0)
        usage_error("--parent, --child, and positive --A are required");
    if (!fs::is_regular_file(options.parent) || !fs::is_regular_file(options.child))
        usage_error("parent and child must be regular files");
    if (fs::equivalent(options.parent, options.child))
        usage_error("parent and child must be different files");
    return options;
}

struct Difference
{
    double maximum = 0.0;
    std::string rank = "zero_body";
    std::size_t channel = 0;
    Eigen::Index row = 0;
    Eigen::Index column = 0;
    double parent_value = 0.0;
    double child_value = 0.0;
};

Difference maximum_difference(const Hamiltonian &parent, const Hamiltonian &child)
{
    const auto &parent_basis = parent.get_jbasis();
    const auto &child_basis = child.get_jbasis();
    if (parent_basis.orbit_number() != child_basis.orbit_number() ||
        parent_basis.j1b_channel_number() != child_basis.j1b_channel_number() ||
        parent_basis.j2b_channel_number() != child_basis.j2b_channel_number())
        throw std::runtime_error("truncated parent and child basis dimensions differ");

    for (std::size_t orbit = 0; orbit < parent_basis.orbit_number(); ++orbit)
    {
        const auto &left = parent_basis.orbit(static_cast<int>(orbit));
        const auto &right = child_basis.orbit(static_cast<int>(orbit));
        if (left.n != right.n || left.l != right.l || left.j != right.j ||
            left.tz != right.tz)
            throw std::runtime_error("truncated parent and child orbit tables differ");
    }

    Difference result;
    result.maximum = std::abs(parent.get_zero_body() - child.get_zero_body());
    result.parent_value = parent.get_zero_body();
    result.child_value = child.get_zero_body();
    for (std::size_t a = 0; a < parent_basis.orbit_number(); ++a)
        for (std::size_t b = 0; b < parent_basis.orbit_number(); ++b)
        {
            const double left = parent.get_j1b(static_cast<int>(a), static_cast<int>(b));
            const double right = child.get_j1b(static_cast<int>(a), static_cast<int>(b));
            const double difference = std::abs(left - right);
            if (difference > result.maximum)
                result = {difference, "one_body", 0,
                          static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b),
                          left, right};
        }

    for (std::size_t channel = 0; channel < parent_basis.j2b_channel_number(); ++channel)
    {
        const auto &left = parent.get_j2b_vmat_of_channel(channel);
        const auto &right = child.get_j2b_vmat_of_channel(channel);
        if (left.rows() != right.rows() || left.cols() != right.cols())
            throw std::runtime_error("truncated parent and child J-channel dimensions differ");
        if (left.size() > 0)
        {
            Eigen::Index row = 0;
            Eigen::Index column = 0;
            const double difference =
                (left - right).cwiseAbs().maxCoeff(&row, &column);
            if (difference > result.maximum)
                result = {difference, "two_body", channel, row, column,
                          left(row, column), right(row, column)};
        }
    }
    return result;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        nucleus::log_init(nucleus::LogLevel::Info);
        const Options options = parse_options(argc, argv);
        Hamiltonian parent;
        parent.read_minipack(options.parent.string(), options.mass_number, 0.0);
        const int parent_emax = parent.get_jbasis().max_2n_l();

        Hamiltonian child;
        child.read_minipack(options.child.string(), options.mass_number, 0.0);
        const int child_emax = child.get_jbasis().max_2n_l();
        if (child_emax >= parent_emax)
            throw std::runtime_error("child emax must be smaller than parent emax");
        if (parent.get_jbasis().get_hbar_omega() !=
            child.get_jbasis().get_hbar_omega())
            throw std::runtime_error("parent and child oscillator frequencies differ");

        parent.truncate(child_emax);
        const Difference difference = maximum_difference(parent, child);
        std::cout << std::setprecision(17)
                  << "parent=" << fs::absolute(options.parent) << "\n"
                  << "child=" << fs::absolute(options.child) << "\n"
                  << "A=" << options.mass_number << "\n"
                  << "hw=" << child.get_jbasis().get_hbar_omega() << "\n"
                  << "parent_emax=" << parent_emax << "\n"
                  << "child_emax=" << child_emax << "\n"
                  << "j_orbits=" << child.get_jbasis().orbit_number() << "\n"
                  << "j_channels=" << child.get_jbasis().j2b_channel_number() << "\n"
                  << "maximum_abs_mev=" << difference.maximum << "\n"
                  << "worst_rank=" << difference.rank << "\n"
                  << "worst_channel=" << difference.channel << "\n"
                  << "worst_row=" << difference.row << "\n"
                  << "worst_column=" << difference.column << "\n"
                  << "parent_value_mev=" << difference.parent_value << "\n"
                  << "child_value_mev=" << difference.child_value << "\n";
        if (difference.maximum != 0.0)
            throw std::runtime_error("child is not an exact minipack restriction of parent");
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "mrimsrg_check_minipack_subset: " << error.what() << "\n";
        return 1;
    }
}
