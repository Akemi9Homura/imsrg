#include "MRReference.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace
{
constexpr std::array<char, 16> kReferenceMagic = {
    'm', 'r', 'i', 'm', 's', 'r', 'g', '_', 'j', 'r', 'e', 'f', '1', 0, 0, 0};
constexpr std::uint32_t kEndianMarker = 0x01020304u;

template <typename T>
T ReadScalar(std::istream &input, const char *description)
{
  static_assert(std::is_trivially_copyable<T>::value,
                "binary bridge scalars must be trivially copyable");
  T value{};
  input.read(reinterpret_cast<char *>(&value), sizeof(value));
  if (!input)
    throw std::runtime_error(std::string("truncated MR reference while reading ") + description);
  return value;
}

std::string ReadDigest(std::istream &input, const char *description)
{
  std::array<char, 64> bytes{};
  input.read(bytes.data(), bytes.size());
  if (!input)
    throw std::runtime_error(std::string("truncated MR reference while reading ") + description);
  const std::string value(bytes.begin(), bytes.end());
  if (!std::all_of(value.begin(), value.end(), [](char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
      }))
    throw std::runtime_error(std::string("invalid lowercase SHA-256 in ") + description);
  return value;
}
}

MRReference::MRReference(ModelSpace &ms, int A_in, int Z_in, int Nrefmax_in,
                         const std::vector<double> &occupations_in)
    : modelspace(&ms), A(A_in), Z(Z_in), Nrefmax(Nrefmax_in),
      occupations(occupations_in), Lambda2(&ms)
{
  if (A < 0 || Z < 0 || Z > A)
    throw std::invalid_argument("MRReference requires 0 <= Z <= A");
  if (Nrefmax < 0)
    throw std::invalid_argument("MRReference requires Nrefmax >= 0");
  if (occupations.size() != modelspace->GetNumberOrbits())
    throw std::invalid_argument("MRReference occupation count does not match ModelSpace");
  Lambda2.SetHermitian();
  NaturalOrbitTransformation.eye(modelspace->GetNumberOrbits(),
                                  modelspace->GetNumberOrbits());
}

MRReference MRReference::ReadBinary(ModelSpace &ms, const std::string &filename,
                                    double validation_tolerance)
{
  std::ifstream input(filename, std::ios::binary);
  if (!input)
    throw std::runtime_error("cannot open MR reference file: " + filename);
  std::array<char, 16> magic{};
  input.read(magic.data(), magic.size());
  if (!input || magic != kReferenceMagic)
    throw std::runtime_error("unsupported MR reference magic in " + filename);
  if (ReadScalar<std::uint32_t>(input, "endian marker") != kEndianMarker)
    throw std::runtime_error("MR reference has unsupported byte order");
  if (ReadScalar<std::uint32_t>(input, "schema version") != 1u)
    throw std::runtime_error("unsupported MR reference schema version");

  const int A = ReadScalar<std::int32_t>(input, "A");
  const int Z = ReadScalar<std::int32_t>(input, "Z");
  const int nrefmax = ReadScalar<std::int32_t>(input, "Nrefmax");
  const int J2 = ReadScalar<std::int32_t>(input, "J2");
  const int parity = ReadScalar<std::int32_t>(input, "parity");
  const int emax = ReadScalar<std::int32_t>(input, "emax");
  const int e2max = ReadScalar<std::int32_t>(input, "e2max");
  const double hw = ReadScalar<double>(input, "hw");
  const std::uint64_t norbits = ReadScalar<std::uint64_t>(input, "orbit count");
  const std::uint64_t nchannels = ReadScalar<std::uint64_t>(input, "channel count");
  const std::string interaction_sha256 = ReadDigest(input, "interaction digest");
  const std::string rdm_sha256 = ReadDigest(input, "RDM digest");
  const std::string wavefunction_sha256 = ReadDigest(input, "wavefunction digest");

  if (norbits != ms.GetNumberOrbits())
    throw std::runtime_error("MR reference orbit count does not match ModelSpace");
  if (nchannels > ms.GetNumberTwoBodyChannels())
    throw std::runtime_error("MR reference has more channels than ModelSpace");

  using OrbitLabel = std::array<int, 4>;
  std::map<OrbitLabel, index_t> model_orbits;
  for (index_t p : ms.all_orbits)
  {
    const Orbit &orbit = ms.GetOrbit(p);
    if (!model_orbits.emplace(OrbitLabel{orbit.n, orbit.l, orbit.j2, orbit.tz2}, p).second)
      throw std::runtime_error("ModelSpace contains duplicate spherical orbit labels");
  }
  std::vector<index_t> file_to_model(norbits);
  std::vector<double> occupations(ms.GetNumberOrbits(), 0.0);
  std::vector<bool> seen_model_orbit(ms.GetNumberOrbits(), false);
  for (std::uint64_t file_index = 0; file_index < norbits; ++file_index)
  {
    const int recorded_index = ReadScalar<std::int32_t>(input, "orbit index");
    const int n = ReadScalar<std::int32_t>(input, "orbit n");
    const int l = ReadScalar<std::int32_t>(input, "orbit l");
    const int j2 = ReadScalar<std::int32_t>(input, "orbit j2");
    const int tz2 = ReadScalar<std::int32_t>(input, "orbit tz2");
    const double occupation = ReadScalar<double>(input, "orbit occupation");
    if (recorded_index != static_cast<int>(file_index))
      throw std::runtime_error("MR reference orbit indices are not contiguous");
    const auto found = model_orbits.find({n, l, j2, tz2});
    if (found == model_orbits.end() || seen_model_orbit[found->second])
      throw std::runtime_error("MR reference orbit table does not map one-to-one to ModelSpace");
    file_to_model[file_index] = found->second;
    seen_model_orbit[found->second] = true;
    occupations[found->second] = occupation;
  }

  MRReference result(ms, A, Z, nrefmax, occupations);
  result.J2 = J2;
  result.parity = parity;
  result.emax = emax;
  result.e2max = e2max;
  result.hw = hw;
  result.interaction_sha256 = interaction_sha256;
  result.rdm_sha256 = rdm_sha256;
  result.wavefunction_sha256 = wavefunction_sha256;
  result.NaturalOrbitTransformation.zeros(norbits, norbits);
  for (std::uint64_t i = 0; i < norbits; ++i)
    for (std::uint64_t j = 0; j < norbits; ++j)
      result.NaturalOrbitTransformation(file_to_model[i], file_to_model[j]) =
          ReadScalar<double>(input, "natural-orbit transformation");

  std::vector<bool> seen_channel(ms.GetNumberTwoBodyChannels(), false);
  for (std::uint64_t file_channel = 0; file_channel < nchannels; ++file_channel)
  {
    const int J = ReadScalar<std::int32_t>(input, "channel J");
    const int channel_parity = ReadScalar<std::int32_t>(input, "channel parity");
    const int Tz = ReadScalar<std::int32_t>(input, "channel Tz");
    (void)ReadScalar<std::uint32_t>(input, "channel reserved field");
    const std::uint64_t npairs = ReadScalar<std::uint64_t>(input, "channel pair count");
    const size_t channel_index = ms.GetTwoBodyChannelIndex(J, channel_parity, Tz);
    if (channel_index >= ms.GetNumberTwoBodyChannels() || seen_channel[channel_index])
      throw std::runtime_error("MR reference has an invalid or duplicate two-body channel");
    seen_channel[channel_index] = true;
    const TwoBodyChannel &channel = ms.GetTwoBodyChannel(channel_index);
    if (npairs != channel.GetNumberKets())
      throw std::runtime_error("MR reference pair count does not match two-body channel");
    std::vector<size_t> file_to_local(npairs);
    std::vector<bool> seen_pair(npairs, false);
    for (std::uint64_t i = 0; i < npairs; ++i)
    {
      const std::uint32_t a_file = ReadScalar<std::uint32_t>(input, "pair first orbit");
      const std::uint32_t b_file = ReadScalar<std::uint32_t>(input, "pair second orbit");
      if (a_file >= norbits || b_file >= norbits)
        throw std::runtime_error("MR reference pair uses an invalid orbit index");
      const index_t a = file_to_model[a_file];
      const index_t b = file_to_model[b_file];
      const size_t local = channel.GetLocalIndex(std::min(a, b), std::max(a, b));
      if (local >= npairs || seen_pair[local])
        throw std::runtime_error("MR reference pair table does not map one-to-one to channel");
      file_to_local[i] = local;
      seen_pair[local] = true;
    }
    arma::mat &matrix = result.Lambda2.GetMatrix(channel_index);
    for (std::uint64_t i = 0; i < npairs; ++i)
      for (std::uint64_t j = 0; j < npairs; ++j)
        matrix(file_to_local[i], file_to_local[j]) =
            ReadScalar<double>(input, "lambda2 matrix element");
  }
  for (size_t ch = 0; ch < seen_channel.size(); ++ch)
    if (ms.GetTwoBodyChannel(ch).GetNumberKets() > 0 && !seen_channel[ch])
      throw std::runtime_error("MR reference omits a nonempty two-body channel");
  if (input.peek() != std::char_traits<char>::eof())
    throw std::runtime_error("MR reference has trailing bytes");
  result.Validate(validation_tolerance);
  return result;
}

std::map<index_t, double> MRReference::ReadOccupationMap(
    ModelSpace &ms, const std::string &filename)
{
  std::ifstream input(filename, std::ios::binary);
  if (!input)
    throw std::runtime_error("cannot open MR reference file: " + filename);
  std::array<char, 16> magic{};
  input.read(magic.data(), magic.size());
  if (!input || magic != kReferenceMagic)
    throw std::runtime_error("unsupported MR reference magic in " + filename);
  if (ReadScalar<std::uint32_t>(input, "endian marker") != kEndianMarker ||
      ReadScalar<std::uint32_t>(input, "schema version") != 1u)
    throw std::runtime_error("unsupported MR reference byte order or schema");
  for (int i = 0; i < 7; ++i)
    (void)ReadScalar<std::int32_t>(input, "reference metadata");
  (void)ReadScalar<double>(input, "hw");
  const std::uint64_t norbits = ReadScalar<std::uint64_t>(input, "orbit count");
  (void)ReadScalar<std::uint64_t>(input, "channel count");
  (void)ReadDigest(input, "interaction digest");
  (void)ReadDigest(input, "RDM digest");
  (void)ReadDigest(input, "wavefunction digest");
  if (norbits != ms.GetNumberOrbits())
    throw std::runtime_error("MR reference orbit count does not match ModelSpace");

  using OrbitLabel = std::array<int, 4>;
  std::map<OrbitLabel, index_t> model_orbits;
  for (index_t p : ms.all_orbits)
  {
    const Orbit &orbit = ms.GetOrbit(p);
    model_orbits.emplace(OrbitLabel{orbit.n, orbit.l, orbit.j2, orbit.tz2}, p);
  }
  std::map<index_t, double> occupations;
  for (std::uint64_t file_index = 0; file_index < norbits; ++file_index)
  {
    const int recorded_index = ReadScalar<std::int32_t>(input, "orbit index");
    const int n = ReadScalar<std::int32_t>(input, "orbit n");
    const int l = ReadScalar<std::int32_t>(input, "orbit l");
    const int j2 = ReadScalar<std::int32_t>(input, "orbit j2");
    const int tz2 = ReadScalar<std::int32_t>(input, "orbit tz2");
    const double occupation = ReadScalar<double>(input, "orbit occupation");
    if (recorded_index != static_cast<int>(file_index))
      throw std::runtime_error("MR reference orbit indices are not contiguous");
    const auto found = model_orbits.find({n, l, j2, tz2});
    if (found == model_orbits.end() || !occupations.emplace(found->second, occupation).second)
      throw std::runtime_error("MR reference orbit table does not map one-to-one to ModelSpace");
  }
  return occupations;
}

bool MRReference::OccupationsMatchModelSpace(double tolerance) const
{
  for (index_t p : modelspace->all_orbits)
  {
    if (std::abs(occupations[p] - modelspace->GetOrbit(p).occ) > tolerance)
      return false;
  }
  return true;
}

double MRReference::MaximumHermiticityViolation() const
{
  double maximum = 0.0;
  for (const auto &entry : Lambda2.MatEl)
  {
    const arma::mat &matrix = entry.second;
    if (matrix.n_rows != matrix.n_cols)
      return INFINITY;
    if (matrix.n_elem > 0)
      maximum = std::max(maximum, arma::abs(matrix - matrix.t()).max());
  }
  return maximum;
}

double MRReference::MaximumContractionViolation() const
{
  double maximum = 0.0;
  for (index_t p : modelspace->all_orbits)
  {
    const Orbit &op = modelspace->GetOrbit(p);
    for (index_t r : modelspace->all_orbits)
    {
      const Orbit &orr = modelspace->GetOrbit(r);
      if (op.l != orr.l || op.j2 != orr.j2 || op.tz2 != orr.tz2)
        continue;
      double contraction = 0.0;
      for (index_t q : modelspace->all_orbits)
      {
        const Orbit &oq = modelspace->GetOrbit(q);
        const int Jmin = std::max(std::abs(op.j2 - oq.j2),
                                  std::abs(orr.j2 - oq.j2)) /
                         2;
        const int Jmax = std::min(op.j2 + oq.j2, orr.j2 + oq.j2) / 2;
        for (int J = Jmin; J <= Jmax; ++J)
        {
          contraction += (2 * J + 1.0) / (op.j2 + 1.0) *
                         Lambda2.GetTBME_J(J, J, p, q, r, q);
        }
      }
      const double expected = p == r ? -occupations[p] * (1.0 - occupations[p]) : 0.0;
      maximum = std::max(maximum, std::abs(contraction - expected));
    }
  }
  return maximum;
}

void MRReference::Validate(double tolerance) const
{
  if (tolerance <= 0.0)
    throw std::invalid_argument("MRReference validation tolerance must be positive");
  if (J2 != 0 || parity != 1)
    throw std::runtime_error("MRReference requires a J=0 positive-parity reference state");
  if (emax >= 0 && emax != modelspace->GetEmax())
    throw std::runtime_error("MRReference emax does not match ModelSpace");
  if (e2max >= 0 && e2max != modelspace->GetE2max())
    throw std::runtime_error("MRReference e2max does not match ModelSpace");
  if (hw > 0.0 && std::abs(hw - modelspace->GetHbarOmega()) > tolerance)
    throw std::runtime_error("MRReference oscillator frequency does not match ModelSpace");
  if (!OccupationsMatchModelSpace(tolerance))
    throw std::runtime_error("MRReference occupations do not match ModelSpace occupations");

  const size_t norbits = modelspace->GetNumberOrbits();
  if (NaturalOrbitTransformation.n_rows != norbits ||
      NaturalOrbitTransformation.n_cols != norbits)
    throw std::runtime_error("MRReference natural-orbit transformation has the wrong dimension");
  if (arma::norm(NaturalOrbitTransformation.t() * NaturalOrbitTransformation -
                     arma::eye<arma::mat>(norbits, norbits),
                 "inf") > tolerance)
    throw std::runtime_error("MRReference natural-orbit transformation is not orthogonal");
  for (index_t p : modelspace->all_orbits)
    for (index_t q : modelspace->all_orbits)
    {
      const Orbit &op = modelspace->GetOrbit(p);
      const Orbit &oq = modelspace->GetOrbit(q);
      if ((op.l != oq.l || op.j2 != oq.j2 || op.tz2 != oq.tz2) &&
          std::abs(NaturalOrbitTransformation(p, q)) > tolerance)
        throw std::runtime_error("MRReference natural orbitals mix different spherical channels");
    }

  double particle_trace = 0.0;
  double proton_trace = 0.0;
  for (index_t p : modelspace->all_orbits)
  {
    const Orbit &op = modelspace->GetOrbit(p);
    const double occupation = occupations[p];
    if (!std::isfinite(occupation) || occupation < -tolerance || occupation > 1.0 + tolerance)
      throw std::runtime_error("MRReference occupation is outside [0,1]");
    particle_trace += (op.j2 + 1.0) * occupation;
    // imsrg++ uses tz2=-1 for proton orbits and +1 for neutron orbits.
    if (op.tz2 < 0)
      proton_trace += (op.j2 + 1.0) * occupation;
  }
  if (std::abs(particle_trace - A) > tolerance)
    throw std::runtime_error("MRReference Tr(gamma1) does not equal A");
  if (std::abs(proton_trace - Z) > tolerance)
    throw std::runtime_error("MRReference proton trace does not equal Z");

  const double hermiticity_violation = MaximumHermiticityViolation();
  if (hermiticity_violation > tolerance)
  {
    std::ostringstream message;
    message << "MRReference lambda2 is not Hermitian; max violation="
            << hermiticity_violation;
    throw std::runtime_error(message.str());
  }
  const double contraction_violation = MaximumContractionViolation();
  if (contraction_violation > tolerance)
  {
    std::ostringstream message;
    message << "MRReference lambda2 contraction is inconsistent with gamma1; max violation="
            << contraction_violation;
    throw std::runtime_error(message.str());
  }
}

double MRReference::ContractLambda2(const TwoBodyME &two_body) const
{
  if (two_body.modelspace != modelspace)
    throw std::invalid_argument("lambda2 contraction uses different ModelSpace objects");
  if (two_body.rank_J != 0 || two_body.rank_T != 0 || two_body.parity != 0)
    throw std::invalid_argument("lambda2 contraction requires a scalar parity-conserving two-body tensor");

  double contraction = 0.0;
  for (const auto &entry : two_body.MatEl)
  {
    const auto lambda_entry = Lambda2.MatEl.find(entry.first);
    if (lambda_entry == Lambda2.MatEl.end())
      throw std::runtime_error("lambda2 and operator two-body channel layouts differ");
    const arma::mat &matrix = entry.second;
    const arma::mat &lambda = lambda_entry->second;
    if (matrix.n_rows != lambda.n_rows || matrix.n_cols != lambda.n_cols)
      throw std::runtime_error("lambda2 and operator two-body block dimensions differ");
    const int J = modelspace->GetTwoBodyChannel(entry.first[0]).J;
    contraction += (2 * J + 1.0) * arma::accu(matrix % lambda.t());
  }
  return contraction;
}

size_t MRReference::DataSize() const
{
  return occupations.size() * sizeof(double) +
         NaturalOrbitTransformation.n_elem * sizeof(double) + Lambda2.size();
}

void MRReference::CheckCompatibleOperator(const Operator &op) const
{
  if (op.modelspace != modelspace)
    throw std::invalid_argument("MRReference and Operator use different ModelSpace objects");
  if (op.GetJRank() != 0 || op.GetTRank() != 0 || op.GetParity() != 0)
    throw std::invalid_argument("MR normal ordering currently supports scalar operators only");
  if (op.GetParticleRank() > 2)
    throw std::invalid_argument("MR normal ordering currently supports NN/NO2B operators only");
  if (!OccupationsMatchModelSpace())
    throw std::runtime_error("MRReference occupations do not match ModelSpace occupations");
}

Operator MRReference::NormalOrder(const Operator &vacuum_operator) const
{
  CheckCompatibleOperator(vacuum_operator);
  Operator result = vacuum_operator.DoNormalOrdering2(+1, modelspace->all_orbits);
  result.ZeroBody += ContractLambda2(vacuum_operator.TwoBody);
  return result;
}

Operator MRReference::UndoNormalOrder(const Operator &mr_operator) const
{
  CheckCompatibleOperator(mr_operator);
  Operator result = mr_operator.DoNormalOrdering2(-1, modelspace->all_orbits);
  result.ZeroBody -= ContractLambda2(mr_operator.TwoBody);
  return result;
}
