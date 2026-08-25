#include "MRReference.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

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
  if (!OccupationsMatchModelSpace(tolerance))
    throw std::runtime_error("MRReference occupations do not match ModelSpace occupations");

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
