#ifndef MRThreeBody_h
#define MRThreeBody_h 1

#include "MRReference.hh"
#include "Operator.hh"

#include <string>

/// What the MR-NO2B normal ordering of the 3N interaction used.
struct MRThreeBodyReport
{
  std::string storage;              // "full" or "no2b"
  bool lambda3_in_reference = false;
  bool lambda3_term_included = false;
  std::string lambda3_reason;       // why the lambda3 zero-body term was left out
  double e3_gamma1_cubed = 0.0;     // (1/36) sum V3 A(g1 g1 g1)
  double e3_gamma1_lambda2 = 0.0;   // (1/36) sum V3 A(g1 lambda2)
  double e3_lambda3 = 0.0;          // (1/36) sum V3 lambda3
};

/// MR-NO2B contribution of the 3N interaction (Hergert 2017, Eqs. 41-43):
///   E3     = 1/36 sum V3 gamma3,  gamma3 = A(g1 g1 g1) + A(g1 lambda2) + lambda3
///   f3     = 1/4  sum V3 gamma2
///   Gamma3 =      sum V3 gamma1
/// with the reference densities; the residual normal-ordered three-body part
/// is discarded.  hbare carries the 3N in the original oscillator basis (any
/// 1B/2B parts are ignored); the result is in the natural-orbit basis, ready
/// to be added to MRReference::NormalOrder of the NN part.
///
/// The contractions are done in the oscillator basis (gamma1 = C n C^T and
/// lambda2, lambda3 rotated back), so the 3N matrix elements are never
/// transformed.  Everything except the lambda3 term uses only the
/// no2b-type matrix elements (ThreeBodyME::GetME_pn_no2b), so "no2b" and
/// "full" files give the same result; the lambda3 term needs general matrix
/// elements and is included only for full storage with a lambda3 reference.
Operator MRNormalOrderThreeBody(const Operator &hbare, const MRReference &reference,
                                MRThreeBodyReport &report, double tolerance = 1e-10);

#endif
