#ifndef MRReference_h
#define MRReference_h 1

#include "ModelSpace.hh"
#include "Operator.hh"
#include "TwoBodyME.hh"

#include <vector>

/// Explicit J_ref=0 multi-reference density context for MR-IMSRG(2).
///
/// The one-body density is diagonal in the natural-orbit basis and is stored as
/// one occupation per spherical orbit.  Lambda2 uses the normalized pair
/// convention of TwoBodyME.  This class does not silently alter ModelSpace:
/// production callers must construct the ModelSpace with matching occupations.
class MRReference
{
 public:
  ModelSpace *modelspace;
  int A;
  int Z;
  int Nrefmax;
  std::vector<double> occupations;
  TwoBodyME Lambda2;

  MRReference(ModelSpace &ms, int A_in, int Z_in, int Nrefmax_in,
              const std::vector<double> &occupations_in);

  bool OccupationsMatchModelSpace(double tolerance = 1e-12) const;
  double MaximumHermiticityViolation() const;
  double MaximumContractionViolation() const;
  void Validate(double tolerance = 1e-10) const;

  /// Contract a scalar two-body tensor with lambda2:
  /// 1/4 sum_abcd O_abcd lambda_abcd.
  double ContractLambda2(const TwoBodyME &two_body) const;

  /// Normal order an NN/NO2B vacuum operator with gamma1 and lambda2.
  Operator NormalOrder(const Operator &vacuum_operator) const;

  /// Convert a gamma1/lambda2-normal-ordered NN/NO2B operator to vacuum form.
  Operator UndoNormalOrder(const Operator &mr_operator) const;

 private:
  void CheckCompatibleOperator(const Operator &op) const;
};

#endif
