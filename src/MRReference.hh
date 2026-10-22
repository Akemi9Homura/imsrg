#ifndef MRReference_h
#define MRReference_h 1

#include "ModelSpace.hh"
#include "Operator.hh"
#include "TwoBodyME.hh"

#include <map>
#include <string>
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
  /// Columns are natural J-orbits expressed in the original spherical basis.
  arma::mat NaturalOrbitTransformation;
  int J2 = 0;
  int parity = 1;
  int emax = -1;
  int e2max = -1;
  double hw = 0.0;
  std::string interaction_sha256;
  std::string rdm_sha256;
  std::string wavefunction_sha256;

  MRReference(ModelSpace &ms, int A_in, int Z_in, int Nrefmax_in,
              const std::vector<double> &occupations_in);

  /// Read the self-describing little-endian mrimsrg_jref_v1 bridge.
  /// The ModelSpace must already carry the file occupations so its channel
  /// occupation weights have been constructed consistently.
  static MRReference ReadBinary(ModelSpace &ms, const std::string &filename,
                                double validation_tolerance = 1e-10);
  /// Write the complete self-describing little-endian mrimsrg_jref_v1 bridge.
  /// Existing files are never overwritten.
  void WriteBinary(const std::string &filename) const;
  /// Read only the orbit occupations so a caller can explicitly rebuild the
  /// ModelSpace occupation-dependent channel lists before ReadBinary().
  static std::map<index_t, double> ReadOccupationMap(
      ModelSpace &ms, const std::string &filename);

  /// Embed this fixed Nmax-truncated reference in a larger spherical model
  /// space. Existing orbit density blocks are preserved, while new orbits
  /// receive zero occupation/cumulant and an identity natural-orbit block.
  /// The target ModelSpace reference occupations are rebuilt explicitly.
  MRReference EmbedInModelSpace(
      ModelSpace &target_modelspace,
      const std::string &target_interaction_sha256,
      double validation_tolerance = 1e-10) const;

  bool OccupationsMatchModelSpace(double tolerance = 1e-12) const;
  double MaximumHermiticityViolation() const;
  double MaximumContractionViolation() const;
  void Validate(double tolerance = 1e-10) const;

  /// Contract a scalar two-body tensor with lambda2:
  /// 1/4 sum_abcd O_abcd lambda_abcd.
  double ContractLambda2(const TwoBodyME &two_body) const;

  /// Bytes in the dense numerical reference payload (occupations, NAT and lambda2).
  size_t DataSize() const;

  /// Normal order an NN/NO2B vacuum operator with gamma1 and lambda2.
  Operator NormalOrder(const Operator &vacuum_operator) const;

  /// Convert a gamma1/lambda2-normal-ordered NN/NO2B operator to vacuum form.
  Operator UndoNormalOrder(const Operator &mr_operator) const;

 private:
  void CheckCompatibleOperator(const Operator &op) const;
};

#endif
