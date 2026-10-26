# BIGSTICK He4 Nrefmax=2 reference (test fixture)

He4 ground state in the emax=2 space with Nmax (W) = 2, N2LO_opt bare,
hw = 20 MeV, beta_cm = 0, computed with BIGSTICK 7.11.4 (OpenMP build,
`lanc_prec = 4`, i.e. single-precision Lanczos vectors).

- Interaction: `TwBME_N2LO_opt_bare_hw20_emax2_e2max4.MFDn` (NuHamil
  `file_convert` of the same me2j source as the minipack used by the
  simpleFCI fixture `prototype/mrimsrg/data/He4_Nrefmax2`), SHA-256
  `dca931389924e2f065c6d1162656bc735d2a13939501cbbebc838e1b9a1a36af`.
- `spectrum.in`: menu `n`, truncation `y` / `2`; writes `he4_nref2.wfn`.
  Ground state E = -17.01349 MeV, dimension 59.
- `basis.in`: menu `ba` with the same space and truncation; writes the ASCII
  `he4_nref2.bas` needed to read the wavefunction.

BIGSTICK's own two-body densities (`2d`) fail in this truncated space
("Problem cannot find sector jumps PP"), which is why
`prototype/mrimsrg/bigstick_jref.py` computes the densities from the
wavefunction.
