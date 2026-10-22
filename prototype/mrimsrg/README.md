# Rapid MR-IMSRG prototype

This directory contains the deliberately small m-scheme prototype described in
`docs/MR-IMSRG-快速结果计划.md`.  The C++ bridge does not implement another
interaction reader or NCSM solver: it builds against the validated
`/home/mengziyan/simple-ncsm` `Hamiltonian` and `simpleFCI` libraries.  This
lightweight path is only for deterministic bridge smoke tests; formal NCSM
spectrum acceptance uses BIGSTICK.

Build the input bridge:

```bash
cmake -S prototype/mrimsrg -B prototype/mrimsrg/build \
  -DSIMPLE_NCSM_ROOT=/home/mengziyan/simple-ncsm \
  -DCMAKE_BUILD_TYPE=Release
cmake --build prototype/mrimsrg/build \
  --target mrimsrg_prepare mrimsrg_validate mrimsrg_export_no2bpack -j2
```

Prepare and check the first He4 reference:

```bash
prototype/mrimsrg/build/mrimsrg_prepare \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --output prototype/mrimsrg/data/He4_Nrefmax0 --Z 2 --N 2 --nrefmax 0
python3 prototype/mrimsrg/densities.py prototype/mrimsrg/data/He4_Nrefmax0
```

Generated Hamiltonians, wavefunctions and densities belong under `data/` and
are ignored by Git.  The bridge refuses to overwrite an existing bundle.
Both C++ entry points compute the fixed interaction's SHA-256 with the standard
`sha256sum` utility before reading it; the Python runner independently repeats
the content check, including when the file has been relocated.

Run the dependency-free unit tests with:

```bash
python3 -m unittest discover -s prototype/mrimsrg/tests -v
```

For the strict single-reference degeneration gate, build the current
`pyIMSRG` binding (the extra denominator methods are read-only diagnostics),
generate a float64 `jcoupled64` file from the `s=0` checkpoint, and run:

```bash
PYTHONPATH=prototype/mrimsrg python3 prototype/mrimsrg/sr_imsrgpp_check.py \
  --nucleus He4 \
  --reference prototype/mrimsrg/data/He4_Nrefmax0_final \
  --jcoupled64 /tmp/He4_s0.jcoupled64
```

The checker verifies the compiled oracle's core source hashes, reconstructs
the complete m-scheme tensors, and compares the production MR path against
actual calls to `Operator::DoNormalOrdering`, `Generator::Update`, and
`Commutator::Commutator`.  It checks the relaxed mask, every selected EN
denominator, `H/eta/RHS` at `s=0`, and the same quantities after a shared
`ds=1e-4` Euler step.  It then compares three shared fixed-step RK4
checkpoints at `s=0.001,0.002,0.003`.  It exits nonzero if any rank exceeds
`1e-10 MeV`.  The RHS check also compares each of the eight named SR
contractions to the corresponding current C++ routine, so cancellation
between independently wrong terms cannot pass the gate.

For a completed fixed-``s`` production flow, add ``--production-flow`` and
set ``--full-flow-ode-tolerance`` to the common production ``rtol=atol``.
The checker then runs the current C++ direct-flow solver to the same ``s`` and
compares final ``H``, ``eta``, RHS, every selected denominator, and the two
explicitly different norm conventions.  Point7 scripts for this mode must be
created by ``gen_flow_job.py --sr-check-flow ...`` just like production flows.

For a correlated-reference fixed-``s`` flow, use the production J-scheme
validator instead. It exports the same NCSM reference to ``jref``, performs
HO→NAT, MR normal ordering, the C++ adaptive direct flow, MR de-normal
ordering, and NAT→HO, then compares final ``H/eta/RHS`` and the ordinary
vacuum Hamiltonian against the Python m-scheme flow:

```bash
PYTHONPATH=prototype/mrimsrg python3 \
  prototype/mrimsrg/mr_imsrgpp_flow_check.py \
  --nucleus He4 --reference prototype/mrimsrg/data/He4_Nrefmax2 \
  --interaction /path/to/fixed.minipack --production-flow /path/to/flow \
  --pyimsrg-dir build/src --ode-tolerance 1e-9 \
  --output-jcoupled64 /path/to/cpp_vacuum.jcoupled64 --json /path/to/report.json
```

The production and C++ tolerances must match and the Python flow must reach
its fixed target ``s`` without an early residual stop. The checker evaluates
the correlated generator and RHS in the temporary NAT basis and transforms
them back to HO for comparison; attempting to evaluate them directly with a
non-diagonal HO-basis ``gamma1`` remains an error.
On point7 this validator must be submitted through the same single-point job
generator, using ``gen_flow_job.py --mr-check-flow FLOW --pyimsrg-dir
build/src``; the generated job deliberately uses system Python matching the
current ``pyIMSRG`` build rather than loading miniconda. If compute nodes do
not share the system NumPy installation, create a separate ABI-matched build
and pass both explicitly, for example ``--pyimsrg-dir build-py312/src
--mr-check-python /opt/library/miniconda-3.12.9/bin/python3``; never load
miniconda around a system-Python extension.

`commutator.py` implements the natural-orbital MR-IMSRG(2) commutator with
`lambda3=0`.  Its equations were regenerated independently with QCombo 0.2.0:
for a random four-orbit correlated reference, direct summation of the QCombo
0B/1B/2B expressions agreed with the numerical implementation to
`1.8e-15`.  The committed tests additionally compare the correlated 1B--2B
zero-body term and the single-reference single/double-excitation blocks with
explicit Fock-space matrix commutators. Direct calls reject a non-diagonal
`gamma1` instead of silently treating its diagonal as natural occupations.
For the required `He4/O16, Nrefmax=2` flows, `flow.py` diagonalizes only the
connected off-diagonal blocks of `gamma1` and evaluates the published equations
and `Delta e` mask in that temporary spherical natural-orbital basis, as stated
in Vobig Sec. 6.5.3.  The natural orbitals inherit the `e` labels of their
output slots. Final and checkpoint Hamiltonians are transformed back to the
unchanged input HO orbit ordering; the exact
orthogonal matrix is saved as `formula_basis_vectors.npy`. This is an internal
basis-covariant evaluation, not a natural-orbital NCSM optimization or a
replacement of the HO quantum-number labels.

`generator.py` implements one generator only: the modified White generator
used in the IM-NCSM literature, restricted by the relaxed `Delta e != 0`
one- and two-body masks. It follows the production `White-NCSM` truncation
defined in Vobig Sec. 6.5.4: all irreducible-density terms are omitted from
the generator numerator and the displayed `O(lambda2)` terms are omitted
from the Epstein--Nesbet denominators. The MR-IMSRG(2) commutator itself still
retains `lambda2`. Denominators with magnitude below `1e-6 MeV` are replaced
by the positive cutoff, matching the current `src/Generator.cc` behavior;
the diagonal one-body values, occupations, and two-body diagonal terms are
first reduced to spherical-orbit scalars.  In particular, the ordered
magnetic-substate average of `Gamma[p,q,p,q]` is exactly the unnormalized
`src/TwoBodyME.cc::GetTBMEmonopole()` convention.  This is essential: using
individual m-scheme diagonals in the denominator makes `eta` non-scalar even
when `H` and the `J=0` reference are scalar.  In the Slater limit the resulting
denominators reduce term by term to `src/Generator.cc`. Independently, the
strict directional `D1/D2=<Psi|H:A:|Psi>` diagnostic retains the published terms linear in
`lambda2` with `lambda3=0`. Its anti-Hermitian combination is monitored as
the strict diagnostic and agrees with `<Psi|[H,:A:]|Psi>` from the QCombo and
explicit-commutator checks.

`flow.py` directly integrates `dH/ds=[eta,H]` with an adaptive DOP853
stepper. Like the existing C++ `IMSRGSolver`, it updates `eta` from every
trial Hamiltonian and monitors the norm of the actual masked White-NCSM
generator in Vobig Eqs. (6.5.28)--(6.5.29). The C++ solver uses an absolute
`Eta.Norm()` threshold; this prototype deliberately applies the plan's
stricter normalized gate and stops when that norm falls to the requested
fraction of its initial value. The ODE state stores only independent
Hermitian one-body and antisymmetric/Hermitian pair-space two-body elements;
this mirrors the established operator storage and prevents redundant tensor
components from developing symmetry-violating numerical modes. The default
`max_step=10` is only an upper bound: DOP853 continues to reduce trial steps
to satisfy `rtol/atol`.

Every trajectory point stores three separately labeled quantities:
`generator_residual_ratio` is the norm of the denominator-weighted,
anti-Hermitian generator and is the formal acceptance condition;
`generator_numerator_residual_ratio` is the unweighted lambda-free numerator;
and `residual_ratio` is the strict lambda2-dependent `D-D^dagger` diagnostic.
The latter two are diagnostics and need not vanish at the approximate
generator's fixed point. This matches the `Eta.Norm()` convergence criterion
of the existing C++ `IMSRGSolver` and avoids claiming that the stricter
quantity vanished when it did not.

Run and materialize one calculation with:

```bash
PYTHONPATH=prototype/mrimsrg python3 prototype/mrimsrg/run_mrimsrg.py \
  prototype/mrimsrg/data/He4_Nrefmax0_final \
  prototype/mrimsrg/data/He4_flow_rtol1e-6 \
  --checkpoint-s 1.0
```

The output directory contains `gamma1/gamma2/lambda2`, initial and final
MR-normal-ordered tensors, and the final ordinary vacuum-normal-ordered
`E0+t+V` tensors.  `metadata.json` records the fixed interaction identity,
reference metadata, `lambda3=0`, generator/mask, ODE settings and every
accepted-step residual.  Existing output directories are never overwritten.
The output also materializes the ordinary vacuum Hamiltonian at `s=0` and at
the requested intermediate `--checkpoint-s` under `checkpoints/`; the root
payload is the final point.  Thus a production run retains the three points
required by the acceptance plan without keeping every large ODE state.
The command returns status 2 after saving if `smax` is reached before the
White-NCSM generator target, allowing a Slurm job to distinguish a diagnostic from an
accepted result.

If a long correlated-reference flow reaches `smax`, continue from its saved
MR-normal-ordered Hamiltonian without restarting or renormalizing the gate:

```bash
PYTHONPATH=prototype/mrimsrg python3 prototype/mrimsrg/run_mrimsrg.py \
  prototype/mrimsrg/data/He4_Nrefmax2 next-flow \
  --resume-from prior-flow --smax 50000 --checkpoint-s 40000
```

The continuation verifies the generator implementation identifier, reference
metadata, and exact `gamma1/gamma2/lambda2` arrays match, rejects legacy
outputs that predate the current spherical-monopole denominator or separate
`Rgen/Rnum` records, and retains the original three residual norms as the
ratio denominators. `smax` and `checkpoint-s` are absolute flow parameters.
The dedicated point7 generator accepts the same `--resume-from` option and
still emits exactly one inspected Slurm job.

Export a saved NCSM reference to the compact production C++ J-scheme input:

```bash
python3 -m prototype.mrimsrg.export_jref \
  --reference prototype/mrimsrg/data/Be8_Nrefmax0_final \
  --output /tmp/Be8_Nrefmax0.jref
```

The `mrimsrg_jref_v1` payload records provenance hashes, physical metadata,
the spherical orbit table, natural-orbit transformation and occupations, and
float64 normalized-pair `lambda2` blocks.  The C++ reader maps by
`(n,l,j2,tz2)` and validates the scalar density and cumulant contraction; it
does not load the full m-scheme density during a production flow. Existing
output files are never overwritten.

Embed that same fixed `Nrefmax` reference in a larger Hamiltonian space with
zero density/cumulant on the added orbits and an identity added NAT block:

```bash
python3 -m prototype.mrimsrg.embed_jref \
  --source /tmp/He4_Nrefmax2_emax2.jref \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax4_e2max8.minipack \
  --output /tmp/He4_Nrefmax2_emax4.jref \
  --pyimsrg-dir build/src \
  --json /tmp/He4_Nrefmax2_emax4.json
```

The command verifies the minipack header and SHA-256, maps all old spherical
orbits by `(n,l,j2,tz2)`, validates the embedded density in C++, writes the
complete larger J-scheme channel layout, and reads it back independently. It
does not claim that the NCSM reference was rediagonalized in a larger
`Nrefmax` space.

Convert the ordinary A-dependent bare minipack to the lossless J64 input read
directly by the production `imsrg++` driver, without constructing m-scheme
matrix elements:

```bash
prototype/mrimsrg/build/mrimsrg_minipack_to_j64 \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax4_e2max8.minipack \
  --output /tmp/He4_NNLOopt_emax4.jcoupled64 --A 4
```

This bridge reuses `shell-model-obs::Hamiltonian::read_minipack()` for the
same A-dependent intrinsic kinetic energy as the NCSM path, copies only its
J-channel matrices, and requires an exact J64 write/read round trip. The J64
file already contains the complete vacuum Hamiltonian, so the IMSRG driver
must use `fmt2=jcoupled64` and must not add `Trel_Op` again.

For a larger HO child space, preserve the bare minipack payload directly
instead of reading an A-dependent Hamiltonian and writing it back:

```bash
prototype/mrimsrg/build/mrimsrg_extract_minipack_subset \
  --parent /path/to/TwBME_N2LO_opt_bare_hw20_emax14_e2max28.minipack \
  --output /path/to/TwBME_N2LO_opt_hw20_emax6_e2max12.minipack \
  --emax 6

prototype/mrimsrg/build/mrimsrg_check_minipack_subset \
  --parent /path/to/TwBME_N2LO_opt_bare_hw20_emax14_e2max28.minipack \
  --child /path/to/TwBME_N2LO_opt_hw20_emax6_e2max12.minipack --A 4
```

The extractor copies the original float32 interaction and optional Hcm/
`p1.p2` records selected by the existing Oslo HO channel layout, verifies
both record counts and rejects trailing data. The checker independently reads
both files with `shell-model-obs`, calls its production
`Hamiltonian::truncate()`, and demands exact 0B/1B/2B agreement. Run the
checker at two different masses to constrain both the bare interaction and
the A-dependent intrinsic-kinetic payload. Do not implement this operation as
`read_minipack(A) -> write_minipack()`: the read Hamiltonian already contains
the A-dependent kinetic term and is not the original bare payload.

Read a materialized ordinary Hamiltonian back into the existing NCSM solver:

```bash
prototype/mrimsrg/build/mrimsrg_validate \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --flow-output prototype/mrimsrg/data/He4_flow_rtol1e-6 \
  --Z 2 --N 2 --nmax 8 --states 3
```

The validator uses the original interaction only to initialize the identical
orbit/channel ordering, replaces its m-scheme matrix elements with the saved
vacuum `E0+t+V`, and then calls the existing `simpleFCI`.  At `s=0` this path
reproduces the complete 3060-dimensional He4 result as
`-20.3388325043 MeV, J=0`. By default the validator requests the three
lowest states (or the complete space when its dimension is smaller) and
prints each energy, excitation energy, and `2J`; `--states` changes that
small readback list without introducing another solver.

Export the same vacuum Hamiltonian in the packed J-coupled format understood
by both NCSM and FCIQMC:

```bash
prototype/mrimsrg/build/mrimsrg_export_no2bpack \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --jcoupled64 /tmp/He4_mrimsrg.jcoupled64 \
  --output He4_mrimsrg.no2bpack --Z 2 --N 2 \
  --diagnostic-jcoupled64 /tmp/He4_mrimsrg_roundtrip.jcoupled64
```

The exporter uses the fixed interaction only for the verified HO orbit table.
For production C++ output, `--jcoupled64` reads the lossless J-scheme file
directly and maps its explicit `(n,l,2j,2tz)` orbit table and
`(a,b,c,d,J)` TBME records to the independent `shell-model-obs` channel
ordering. It does not assume that the two programs enumerate channels or
pairs identically. The older `--flow-output` alternative projects a final
vacuum m-scheme tensor onto scalar J-coupled matrix elements with the same
`shell-model-obs` Clebsch--Gordan, pair normalization and phase conventions;
that path reconstructs the complete m-scheme tensors and refuses to write
unless the maximum discrepancy is below `--scalar-tolerance` (default
`1e-9 MeV`). A checkpoint directory containing `vacuum_mscheme.bin` can be
passed in place of the final flow directory. Select exactly one of
`--jcoupled64` and `--flow-output`.

The optional `--diagnostic-jcoupled64` output stores the same J-coupled OBMEs
and TBMEs as float64, with an identifying magic header. It is an acceptance
format, not a downstream interaction format and must not be passed to a
`no2bpack` reader. Read it back, reconstruct m-scheme matrix elements from the
serialized J-coupled values, and diagonalize with the same NCSM solver using:

```bash
prototype/mrimsrg/build/mrimsrg_validate \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --jcoupled64 /tmp/He4_mrimsrg.jcoupled64 \
  --Z 2 --N 2 --nmax 8 --states 3
```

As an end-to-end rotational-symmetry regression, a fresh He4 flow to
`s=0.001` reconstructs both the one- and two-body m-scheme tensors from the
projected J-coupled data with maximum errors `3.55e-15 MeV`.  At `Nmax=8`,
the direct double-precision dense path gives `-20.3396323958 MeV` and the
packed reader gives `-20.3396333376 MeV`; their `9.42e-7 MeV` difference is
`0.000942 keV` (`0.942 eV`), the measured float32 OBME/TBME packing effect.
The float64 J-coupled readback differs from the direct dense energies by at
most `1.8e-14 MeV` across the three printed states. Thus the
Clebsch--Gordan phases, identical-pair normalization, channel ordering and
inverse coupling pass the spectral check at double-precision numerical noise;
production output remains the existing float32 `no2bpack`.
The independent `/home/mengziyan/fciqmc/fciqmc-mpi` production reader also
completed a one-step He4 FCIQMC initialization from this exact file with
`int_format=no2bpack`.  Its native NCSM executable reproduced the same three
levels as `-20.3396333376`, `-16.6810264996`, and `-16.6581354522 MeV`
(`2J=0,0,4`).

Read the exported file through the independent native `no2bpack` reader and
the same NCSM solver with:

```bash
prototype/mrimsrg/build/mrimsrg_validate \
  --no2bpack He4_mrimsrg.no2bpack \
  --Z 2 --N 2 --nmax 8 --states 3
```

Use `beta_cm=0` in downstream NCSM/FCIQMC calculations. The format stores its
zero-body term as `double`, but OBMEs and TBMEs as `float`; retain the original
NPY/dense payload as the internal double-precision result and use spectral
readback to quantify the expected packed-format rounding. Every nucleus has
an A-dependent intrinsic Hamiltonian, so files exported for different nuclei
must not be interchanged.

The bridge requires NumPy and SciPy on the Python side.  The first checked
integration run used `shell-model-obs` revision `1687f16` and reproduced the
full He4 `emax=2` ground-state benchmark as `-20.3388325043 MeV`; contraction
of the independently constructed `gamma1/gamma2` gave the same value.
The dense m-scheme readback interface is revision `4b10afa` of that dependency.

On point7, generate and inspect exactly one Slurm job at a time with:

```bash
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 --submit
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 --nrefmax 2
```

The dedicated generator is for this Python prototype; it deliberately does
not alter the existing production `gen_job.py` used by `imsrg++`.  It fixes the
documented point7 paths, one node/task, 64 cores, QOS/log conventions, sources
`sourceme.sh`, and prints the complete generated script before an optional
submission. Its single-point defaults are the accepted emax2 production window
`smax=25000`, `checkpoint_s=40`, `rtol=1e-6`, and `atol=1e-8`; the flow stops
early when its generator target is reached.  The longer upper bound leaves
room for the successive `2e-4`-scale occupation modes measured in the He4
correlated-reference case; it is only a bound, not a request to integrate past
convergence.
The accepted-step safety cap is 8000; neither upper bound changes the ODE when
the generator criterion is reached earlier. Use explicit
smaller values only for a labeled diagnostic. Reinvoke it for one different
nucleus or tolerance at a time;
there is no scan or batch-generation mode.
Use an explicit single-run label such as `--label smax5` when preserving a
short diagnostic and generating a longer flow with otherwise identical
tolerances; an existing materialized output is still never overwritten.

Reference metadata records the original development-machine interaction path.
When the exact fixed file is relocated, pass `--interaction <path>`; the runner
does not trust the path alone and still requires the frozen SHA-256 before it
starts a flow.  The point7 generator uses the task-specific copy under
`/tns/mengziyan/mr-imsrg-inputs/` because the shared force tree has larger
N2LOopt spaces but not this exact emax2 file.
The production path is compatible with point7's system Python 3.9 and the
generated batch script initializes Environment Modules explicitly before it
sources the machine-aware repository environment.
Compute nodes do not carry the login node's `/usr/local` NumPy/SciPy install,
so this standalone prototype calls the shared
`/opt/library/miniconda-3.12.9/bin/python3` interpreter explicitly and logs its
versions.  It does not activate or load miniconda and does not import the
system-Python pyIMSRG module.
