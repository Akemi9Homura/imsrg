# Repository Notes

## `fmt2=no2bpack`

This repository can read the packed binary NO2B interaction format produced by the `normal-order` code.

Use it by setting:

```bash
fmt2=no2bpack 2bme=<normal-order-output.bin> 3bme=none
```

The reader is `ReadWrite::Read_no2bpack()`. It reads the binary layout written by `normal-order`'s `Write_minipack()` implementation:

- oscillator frequency and `emax`
- orbit table in `(n, l, 2j, 2tz)`
- zero-body term
- upper-triangular one-body matrix elements
- packed J-coupled two-body matrix elements
- optional center-of-mass TBME payloads, which are consumed and ignored

The reader remaps file orbit indices into the active `ModelSpace` using `(n, l, 2j, 2tz)`, so it does not depend on matching raw orbit numbering between executables.

Because `no2bpack` files already contain the normal-ordered Hamiltonian pieces from `normal-order`, the main IMSRG driver does not add another `Trel_Op` for this format.

## wm2 run conventions

Use `source ./sourceme.sh` before building or running on wm2. This checkout is built with the executable and shared library under `build/`, so generated Slurm scripts should call:

```bash
/lustre/home/2401110128/imsrg/build/imsrg++
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/lustre/home/2401110128/imsrg/build"
```

On wm2, use the Slurm partition names exactly as:

```bash
#SBATCH --partition=C064M1024G   # 64 CPU cores, about 1 TB memory per node
#SBATCH --partition=C064M0256G   # 64 CPU cores, about 256 GB memory per node
```

For the production runs here, use `#SBATCH --qos=low` and request all 64 cores on one node:

```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
```

Do not include `%N` in new log names unless the node name is intentionally needed. Prefer:

```bash
#SBATCH -o <result-dir>/log_<prefix>_%j.txt
```

## `gen_job.py` interaction parsing

Keep the result path naming independent of the input file format. The top-level result directory should be the interaction name, for example:

```text
result/<interaction>/<valence>_<reference>_hw<hw>_emax<emax>_e3max<e3max>/
```

Do not add format labels such as `no2bpack`, `no2b`, or the reference nucleus to the interaction flag unless the physics interaction name itself includes them.

For `fmt2=me2j`, keep the original filename parsing rule:

```python
_RE_E2MAX = re.compile(r'_emax(\d+)_e2max(\d+)\.')
emax_nn, e2max_nn = extract_emax_e2max(params["2bme"])
params["file2e1max"] = emax_nn
params["file2e2max"] = e2max_nn
```

For a separate 3BME file, keep the original 3B filename parsing rule:

```python
_RE_E3MAX = re.compile(r'_emax(\d+)_e2max(\d+)_e3max(\d+)\.')
emax_3n, e2max_3n, e3max_3n = extract_emax_e2max_e3max(params["3bme"])
params["file3e1max"] = emax_3n
params["file3e2max"] = e2max_3n
params["file3e3max"] = e3max_3n
```

For `fmt2=no2bpack`, use a separate parser. A packed file has `hw`, `emax`, and usually `e3max`, but it does not have `e2max`:

```python
_RE_NO2BPACK = re.compile(r'_hw(\d+)_emax(\d+)_e3max(\d+)\.')
hw, emax_nn, e3max_nn = extract_no2bpack_hw_emax_e3max(params["2bme"])
params["hw"] = hw
```

Use the parsed `hw`, `emax`, and `e3max` only for script naming and consistency checks. Do not run the normal me2j `_emax..._e2max...` parser on no2bpack files, and do not invent `file2e1max/file2e2max` values for this format.

For `fmt2=no2bpack`, set:

```python
params["3bme"] = "none"
```

Do not add `file2e1max`, `file2e2max`, `file3e1max`, `file3e2max`, or `file3e3max` to the generated command for this format. Those limits are meaningful for ordinary 2B/3B input readers, not for `Read_no2bpack()`.
