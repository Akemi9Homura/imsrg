# C++ J-scheme MR-IMSRG 仓库说明

## 最高优先级强制规则

1. **及时 commit，不能丢代码。** 每完成一个边界清楚、已经验证的阶段就立即提交；提交前检查 `git status` 和 `git diff`，只提交本任务相关文件，不夹带工作树中原有的其他改动。不得让已完成的实现长期只留在未提交工作树中。
2. **遇到问题先自行调研并作出决定。** 在用户已经授权的任务范围内，主动上网查询官方文档、开源代码、论文及其补充材料；论文必须尽量下载并阅读原文 PDF，不能只看网页摘要。结合代码和证据选择方案并继续推进，不因一般技术选择或实现偏好停下来询问用户。只有缺少必要权限、任务范围发生实质变化或操作具有不可恢复风险时才请求用户决定。
3. **尽量复用已有代码，不重造轮子。** 开始实现前先查找并理解仓库中已有的程序、接口、脚本和测试路径；已有能力可以直接完成任务时必须复用。只有现有实现确实不能满足任务需求时才新增代码，并把新增范围控制在必要的最小边界内，避免复制默认程序已经负责的初始化、运行或验证流程。
4. **仔细核对 `refs/` 中的代码与公式，正确性优先。** 实现 MR 正规序、RDM/cumulant、MR-IMSRG(2) 对易子、生成元或真空正规序转换之前，必须阅读 `refs/` 中对应的论文 PDF、公式笔记和 QCombo 代码/notebook，核清指标顺序、反对称化、组合系数、相位及采用的密度截断。不能凭记忆或只依赖 PDF 文本提取/OCR 抄写长公式；关键公式必须使用 QCombo、显式小 Fock-space 矩阵或另一份独立实现做符号或数值交叉验证，未通过单参考极限、Hermiticity/anti-Hermiticity、RDM 恒等式和随机小模型测试的实现不得用于四核结果。
5. **reader、NCSM 和下游接口优先复用家目录中的已验证实现。** 涉及相互作用 reader、m-scheme/J-scheme 映射、NCSM/壳模型基与波函数、RDM/观测量或 FCIQMC 输入输出时，必须先阅读并优先复用 `/home/mengziyan/BigstickPublick`、`/home/mengziyan/simple-ncsm`、`/home/mengziyan/shell-model-obs` 和 `/home/mengziyan/fciqmc` 中的实现、格式说明与测试；进入这些仓库前先阅读其 `CLAUDE.md`/`AGENTS.md`。MR-IMSRG 专用的 dense m-scheme readback 放在 `simple-ncsm`；小空间必须让 `simple-ncsm` 与 BIGSTICK 直接读取同一 Hamiltonian，并用完整谱、行列式集合和波函数残差交叉验收。大空间正式 NCSM 谱仍优先使用 BIGSTICK，不得用 `simpleFCI` 单独替代。`shell-model-obs` 保持上游项目范围，不再承载本项目定制修改。只有确认这些实现和本仓库现有代码均不能满足需求后，才允许新增最小替代实现。
6. **本机资源不足时转到远程 Slurm 集群。** 远程地址固定为 `mengziyan@162.105.151.37`。如果预计内存、CPU 时间或并行规模超出本开发机可承受范围，不得通过削弱物理空间或跳过验收来规避；应连接远程集群，在对应 checkout 中先执行 `source ./sourceme.sh`，再按仓库的 Slurm 约定生成、检查并提交作业。禁止在登录节点直接运行重计算；提交后必须检查队列和日志，确认输入、资源、输出及实际启动状态。只同步任务必要的代码、输入和小型结果，不提交或误同步大波函数、Hamiltonian、RDM 与日志。
7. **论文公式与现有 `imsrg++` 代码必须双重对照。** MR 新增的 cumulant 收缩要以论文和 QCombo 为依据；与 SR-IMSRG 共通的张量约定、对易子极限、生成元、能量分母、对角/非对角划分、ODE 控制、对称性维护及输出流程，必须同时阅读并优先复用本仓库当前 `src/` 中的成熟实现和测试，不能只看文章另写一套。若原代码与论文写法不同，先查明是基底、指标、截断或数值实现差异，再选择并记录与本原型一致的做法。

## 当前唯一目标

当前目标是在现有 `imsrg++` 生产路径中实现高效、球形的 C++
J-scheme MR-IMSRG(2)，供 IM-NCSM 与后续 NCSM/FCIQMC 生产计算使用。
已验收的 `prototype/mrimsrg/` Python m-scheme 实现不再是交付主路，
而是相关参考态的独立物理 oracle；现有 C++ SR-IMSRG(2) 则是单
Slater 极限的生产 oracle。

执行计划以 `docs/MR-IMSRG-C++-J-scheme实现计划.md` 与根目录
`TODO.md` 为准。禁止为了通过验收另写一套测试专用的 MR 或 SR
commutator/generator/flow；验收必须运行将来真正用于生产计算的同一
C++ 入口。

## 固定物理范围

- 相互作用：NNLOopt（本机文件名采用 `N2LO_opt`）。
- `hw = 20 MeV`。
- 单粒子截断 `emax = 2`，二体截断 `e2max = 4`。
- 只用 NN；不加入显式 3N，不做初始 3N MR-NO2B。
- 目标核：`He4`、`Be8`、`C12`、`O16`，均取最低 `J=0` 正宇称态。
- 第一批统一使用 `Nrefmax=0`：`Be8`、`C12` 是非平凡多参考测试；`He4`、`O16` 是单参考极限对照，不能宣称为非平凡 MR。
- 第一批管线通过后，仅把 `He4`、`O16` 升级到 `Nrefmax=2`，获得相关参考态版本。

固定核力文件：

```text
/home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack
```

该文件 SHA-256：

```text
76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84
```

`/home/mengziyan/ptqmc/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack` 是字节相同的副本。普通 NN `minipack` 与 `normal-order` 产生的 `no2bpack` 是不同格式，禁止混用 reader 或语义。

基准采用 HO 基、`BetaCM=0`、无库伦、无额外核子质量修正，和现有 N2LO_opt emax2 的 NCSM/FCI 基准保持一致。读取 `minipack` 时仍必须按目标质量数 `A` 构造内禀动能，因此四个核的 Hamiltonian 是四份 A-dependent 输出，不能互相复用。

## 当前算法边界

- 生产实现是 `Jref=0` 标量参考态下的球形 J-scheme
  MR-IMSRG(2)；不通过展开完整 m-scheme 张量实现生产 commutator。
- 必须复用 `Operator`/`TwoBodyME`/`ModelSpace` 的 J-coupled 块、现有
  SR contractions、`Generator` 和 `IMSRGSolver`；MR 代码只新增不可约化为
  现有 SR 项的 reference-density/cumulant 数据与收缩。
- Python m-scheme 原型仅作独立 oracle 和小模型诊断，不能成为 C++
  生产 driver 的运行时后端，也不能在 C++ 内部无条件展开整个
  m-scheme 空间。
- 参考态由 `Nrefmax` 截断 NCSM 波函数产生。
- 至少输入 `gamma1`、`gamma2`，构造并真正保留 `lambda2`。
- 当前版本设 `lambda3=0`；这一近似必须写入输出元数据。
- Hamiltonian 和生成元始终截断到参考态正规序 0B/1B/2B。
- 先使用现有直接积分 `dH/ds = [eta,H]_(0,1,2B)` 打通生产验收；
  Magnus/BCH 在直接流通过后复用现有框架，不与首个 MR commutator
  同时开发。
- 使用 IM-NCSM 的放松解耦：只处理 `Delta e != 0` 的 1p1h 与 2p2h 通道，保留同 HO 量子数参考空间内部耦合。
- 首选一个有文献公式和 QCombo 输出可核对的生成元；第一批只实现一种。不要同时开发多种生成元。
- 当前不做显式 3N、`lambda3`、奇核、`Jref!=0` 的非标量密度、
  一般张量观测算符或 MPI 重构。这里不禁止 Vobig Sec. 6.5.3
  要求的球形自然基：相关参考态的流方程、生成元和 `Delta e`
  掩码必须在该基中自洽求值，最终 Hamiltonian 再变回下游约定。

## 输入输出硬要求

内部中间格式至少保存：

- 完整轨道表与指标顺序；
- `A,Z,hw,emax,e2max,Nrefmax`；
- NCSM 参考态标识以及 `gamma1/gamma2`；
- 初始与演化后的 `E,f,Gamma`；
- 生成元、流参数、ODE 容差、`lambda3=0` 和解耦掩码；
- 核力路径及 SHA-256。

交给 NCSM/FCIQMC 的文件必须转换回普通真空正规序的

```text
E0 + one-body + two-body
```

不能直接把相对于相关参考态正规序的 `E,f,Gamma` 当成普通矩阵元输出。零体常数必须保留并由下游计入总能量。

优先输出现有下游程序可读的普通 `minipack`；如第一版先用自描述 NPZ/HDF5，则必须同时提供到 NCSM reader 的转换器，不能只停在内部张量文件。

## 开发顺序

1. 锁定已通过的 Python m-scheme 与 C++ SR 基准、轨道/相位/归一化约定。
2. 建立球形 `gamma1/lambda2` J-coupled 输入与恒等式/J 重构测试。
3. 在现有 `Commutator` 中复用全部 SR 项，只新增并逐项验证
   `lambda2` 收缩。
4. 在现有 `Generator` 和 `IMSRGSolver` 中加入 White-NCSM 和
   `Delta e != 0` 路径，不复制 ODE 框架。
5. 先让 `He4/O16, Nrefmax=0` 的同一生产入口直接退化到 SR，
   再让 `Be8/C12, Nrefmax=0` 及 `He4/O16, Nrefmax=2` 逐项复现
   Python m-scheme 原型。
6. 复用现有真空正规序、`no2bpack` 和 NCSM 读回路径完成谱验收。
7. 代数、短流、完整流和性能门禁通过后，才扩展到更大空间。

## 最低验收门槛

- `Tr(gamma1)=A`，RDM 的 Hermiticity、反对称性和收缩关系通过。
- 单 Slater 参考态给出 `lambda2=0`，并复现现有 SR-IMSRG(2) 极限。
- 单 Slater 极限不允许调用测试专用分支；生产 MR commutator 在
  `lambda2=0` 时必须调用同一批现有 SR contractions，并对
  `E/f/Gamma -> denominator -> eta -> RHS -> flow` 逐层通过。
- 相关参考态的 C++ J-scheme 结果必须展开成 m-scheme，与已验收
  Python 原型的正规序、生成元、每个 commutator 项、RHS、短流和
  最终真空 Hamiltonian 逐项比较；只比较最终能量不算通过。
- `E=<Psi_ref|H|Psi_ref>`；真空正规序与 MR 正规序往返误差不高于 `1e-10` 相对量级。
- 每步保持 `H` Hermitian、`eta` anti-Hermitian 和二体反对称性。
- `s=0` 导出的 Hamiltonian 由下游读回后复现原始谱。
- 只在 `Delta e != 0` 掩码内的解耦残差显著下降；目标相对初值至少 `1e-6`。
- ODE 容差缩小十倍后，测试能量变化小于 `1 keV`。
- 不以零体项 `E(s)` 代替后 NCSM 对角化结果。
- 已知 He4 基准：上述设置下全 emax2 FCI 基态为 `-20.33883250 MeV`；首个读写闭环必须复现它。

## 仓库工作约定

- 新的项目设计文档放在 `docs/`；上游 Doxygen 文件留在 `doc/`。
- `refs/` 是本地文献库，已被 `.gitignore` 排除，不要提交 PDF、提取文本或 QCombo 克隆。
- 构建前执行 `source ./sourceme.sh`。
- 生产 MR 必须扩展现有 `imsrg++` 路径，但默认 SR/VS-IMSRG 的物理
  与参数语义必须保持不变；新 MR 行为只能由显式参数启用。
- 快速 Python MR-IMSRG 原型的集群作业必须由专用的 `prototype/mrimsrg/gen_flow_job.py` 单参数生成器产生；同样禁止手写 Slurm 或 `sbatch --wrap`。生成后必须按本文的 point7 规则检查完整脚本并通过 `sbatch --test-only` 再提交。
- 不要提交生成的 Hamiltonian、波函数、RDM、大日志或结果目录；提交小型测试 fixture 时必须说明来源与校验和。
- 工作区可能已有用户修改。只提交当前任务明确生成或修改的文件，不顺手清理、移动或覆盖其他改动。

## `fmt2=no2bpack` 输入约定

本仓库能够读取 `normal-order` 程序产生的 packed binary NO2B interaction：

```bash
fmt2=no2bpack 2bme=<normal-order-output.bin> 3bme=none
```

reader 是 `ReadWrite::Read_no2bpack()`，其二进制布局与 `normal-order` 的
`Write_minipack()` 完全对应，依次包含：

- oscillator frequency 和 `emax`；
- `(n,l,2j,2tz)` 轨道表；
- zero-body term；
- upper-triangular one-body matrix elements；
- packed J-coupled two-body matrix elements；
- 可选的 center-of-mass TBME payload（reader 消费但忽略）。

reader 通过 `(n,l,2j,2tz)` 把文件轨道映射到当前 `ModelSpace`，不能依赖两个程序的原始轨道编号相同。因为这种文件已经包含 `normal-order` 产生的正规序 Hamiltonian，主 IMSRG driver 对 `fmt2=no2bpack` 不得再次加入 `Trel_Op`。

不得混淆以下两条路径：

- `fmt2=no2bpack`：`2bme` 是 `normal-order` 产生的 packed binary normal-ordered Hamiltonian，必须使用 `3bme=none`。
- `3bme_type=no2b`：读取普通二体相互作用（例如 `fmt2=me2j`），并通过 `3bme=<...me3j.gz>` 读取单独的三体 NO2B 文件。

质量修正和其他 bare-Hamiltonian 附加项属于普通 `fmt2=me2j` + `3bme_type=no2b` 路径。直接 `me2j` 计算为和 `normal-order` 产生的 packed 文件保持一致，通常应设置：

```python
params["nucleon_mass_correction"] = "true"
```

`fmt2=no2bpack` 的文件中这些效应应已在生成阶段加入；driver 对该格式忽略 `nucleon_mass_correction=true`，不得重复修正。

## `gen_job.py` 输入解析与命名

结果路径命名必须独立于输入文件格式。顶层目录使用物理 interaction 名称：

```text
result/<interaction>/<valence>_<reference>_hw<hw>_emax<emax>_e3max<e3max>/
```

不得仅因文件格式而把 `no2bpack`、`no2b` 或参考核名称拼入 interaction flag，除非它本来就是物理相互作用名称的一部分。

对 `fmt2=me2j` 保留原有文件名解析：

```python
_RE_E2MAX = re.compile(r'_emax(\d+)_e2max(\d+)\.')
emax_nn, e2max_nn = extract_emax_e2max(params["2bme"])
params["file2e1max"] = emax_nn
params["file2e2max"] = e2max_nn
```

单独的 3BME 文件保留原有解析：

```python
_RE_E3MAX = re.compile(r'_emax(\d+)_e2max(\d+)_e3max(\d+)\.')
emax_3n, e2max_3n, e3max_3n = extract_emax_e2max_e3max(params["3bme"])
params["file3e1max"] = emax_3n
params["file3e2max"] = e2max_3n
params["file3e3max"] = e3max_3n
```

`fmt2=no2bpack` 必须使用独立 parser。packed 文件名包含 `hw`、`emax`，通常还包含 `e3max`，但不包含 `e2max`：

```python
_RE_NO2BPACK = re.compile(r'_hw(\d+)_emax(\d+)_e3max(\d+)\.')
hw, emax_nn, e3max_nn = extract_no2bpack_hw_emax_e3max(params["2bme"])
params["hw"] = hw
```

解析出的 `hw/emax/e3max` 仅用于脚本命名和一致性检查。不得对 no2bpack 文件调用普通 me2j 的 `_emax..._e2max...` parser，也不得为该格式虚构 `file2e1max/file2e2max`。对 `fmt2=no2bpack` 必须设置 `params["3bme"] = "none"`，生成的命令中不得添加 `file2e1max`、`file2e2max`、`file3e1max`、`file3e2max` 或 `file3e3max`。

## Model-space 选择

在 `gen_job.py` 中，`params["valence_space"]` 是结果路径和输出 prefix 使用的短名称。`p-shell`、`sd-shell` 等内置 IMSRG 名称也可直接交给 `imsrg++` 解析。

非内置或裁剪空间必须用 `params["custom_valence_space"]` 写出真正的物理定义：

```python
params["valence_space"] = "<short-name-for-output>"
params["custom_valence_space"] = "<core>,<orbit>,<orbit>,..."
```

例如以 `He4` 为 core，包含质子和中子的 `p`、`d5/2`、`s1/2` 轨道：

```python
params["valence_space"] = "pd5s1-shell"
params["custom_valence_space"] = "He4,p0p3,n0p3,p0p1,n0p1,p0d5,n0d5,p1s1,n1s1"
```

设置 `custom_valence_space` 后，它才是 `imsrg++` 使用的实际 `ModelSpace`，原 `valence_space` 仅是名称。生成或提交生产作业之前，必须用物理知识核对短名称和定义是否描述同一空间；若不一致，必须向用户说明疑点并给出建议的名称或轨道表，且必须获得用户明确确认后才能运行或提交。

使用内置空间时，注释而不是删除 custom definition。例如 `fp-shell` 是内置空间（`Ca40` core，质子和中子均含 `f7/2,p3/2,f5/2,p1/2`）：

```python
params["valence_space"] = "fp-shell"
# params["custom_valence_space"] = "Ca40,p0f7,n0f7,p1p3,n1p3,p0f5,n0f5,p1p1,n1p1"
```

## wm2 运行约定

在 wm2 构建或运行前必须执行 `source ./sourceme.sh`。该 checkout 的 executable 和 shared library 在 `build/`，生成的 Slurm 脚本必须使用：

```bash
/lustre/home/2401110128/imsrg/build/imsrg++
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/lustre/home/2401110128/imsrg/build"
```

wm2 相互作用文件位于 `/lustre/home/2401110128/Forces`；point7 对应共享目录为 `/tns/public/Forces`。从 point7 检查 force 时用后者；生成 wm2 生产脚本时必须用前者。

wm2 partition 名称必须精确使用：

```bash
#SBATCH --partition=C064M1024G   # 64 cores, ~1 TB/node
#SBATCH --partition=C064M0256G   # 64 cores, ~256 GB/node
```

这里的生产计算使用 `#SBATCH --qos=low`，单节点、单 task、64 cores：

```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
```

新日志名除非确实需要 node name，否则不得包含 `%N`，优先采用：

```bash
#SBATCH -o <result-dir>/log_<prefix>_%j.txt
```

生产计算必须由 `gen_job.py` 生成 Slurm 脚本；不得手写独立脚本或用临时 `sbatch --wrap`。生成后、提交前必须检查 interaction path、`fmt2/fmt3`、truncations、`reference`、`valence_space`、operators、`BetaCM`、输出名、partition、QOS 和线程数。

`gen_job.py` 必须保持一次只生成一组参数。除非用户明确要求重新设计 batch generator，否则禁止在其中加入 beta、frequency、truncation、reference、model space、interaction 或其他扫描维度的循环/列表/批量生成逻辑。多个参数点必须手动设置一个值、生成/检查/提交，然后再改下一个单值。

## point7 运行约定

point7 checkout 位于 `/tns/mengziyan/imsrg`，executable 和 shared library 位于 `/tns/mengziyan/imsrg/build`；生成脚本的 executable、`cd` 和 `LD_LIBRARY_PATH` 必须使用这些路径。相互作用文件位于 `/tns/public/Forces`。

point7 partitions：

```text
c128m1024   128 CPU cores, ~1 TB/node（默认，常 busy/alloc）
c128m512    128 CPU cores, ~512 GB/node
compute_C    96 CPU cores/node
compute_A    28 CPU cores/node
```

在 point7 运行时，`gen_job.py` 模块级 `partition` 设置为对应名称（例如 `c128m512`）。其余 Slurm header 与 wm2 相同：一节点、一 task、`--cpus-per-task` cores、`OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-<cpus>}` 和 `-o <result-dir>/log_<prefix>_%j.txt`。

point7 使用 `accounting_storage/none` 且 `AccountingStorageEnforce=none`，不存在 slurmdbd/QOS object；`--qos=<name>` 会被静默忽略。因此 point7 脚本仍保留 `#SBATCH --qos=low`，它在那里无作用但不会报错，也便于和 wm2 共用设置。

### point7 调度与内存风险

point7 使用 `select/cons_tres`、`CR_CORE`、`OverSubscribe=NO`，CPU core 独占。`c128m512` 是双 AMD EPYC 7763，Slurm 计 128 个物理 core；两个 64-core job 可在同一节点运行且 CPU 不争用。

但是 point7 **完全不跟踪、不限制内存**：`DefMemPerNode=UNLIMITED`、`MaxMemPerNode=UNLIMITED`、`ConstrainRAMSpace=no`。`--mem` 不会保护作业，节点甚至在重负载时仍报告 `AllocMem=0`。必须据峰值 RSS 而不是当前 RSS 规划；IMSRG Magnus 的 `Omega` 会随 flow 累积，`emax=14,e3max=24` 的 fp-shell 运行峰值远超 80 GB。共址作业总 RSS 超过节点约 515 GB 时，Linux OOM killer 会任意杀进程，Slurm 不提供保护。

需要独占节点时，在 `gen_job.py` 使用单值模块级 `nodelist`（例如 `"node2"`）写入 `#SBATCH --nodelist=<node>`；设为 `None`/`""` 则由 Slurm 选择。`nodelist` 不能改成扫描列表。多个作业必须逐个改一个 node、生成、检查、提交。共址前也可用 `scontrol show node <node> | grep FreeMem` 检查，但要记住该值不等于 Slurm 的内存约束。

### machine-aware 环境

`sourceme.sh` 必须根据机器分支，并在每次构建或运行前 source：

- wm2：`cmake/3.31.9`、`OpenBLAS/0.3.17`、`gsl/2.7.0`、`boost/1.83.0`，并设置 wm2 GSL 所需的 `GSL_ROOT_DIR/CMAKE_PREFIX_PATH`。
- point7：`cmake/3.25.2`、`openblas/0.3.10-single`、`gsl/2.7.1`、`boost/1.81.0`。

point7 的 `imsrg++` 链接 `libopenblas.so.0`（来自 `openblas/0.3.10-single`）和 `libgsl.so.27`（来自 `gsl/2.7.1`）。若作业立即报 `ERROR: Unable to locate a modulefile for ...`，说明 wm2 module name 泄漏到 point7 分支，`set -e` 会使脚本在启动 `imsrg++` 前退出。source 后用 `ldd build/imsrg++ | grep "not found"` 检查运行库。

默认不得加载 `miniconda`：pyIMSRG 针对 system Python 构建，加载 miniconda 会改变 `python3` 并破坏已有 module import。

### 生成、提交和检查

```bash
python gen_job.py --generate-only   # 生成脚本并打印路径；提交前检查
python gen_job.py --submit          # 生成并通过 sbatch 提交
python gen_job.py --smoke-test      # 运行 imsrg++ help，轻量检查构建/运行环境
python gen_job.py                   # 生成后询问 submit job: y/n；n 表示本地 bash 运行
```

提交前检查 interaction paths、`fmt2/fmt3`、truncations、`reference`、`valence_space`、operators、`BetaCM`、`denominator_delta`、输出名、partition、QOS 和 thread count。提交后必须检查 `squeue -u $USER` 以及结果目录中的 `log_<prefix>_<jobid>.txt`，确认已越过 module loading 并开始读取 interaction。

对单一 major shell 且没有 cross-shell coupling 的 valence-space target，可以关闭 intruder suppression 和 CM terms：

```python
params["BetaCM"] = 0.0
params["denominator_delta"] = 0.0
```
