# C++ J-scheme MR-IMSRG(2) 生产实现计划

## 1. 目标与非目标

当前任务是在现有 `imsrg++` 中实现 `Jref=0` 球形参考态的
C++ J-scheme MR-IMSRG(2)。它必须是真正的生产路径，而不是调用
Python 原型或为测试复制一套 SR 公式的 wrapper。

两个独立 oracle 同时约束它：

1. `lambda2=0` 的单 Slater 极限必须直接调用并逐项复现现有
   C++ SR-IMSRG(2)。
2. `lambda2!=0` 的相关参考态必须展开成 m-scheme，逐项复现已
   通过 QCombo、显式 Fock-space 和 NCSM 验收的 Python 原型。

当前范围不包含显式 3N、`lambda3`、`Jref!=0` 非标量密度、奇核、
一般张量观测算符和 MPI 重构。首个生产闭环仍使用 NNLOopt、
`hw=20 MeV`、`emax=2`、`e2max=4`、NN-only 和已验收的六个参考态。

## 2. 冻结的物理定义

- 参考态是最低 `J=0+` NCSM 态，流在球形自然轨道基中进行。
- `gamma1` 是 rank-0 一体密度，在自然基中由每个 J-orbit 共同的
  magnetic-substate 占据 `n_p` 表示。
- `lambda2=gamma2-A(gamma1 gamma1)` 是 rank-0 反对称二体 cumulant；
  生产存储必须与 `Operator::TwoBody` 的 normalized pair/TBME 约定一致。
- `lambda3=0`，Hamiltonian 和生成元截断为参考态正规序 0B/1B/2B。
- 生成元是 Vobig Sec. 6.5.4 的 `White-NCSM`：分子和 EN 分母舍去
  `lambda[2,3,...]` 项，但 MR commutator 保留全部线性 `lambda2` 项。
- 仅解耦 `Delta e != 0` 的 1p1h/2p2h 通道，保留同 HO 量子数
  参考空间内部耦合。
- direct flow 已验收并冻结为 oracle；当前门禁是复用现有生产 Magnus/BCH
  框架实现 scalar MR-Magnus(2)。

## 3. 实现前必须完成的公式审计

代码开始前要建立一张可追溯的 contraction 表，每行至少记录：

| 输出 rank | m-scheme 来源 | J-scheme 归约 | 相位/归一化 | C++ 入口 | 独立验收 |
|---|---|---|---|---|---|
| 0B | Hergert 2016 Eq. 49 / QCombo | `1/4 C2·lambda2` 的 J trace | identical-pair normalization | MR 0B addon | J-to-m random oracle |
| 1B | Hergert 2016 Eq. 50 / QCombo | 四类 2B--2B--`lambda2` 收缩 | 6j/phase/channel | MR 1B addon | named Python terms |
| 2B | Hergert 2016 Eq. 51 | 现有 fractional-occupation SR J contractions | 现有约定 | 直接复用 | `lambda2=0` 与 MR random |

审读材料至少包括：

- Hergert 2016 MR-IMSRG 论文 PDF 及 TeX 源，尤其 MR-IMSRG(2)
  flow Eqs. 49--51 和 generalized Wick 附录；
- Gebrerufael 2017、Vobig 2020、Mongelli 2022 学位论文的 IM-NCSM
  程序、球形密度、生成元和放松解耦章节；
- MR-IMSRG/IM-NCSM 系列论文与 `refs/qcombo/examples/MR_IMSRG2.ipynb`；
- `src/Operator.cc`、`TwoBodyME.cc`、`TwoBodyChannel.cc`、`ModelSpace.cc`、
  `Commutator.cc`、`Generator.cc`、`IMSRGSolver.cc`、`UnitTest.cc` 和
  `ReferenceImplementations.cc`。

不允许从 OCR 文本直接抄长公式。J-scheme 归约要由 Clebsch--Gordan
展开数值重建，再与 QCombo/m-scheme 逐元素核对。

审计结果已冻结在 `docs/MR-IMSRG-J-scheme公式与代码映射.md`。特别是
Gebrerufael 未正规化 pair 矩阵元、`TwoBodyME` normalized-pair 存储和
rank-0 TBTD 三种数值不能混用。

## 4. 软件边界

### 4.1 MR reference 对象

新增一个小型、显式传递的 `Jref=0` reference 对象，包含：

- ModelSpace/轨道顺序 fingerprint；
- `A,Z,Nrefmax` 和参考态/RDM 校验和；
- 每个 orbit 的 `n_p` 与 `1-n_p`；
- J-coupled `lambda2` 块；
- `lambda3=0`、natural-basis 和 `Jref=0` 标记。

密度不得依赖未记录的全局单例，也不得通过改写普通 SR
occupation 来丢失 `lambda2`。若复用 `Operator::TwoBody` 存储 `lambda2`，
必须由类型/接口明确它是密度而非 Hamiltonian。

### 4.2 MR commutator

生产入口的结构为：

```text
existing SR/fractional-occupation J-scheme contractions
  + MR 1B lambda2 contraction
  + MR 0B contraction of the completed 2B commutator with lambda2
  = MR-IMSRG(2), lambda3=0
```

2B flow 方程不显式含 `lambda2`，应完全复用现有 J-scheme SR 收缩。
0B 的 `C2·lambda2` 必须收缩已经完整累加的 2B commutator，不能只收缩
2B--2B 子项。1B 新项保持独立命名，以便与 Python/QCombo 的四类
`lambda2` 收缩逐项比较。

### 4.3 generator/solver/driver

- `Generator` 新增显式 `white_ncsm` 模式，不改变现有 generator 名称语义。
- 分母复用 `GetTBMEmonopole()` 的 `(2J+1)` 加权 monopole 和当前
  `1e-6 MeV` cutoff 行为。当前 SR 代码对小分母一律置为正 cutoff；首版
  MR 为严格退化必须原样复用，不在 MR 分支单独改变符号规则。
- `Delta e` 掩码由 orbit `e=2n+l` 和 pair total `e` 决定，不用
  core/valence/q-space 标签猜测。
- `IMSRGSolver` 继续负责 ODE、对称性维护、checkpoint 和停止；MR 参考
  只通过显式上下文传给 commutator/generator。
- `imsrg++` driver 新增显式 MR 开关和 reference-density 输入；未启用时的
  SR/VS-IMSRG 输入、数值和性能必须保持不变。

## 5. 验收梯子

### Gate 1：密度与 J-coupling

- RDM 恒等式在 J-scheme 中直接通过；
- 随机 scalar `lambda2` 的 J-to-m-to-J 误差 `<=1e-12`；
- 真实参考态另报告输入 m-scheme RDM 的标量投影残差，拒绝阈值先取
  `1e-10`；coupled block 自身闭环仍要求 `<=1e-12`。

### Gate 2：命名 contraction

- 随机标量 Hamiltonian/生成元/`lambda2` 下，每个新 J-scheme MR 项
  展开到 m-scheme 后 max-abs `<=1e-10 MeV`。
- `lambda2=0` 时新项为精确零，生产 dispatcher 的 SR 结果不变。
- 对易子满足预期的 Hermitian/anti-Hermitian 类型和 pair 反对称性。

### Gate 3：真实体系 `s=0`

- `He4/O16, Nrefmax=0` 对当前 SR `imsrg++` 比较
  `E/f/Gamma -> denominator -> eta -> each RHS term -> total RHS`。
- `Be8/C12, Nrefmax=0` 和 `He4/O16, Nrefmax=2` 对 Python m-scheme
  比较同一整条链。
- 代数级 max-abs 门槛 `<=1e-10 MeV`，J/m 转换本身 `<=1e-12`。

### Gate 4：共同短流

- 共同 RHS 的 `ds=1e-4` Euler 一步；
- 共同固定步 RK4 的 `s=0.001,0.002,0.003` checkpoints；
- 每点比较完整 `H/eta/RHS`，不只看 0B 或能量。

### Gate 5：完整流、物化与谱

- 同一生成元/cutoff/停止条件下比较完整 flow；
- 容差收紧十倍后后 NCSM 能量改变 `<1 keV`；
- 真空 `E0+t+V` 的 float64 J-coupled readback 先通过，再做 float32
  `no2bpack` 格式验收；
- NCSM 读回逐态比较，不以 MR 0B 代替后对角化能量。

### Gate 6：生产性能

- RHS 不保存完整 m-scheme 张量；
- 主循环只遍历对称允许的 orbit/J/channel 块；
- 记录新 MR 收缩的 profiler 名称、峰值 RSS 和 wall time；
- 在至少 `emax=2` 六个已验收参考态上，C++ 的时间和内存显著优于
  Python m-scheme，且并行执行不改变数值。

### Gate 7：真实较大空间运行

- 先把已经验收的 `Nrefmax` 截断参考态嵌入更大的单粒子空间：旧轨道上的
  占据、自然变换和 `lambda2` 必须逐元素保持；新增轨道取零占据、零
  cumulant 和单位自然变换。该操作只扩大 Hamiltonian/流空间，不冒充在
  更大参考空间重新对角化得到的新 NCSM 波函数。
- 首个真实门禁固定为本机已有且来源明确的 NNLOopt
  `hw=20, emax=4, e2max=8` 和 `He4 Nrefmax=2` 嵌入参考。至少执行完整
  RHS 与固定短流，并记录峰值 RSS、wall time、profiler 分项、Hermiticity、
  anti-Hermiticity、线程复现误差和 checkpoint 大小。
- emax4 运行过程中不得创建 m-scheme 四指标张量；若任一 scratch 的实测
  峰值偏离解析上界，先定位具体 contraction，再改成 channel/block 流式
  构造。
- emax6 只在相互作用来源、header 与 SHA-256 验证后进入验收。文件名含
  `candidate` 的输入只能用于诊断，不能写进正式结果表。
- `emax=14` 的解析存储表仍是容量规划，不是已运行证据；计划和验收文档
  必须始终分别报告“解析外推”与“实际运行”的最大空间。

emax6 输入前置现已建立。冻结的 emax14/e2max28 母 minipack SHA-256 为
`118c4c27ab7e2b1f3ae8df8ea24beb2784df681dd286bdb3ae5f42dbe09b9d7d`；
流式保留原始 float32 interaction/Hcm/`p1.p2` payload 的 child SHA-256 为
`199d5a7a3e060427582d112bd2ab36688db199629ddebb6e097a4676980b913b`。
extractor 先用 emax4→emax2 正式文件完成逐字节校准，再由
`Hamiltonian::truncate(6)` 在 A=4 和 A=16 下独立读回，0B/1B/2B 均严格
零差。旧 `candidate` 在两种 A 下同一 2B 元均差 `3.814697265625e-6 MeV`，
只判为额外 float32 roundtrip 诊断品，不能作为正式输入。下一步可在不改变
参考态定义的前提下进入真实 emax6 RHS/短流。

该 emax6 实跑现已完成。He4 Nrefmax=2 固定参考从 emax2 嵌入到 56 个
J-orbit，`lambda2` 收缩误差仍为 `1.85e-15`；canonical minipack 经
A=4 reader 转为 467032 条 TBME 的 lossless J64 并零误差读回。固定
`ds=smax=1e-4` RK4 一步单线程墙钟 `11.70 s`、峰值 RSS `281132 KiB`，
流后 `E=-17.016801664428 MeV`、`||H||=10983.182992216 MeV`，真空导出零体
一致性严格为零。profiler 中 V/VI 为 `7.751/0.231 s`，V 又分为 build
`5.855 s`、BLAS `1.873 s`。1/64 线程 0B/1B 相同、2B 最坏
`8.88e-16 MeV`；64 线程在 10-core 开发机上用时 `14.20 s`，不宣称并行
加速。Gate 7 因此已覆盖真实 emax4 和 emax6 RHS/短流，但仍不等于 emax6
完整收敛流。

根据 emax6 profiler，下一轮只对 V 项固定不变的 `lambda2` 支撑做精确
降维。定义活跃轨道集合为所有非零标准耦合 cumulant 元素实际出现的轨道；
Pandya 变换后的 `Lambda[(a,b),(c,d)]` 若任一轨道不在该集合中必为零，故

```text
X Lambda Y = X(:,A) Lambda(A,A) Y(A,:)
```

是恒等变换而不是 density threshold。一般全空间相关参考态令 `A` 等于完整
有序 pair 集，自动恢复原实现；固定 emax2 参考嵌入 emax4/6 时才获得降维。
实现 commit `cef4fce6` 已通过随机全活跃参考、SR 精确退化、emax2 Python
m-scheme (`3.106e-12 MeV`) 及 emax4/6 短流。两份短流 J64 均与优化前
逐字节相同；emax4 墙钟 `0.602 -> 0.408 s`，emax6
`11.70 -> 6.81 s`，后者 V/build/BLAS 分别为
`2.851/2.765/0.072 s`，峰值 RSS `280880 KiB`。emax6 1/64 线程的
0B/1B 差为零、2B 最坏 `8.88e-16 MeV`。随后 point7 generator jobs
`100385/100387/100389/100391` 对该 commit 重跑四核
`s=100, rtol=atol=1e-10`：C++/Python 全流最坏差依次为
`2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`。相对前一版 Pandya 复用，
流后 J64 的全局 max-abs `8.384e-13 MeV`；现有 `simpleFCI` 按
`Nmax=8/0/0/2` 重算 12 条低能谱，最大变化 `1.741e-13 MeV`。因此该
活跃支撑优化已通过完整流和 NCSM 门禁；下一步才是根据 emax6 外推规划
更大空间的完整流、Magnus 或重启策略。

对 `cef4fce6` 的 profiler 再检查发现，`Lambda(A,A)` 已经很小，但 X/Y
仍先按完整 `d*d` Pandya 块构造再取活跃行列。commit `dea9125b` 把 partial
support 分支改为直接形成 `X/Y(:,A)` 与 `X/Y(A,:)`；只有 `A=d` 的一般
全空间相关参考才继续一次性形成完整 X/Y，因而不会用四个 skinny 矩阵把
全活跃工作翻倍。emax4/emax6 短流与上一版 J64 逐字节相同，墙钟分别
`0.408 -> 0.312 s` 和 `6.81 -> 4.15 s`；emax6 V/build 降为
`0.170/0.111 s`，相对上一版快 `16.76/25.00` 倍。随机全活跃参考和四核
checkpoint 的长 CTest 在 `352.52 s` 后通过。此时 emax6 的
MR `lambda2` addon 为 `0.717 s`，既有 White-NCSM generator update 为
`2.191 s`，主瓶颈已经跨到公共生成元路径。该结论在 point7 对
`dea9125b` 的四核 full-flow/NCSM 回归后仍成立。专用生成器提交的 jobs
`100393/100395/100397/100399` 均完成 `s=100, rtol=atol=1e-10` 检查；
He4/Be8/C12/O16 的 C++/Python 全流最坏差依次为
`2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`。新 J64 相对
`cef4fce6` 的全局 max-abs 为 `8.527e-13 MeV`，由现有 `simpleFCI`
按 `Nmax=8/0/0/2` 重算的 12 条低能谱最大变化为 `1.990e-13 MeV`。
因此 skinny Pandya 构造已通过完整流、真空物化和后 NCSM 门禁。下一步
不再继续改 MR V 项，而是先剖析和复用现有 White-NCSM generator 路径，
再根据 emax4/6 实测外推决定大空间完整流及 Magnus/重启策略。

这一步现已在 commit `6e02676c` 完成。Vobig Eqs. (6.5.31--34) 只依赖
随当前 Hamiltonian 更新、但在一次生成元构造内反复复用的对角 monopole
`Gamma_ijij`；实现每个 RHS 先构造一次完整有序轨道对表，再让原分母公式
按原算术次序查表，不改变 `Delta e` 掩码、占据权重、cutoff 或分子。
emax4/emax6 短流相对 `dea9125b` 均逐字节相同，generator update 从
`0.1112/2.191 s` 降到 `0.00324/0.0426 s`，总墙钟从 `0.312/4.15 s`
降到 `0.222/1.98 s`。随机全活跃参考长门用时 `354.53 s` 并通过。
point7 jobs `100401/100403/100405/100407` 又完成四核全流；C++/Python
最坏差保持为 `2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`，新旧 J64
全局 max-abs `9.308e-13 MeV`，12 条 NCSM 低能谱最大变化
`1.137e-13 MeV`。该瓶颈已消除；emax6 当前主要耗时转为 MR addon
`0.725 s` 和现有 SR commutator `0.611 s`，足以进入大空间完整流与
积分/Magnus/重启策略的实测门禁。

### 5.1 当前下一阶段：大空间 direct flow 与可恢复运行

下一阶段不直接跳到 Magnus。先把现有、已经在 emax2 验收的 direct ODE
路径变成一条可复现的单参数生产作业，并用 emax4 实测回答是否真的需要
Magnus。执行顺序固定为：

1. 在当前环境重新执行 QCombo `MR_IMSRG2.ipynb`，把输出的 0B/1B/2B
   contractions 与生产 `MRCommutator.cc`、慢速 J-scheme reference 和
   Python m-scheme oracle 三方核对。该步骤只重生独立证据，不修改已经
   通过数值门禁的公式来迎合 notebook 输出。
2. 根目录 `gen_job.py` 增加独立 `--mr-jscheme` 模式；每次只接受一个核、
   一个 `Nrefmax`、一个 `emax` 和一个累计目标 `s`。脚本必须记录 git
   commit、输入/reference SHA-256、求解器及 ODE 参数来源、point7 资源
   以及两种物化输出，且拒绝覆盖现有结果目录。生产入口固定使用原有
   `method=flow`，不另设或暴露 `ds_0/dsmax/ode_tolerance`，直接继承
   `imsrg++` 运行时默认值；固定步 `flow_RK4` 只留在底层测试入口。
3. direct flow 的恢复语义冻结为 Hamiltonian 分段：上一段导出的 vacuum
   J64 作为下一段 `fmt2=jcoupled64` 输入；同一 reference 重新执行
   HO→NAT→MR 正规序后积分 `target_s-start_s`。`start_s/target_s` 只记录
   累计物理流参数，solver 收到的 `smax` 是本段长度。先在 emax2，再在
   emax4 比较 `0→S` 与 `0→S1→S` 的完整 0B/1B/2B，而不是只比能量。
4. 通过连续/分段等价后，由 `gen_job.py` 在 point7 提交
   `He4 Nrefmax=2, emax=4` 完整 direct flow。提交前执行
   `sbatch --test-only` 并逐项检查 interaction/reference、commit、
   `fmt2/3bme`、截断、生成元、`BetaCM`、ODE 容差、partition、线程和
   输出路径；提交后检查队列及日志已经越过 module/reader 初始化。
5. emax4 终点必须同时通过 lossless J64 和 float32 no2bpack 的独立 NCSM
   读回，并重复一次十倍更严 ODE 容差。只有收敛残差、谱稳定性、wall、
   峰值 RSS 和 checkpoint 体积都有实测记录，才外推 emax6 完整流。
6. 若 emax6 direct flow 可由分段 checkpoint 安全完成，则暂不实现 Magnus；
   若 wall/checkpoint 证明确有必要，再把 MR Magnus 当成新门禁，从同一
   RHS、短流、真空物化和后 NCSM 谱逐层验收，不能只比较零体项。

这条路径不改变既有 SR/VS `gen_job.py` 默认参数，也不把 J64 checkpoint
当成严格的 ODE 内部状态恢复：它是从物化后的 Hamiltonian 开始下一段自主
流。连续/分段差异必须作为 ODE 误差单独量化并随结果保存。

前 3 步现已完成。commit `fee05bdf` 的 `gen_job.py --mr-jscheme` 只生成
一个参数点，冻结 commit、executable/input/reference SHA-256、累计流参数、
point7 环境和 J64/no2bpack 输出，并由 CTest 与真实 emax4 脚本预演覆盖。
随后以 `rtol=atol=1e-10`、`eta_criterion=1e-12` 比较
`0→0.02` 和 `0→0.01→0.02`：emax2 的 vacuum 0/1/2B 最坏差为
`3.29e-14/2.13e-14/6.15e-14 MeV`，emax4 为
`2.31e-14/2.13e-14/5.92e-14 MeV`；终点 MR `E`、`||H||` 和
`||eta1,2||` 在 flowfile 打印精度相同。因此 J64 已验收为 Hamiltonian
分段恢复点，但仍不称为 adaptive stepper 内部状态 checkpoint。完整哈希、
wall/RSS 和最坏 TBME 见 `MR-IMSRG-Jscheme-large-space.json`；下一项是
point7 的 emax4 完整收敛流。

point7 的第一轮实测也已完成，但暴露了必须先澄清的停止尺度。commit
`6bd0c200`、64 核 `compute_C` 上，job `100413` 从累计 `s=0` 积分到
`100`，墙钟 `159.38 s`、峰值 RSS `166.99 MB`；job `100423` 从其 J64
Hamiltonian checkpoint 继续到累计 `s=1000`，墙钟 `187.56 s`、峰值
RSS `167.19 MB`。终点 `Rgen=||eta||` 分别为 `1.10716e-2` 和
`1.41309e-3`，相对同一初始值 `1.07854` 为 `1.02654e-2` 和
`1.31019e-3`，因此两段虽然正常完成，却都没有通过本项目的相对 `1e-6`
严格门槛。两段 MR 反正规序后的零体核对在打印精度为零，J64 SHA-256
分别为 `c45495a0...fc1636`、`67d8f171...3c300f`，no2bpack SHA-256
分别为 `7640dd30...71f85`、`dcf8ef57...6e05`。

这不是生成元构造慢或 ODE 容差造成的假尾。逐通道工具
`prototype/mrimsrg/diagnose_white_ncsm_tail.py` 直接复用生产
`Generator::ConstructGenerator_WhiteNCSM`，并把 `s=0/100/1000` 的
21,840 个 `Delta e != 0` 1B/2B J-scheme 通道拆成未加权 Hamiltonian、
White-NCSM 分子 `Rnum`、正反方向 EN 分母和 `Rgen`。它逐项重构的
`eta` norm 与 `Operator::Norm()` 在 `1e-12` 内相同。结果显示：

- `Rnum/Rnum(0)=1, 3.47839e-3, 4.40316e-4`；
- `Rgen/Rgen(0)=1, 1.02654e-2, 1.31019e-3`；
- `s=1000` 主导通道的分母仍约 `24.1 MeV`，没有接近 `1e-6 MeV`
  cutoff 的小分母；
- 参考态的非零小占据包括 `2.03e-4, 2.20e-4, 2.53e-3, 2.59e-3,`
  `1.76e-2, 1.77e-2`。这些占据乘积直接进入 Vobig 6.5.28--34 的方向性
  分子，使快通道先消失、越来越小权重的通道随后接管 norm，形成慢包络；
- 未加占据权重的掩码 Hamiltonian norm 反而从 `1728.85` 增到
  `1844.71 MeV`，因为 White-NCSM 的已发表固定点不是把所有
  `Delta e != 0` 裸 TBME 逐元素压为零，不能用这个量替换 `Rgen`。

Vobig Fig. 6.7 的 White-NCSM 曲线只展示到 `s≈100`，其实现说明使用
二阶能量修正相对零体项约 `1e-4--1e-3` 的停止判断，并没有给出
`Rgen/Rgen0=1e-6` 的证据。因此已经对 `s=0,0.02,0.1` 做
`Nmax=2,4,6,8,16` 后 NCSM，并对 `s=1,10,100,1000` 做 `Nmax=8`。
裸 Hamiltonian 的 emax4 全空间 (`Nmax=16`) 基态为
`-25.2912232776 MeV`；`s=0.02` 把 Nmax2/4/6/8 的同终点截断误差分别改善
`412.148/147.117/39.306/9.869 keV`，但全空间谱相对裸值漂移
`-45.476 keV`，而且 `Rgen/Rgen0=0.9659`，只消除了约 3.4% 的生成元
norm。`s=0.1` 的收敛改善更大，但全空间漂移已为 `-187.528 keV`；
`s=100/1000` 的 Nmax8 基态达到 `-30.588/-42.621 MeV`，明确显示长流
IMSRG(2) 截断漂移，不能作为生产输出。

所以当前冻结两类互不替代的结论：严格 `Rgen/Rgen0<=1e-6` 脱耦仍未通过；
`s=0.02` 仅是有限短流 IM-NCSM/NCSM 收敛加速器候选。该点的 J64 与
float32 no2bpack 三态差最大 `0.000726 keV`；ODE 容差从 `1e-9` 收紧到
`1e-10` 后 J64 0/1/2B 最坏变化 `3.675e-14 MeV`，NCSM 能级变化约
`3e-11 keV`，格式和积分门禁通过。下一步不是继续盲目追 `eta`，也不是
假设 Magnus 会消除物理截断误差，而是先冻结 finite-s downstream stopping
rule：同时约束最大可算 Nmax/全空间谱漂移、Nmax 收敛改善以及 ODE/格式
误差。只有这条规则写清后才把同一短流窗扩到 emax6；严格门槛继续单列为
失败，若两种目标在 MR-IMSRG(2) 下不兼容就如实报告不兼容。

该 operational rule 现冻结为 `prototype_downstream_stability_v1`：在已测点
中选择通过全部条件的最大 `s`，要求全空间基态相对裸 Hamiltonian 漂移
`<=100 keV`、每个已测 Nmax 相对同一流后全空间的截断误差均改善或不变、
J64/no2bpack packing 与十倍 ODE 容差误差各 `<1 keV`，且生成元 norm
必须下降。100 keV 是本快速原型在进入 emax6 前显式选定的工程预算，不是
Vobig 给出的阈值；文献依据仅是 Sec. 6.5.5 用后 NCSM 平台区分有效预演化
与诱导高体项失控。`select_downstream_flow_window.py` 从机器 JSON 重算
每项判据，当前选择 `s=0.02`，而严格脱耦仍独立返回失败。

### 5.2 当前执行边界：emax6 finite-s pilot

emax4 已证明长流会放大 IMSRG(2) 截断误差，因此 emax6 不再先追逐
`eta/eta(0)<=1e-6`。执行对象冻结为同一个 He4 `Nrefmax=2` emax2 NCSM
参考态的严格 embedding、canonical NNLOopt `emax=6,e2max=12` 输入和
`s=0.02` 短流。该 embedding 只是隔离 C++ J-scheme 容量与有限流行为，
不是在 emax6 参考空间重新求解 NCSM。

当前步骤和判定顺序如下：

1. 输入必须追溯到已冻结的 emax14 母文件并逐 channel 校验；当前使用文件
   SHA-256 为 `199d5a7a...b913b`。旧 candidate 与母文件相差
   `3.8147e-6 MeV`，不得用于结果。
2. point7 只通过根目录 `gen_job.py --mr-jscheme` 运行 direct flow。
   job `100432` 已完成 `s=0->0.02`（`10.48 s`, `886.6 MB`）；终点
   `||eta||/||eta(0)||=0.96099`，所以只进入 downstream pilot，不称为
   严格脱耦。job `100434` 已完成十倍更严 ODE 容差复算。
3. J64 与两条独立 no2bpack 写路径必须由同一 NCSM reader 对角化。
   当前 `Nmax=8` 三态谱在两条 no2b 写路径间为 double 数值噪声；
   `1e-9` 对 `1e-10` 的 J64 0/1/2B 最坏差为
   `1.18e-13/1.13e-13/4.16e-13 MeV`，谱差约 `1.2e-10 keV`。
4. 裸与流后 Hamiltonian 使用完全相同的 NCSM 截断序列。`Nmax=2,4,6,8`
   的相邻 gap 均得到改善；Nmax8 流/裸基态差为 `-102.717 keV`，因此继续
   到 Nmax10。jobs `100441/100442` 用带输入、executable 和环境哈希的
   `--mr-ncsm-readback` 单点脚本运行；计算节点不提供 `/usr/bin/time`，
   生成器使用 `/proc/<pid>/status` 采样 VmHWM。
5. 两个 Nmax10 job 均以 exit `0:0` 完成，维数 `183866`；裸/流后基态为
   `-26.790750086/-26.860681522 MeV`，`Nmax8->10` gap 从
   `618.748` 降为 `585.962 keV`，Nmax10 流/裸差为 `-69.931 keV`。
   所以“emax6 最大可算 Nmax 代理”通过。由于全空间是 Nmax24，这不等于
   emax6 完整 `prototype_downstream_stability_v1` 通过，也不提供全空间
   谱漂移结论。
6. 上述物理/格式/ODE 门均已固化，下一步才评估 emax8。Magnus 不是谱漂移
   修正项；是否实现只由 direct-flow 内存、checkpoint 和 wall 实测决定，
   并仍需重复 RHS、物化与后 NCSM 验收。

### 5.3 当前执行边界：emax8 实测 finite-s pilot

emax8 输入与资源前置已经完成。canonical NNLOopt `emax=8,e2max=16`
minipack 从同一冻结 emax14 母文件直接抽取，SHA-256 为
`3fd1a0038e8c6f00cd93ec26113f02fcb19dd51268a95238da67b17e12fd4bda`；
独立 reader 在 A=4/16 下逐 channel 均为严格零差。bare J64 含 90 个
J-orbit、3526624 条 TBME，SHA-256 为 `cd14cd4b...eb9cce`；嵌入参考 SHA-256
为 `9af7a6e0...a68b3`，保持原 emax2 `Nrefmax=2` 参考的占据和 cumulant，
新增块为空。

point7 job `100444` 已完成 `ds=smax=1e-4` 固定 RK4 一步。4 次完整 RHS
墙钟 `17.28220 s`、峰值 RSS `2070.297 MB`，初末 MR 零体项为
`-17.0135029683/-17.0170919053 MeV`，`||eta2||` 从 `1.214269979` 降至
`1.214014842`；真空导出零体一致性为零。该实测证明 pure J-scheme direct
flow 可在 emax8 运行。相反，660 个 m-orbit 的稠密四指标 double 张量需要
`1517978880000 bytes`，所以会 `bad_alloc` 的旧 lossless J64→dense-m
validator 被明确排除；它不是生产 J-scheme 流失败，也不应成为 emax8 格式
验收的依赖。

隔离的 `simple-ncsm` native no2bpack reader 的 Nmax8 jobs
`100447/100448` 均完成，维数
44838。bare 与固定一步流的三态谱分别为
`[-26.2170529924,-24.7096273755,-24.7096273017]` 和
`[-26.2177863230,-24.7099878250,-24.7099679677] MeV`，峰值 RSS 均约
2.28 GiB。基态微移 `-0.733331 keV` 与正常读回只构成轻量格式/资源 smoke
test；该简化求解器不再作为正式 NCSM 谱验收。

这不再意味着 `simple-ncsm` 的正确性问题被搁置。其 emax8/Nmax2
sanitizer 故障已定位为用单个 `uint64_t` 编码 90 个 J-orbit occupation
导致的移位未定义行为，现改为 collision-free occupation vector key。用原
emax8、660 m-orbit、59 维输入复跑后，三态逐位复现且进程正常 `exit 0`。
此外新增 `ncsm_bigstick_no2_compare`：在 He4 N2LOopt
`hw20/emax2/e2max4/Nmax2` 上让两个程序直接读取同一 SHA-256 为
`5cadb860...e01de` 的 no2bpack，比较全部 59 个行列式、59 条谱和每个
BIGSTICK 波函数残差。最大谱差/最大残差为
`4.165e-6/4.137e-6 MeV`，通过 `1e-5 MeV` 门禁。差值来自 BIGSTICK
version-1 `.wfn` 的 float32 能量/振幅；用该振幅在 simple-ncsm double
Hamiltonian 上重算的 Rayleigh 商相对 exact 谱最坏仅
`3.340e-7 MeV`。因此小空间两个 NCSM 实现已严格对齐；大空间 finite-s
正式谱仍必须由 BIGSTICK 完成，不能以 simple-ncsm 单边结果替代。

下一步固定为：

1. 用 `gen_job.py --mr-jscheme` 单点运行 `s=0→0.02`、
   `rtol=atol=1e-9`，仍用同一 White-NCSM RHS、J64 checkpoint 和
   no2bpack 物化；已由 job `100451` 完成；
2. 用 `1e-10` 独立复算并比较 lossless 0/1/2B；job `100452` 已通过，
   Nmax8 jobs `100458/100456` 的 `<1 keV` 谱一致性只保留为 lightweight
   smoke test，不替代下一项 BIGSTICK 门禁；
3. 由 BIGSTICK 对 bare/flow 使用完全相同的 Nmax 序列，从 2/4/6/8
   逐点增加到实测
   资源允许的最大值；在没有 emax8 全空间（四粒子 Nmax32）时，只能按预先
   声明的最大可行 Nmax proxy 报告，不能声称完整
   `prototype_downstream_stability_v1` 通过；
4. 只有 direct flow 的内存、wall 或 checkpoint 实测成为阻碍时才实现
   Magnus；不得用 Magnus 解释或掩盖 IMSRG(2) 谱漂移。

两条完整流的 25/31 次 scalar commutator 分别用时 `67.67063/85.03911 s`，
峰值 RSS `5194.719/5227.641 MB`；终点 `E=-17.70238931537 MeV`、
`||eta||/||eta(0)||=0.959464`，所以仍只称 finite-s point。双容差 J64 的
0B/1B/2B 最坏差为 `4.80e-14/6.53e-14/1.90e-13 MeV`。native no2bpack
Nmax8 三态在两套容差间最大差 `4.97e-11 keV`。松容差第一次读回 job
`100454` 有完整谱但子进程返回非零，正式重跑 job `100458` 以
`COMPLETED 0:0` 得到同一谱；验收只采用后者。现在唯一未完成的是步骤 3 的
匹配 Nmax 序列与最大可行空间判定。

这一阶段继续保留两个独立结论：C++ J-scheme 公式与实现门禁已经通过；
有限流能否作为更大空间下游 Hamiltonian，则由后 NCSM 数据单独决定。

## 6. 实现顺序与 commit 边界

1. 文献/C++ 调研表与 J-coupled density 约定；
2. MR reference 读入及 J/m 密度闭环；
3. 单个命名 `lambda2` contraction 及其随机 oracle；
4. 完整 MR commutator 与 SR 退化；
5. White-NCSM generator 与 `Delta e` 掩码；
6. solver/driver 集成和真实体系 `s=0`/short-flow 验收；
7. 完整流、物化、NCSM 与性能验收。
8. reference embedding、真实 emax4/6/8 RHS/短流和 profiler 驱动优化。
9. emax8 `s=0.02` 双容差流、真空物化及最大可行 Nmax 后 NCSM pilot。
10. QCombo 再审计、单点 Slurm 生成器、direct-flow 分段等价和 emax4
   完整流/NCSM；依据实测再决定 emax6 direct flow 或 MR Magnus。

每个边界在自身测试通过后立即独立 commit。任何一个未通过 Gate 1--4
的中间实现都不能用于新物理结果。

当前进度（2026-08-25）：步骤 1 完成；步骤 2 的显式 `MRReference`、
J-coupled 恒等式检查、正规序往返及 `mrimsrg_jref_v1` 文件 reader 已完成；步骤 3--4
的慢速参考、块矩阵生产收缩、完整 0B/1B/2B MR commutator 及精确 SR
退化门禁已完成。步骤 5 的 White-NCSM 方向分子、MR EN 分母和
`Delta e != 0` 掩码也已完成并对 Python m-scheme 验证。步骤 6 的 solver
显式 MR context、直接 RK4/ODE RHS、`imsrg++` driver 与 He4 真实
`s=0`/short-flow 首项验收均已接通。HF/HFMBPT 原有 normalized-pair
变换抽成 `Operator::TransformOneAndTwoBody()`：它验证正交性与球形 block，
拒绝 3B/非粒子数守恒输入，并由 HF、HFMBPT 与 MR 共用。真实
`He4 Nrefmax=2` 自然轨道下，随机标量 NN/NO2B 算符的 HO→NAT→HO 总误差
`<2e-11`，NAT 一体和完整 m-scheme 二体张量分别以 `2e-11`、`3e-11`
绝对容差通过独立四指标协变对照。当前不得把 He4 首项通过外推成全部体系；
步骤 6 剩余边界是 He4/O16 SR 退化与 Be8/C12/O16 的同层真实门禁。

严格门禁不能以 float32 `no2bpack` 作为内部真值，因此步骤 6 同时复用
原型验收的 `mrimsrg_j64_v1`，在 `ReadWrite` 增加 float64 scalar J-coupled
reader/writer。它按 `(n,l,2j,2tz)` 映射轨道并验证完整 channel/pair 覆盖；
随机 Hermitian NN/NO2B 写读误差 `<1e-13`，读取既有真实
`He4_s0_srcheck.jcoupled64` 与 Python reader 构造的 `Operator` 逐元素误差
为零。该格式只作无损生产输入和验收检查点，最终 NCSM/FCIQMC 文件仍由
既有 `no2bpack` writer 产生。

`imsrg++` 的显式 `mr_reference_file=<...jref>` driver 现已接通。入口强制
`basis=oscillator`、单步 direct flow、`white-ncsm`、NN/NO2B、`BetaCM=0`、
无额外 flowing operators，并拒绝 Magnus、3B、二次 model-space truncation
与 valence-space 重正规序。执行序列为

```text
vacuum H(HO) -> H(NAT) -> MR normal order -> direct flow
             -> undo MR normal order -> vacuum H(NAT) -> vacuum H(HO).
```

真实 NNLOopt `He4 Nrefmax=2` 的自然基 `E/f/Gamma` 对 Python m-scheme
max-abs 为 `2.95e-13/3.36e-11/9.95e-12 MeV`。端到端 driver 在 `s=0`
和 `ds=1e-4` 单步 RK4 后的 float64 0/1/2B 输出，均与直接调用同一 C++
reference/generator/solver 路径逐元素为零；`s=0` 对原 HO vacuum Hamiltonian
的最坏往返误差为 `1.15e-14/1.07e-14/8.88e-15 MeV`。回归测试同时确认
两点都能生成非空下游 `no2bpack`。当前步骤 6 剩余工作转为 He4/O16 SR
退化与 Be8/C12/O16 真实参考的共同 checkpoints，而不再是 driver 接线。

随后 `He4/O16, Nrefmax=0` 的真实 driver SR 门禁也已通过：MR 入口使用
`white-ncsm` 和显式零 `lambda2` reference，原生入口使用现有 `white` 且
不设置任何 MR context；二者在 `s=0` 及共同 `ds=1e-4` RK4 后导出的
vacuum 0/1/2B float64 矩阵元全部逐元素为零。这一回归通过 reference 文件
内容判断零 cumulant，不含 He4/O16 核名分支。下一层仍需单独报告真实
denominator/eta/named RHS 的误差表，并扩展到相关 Be8/C12/O16 reference。

相关 reference 的同层真实门禁现已扩展到 `He4/O16 Nrefmax=2` 与
`Be8/C12 Nrefmax=0`。相对 Python m-scheme，四体系自然基 `E/f/Gamma`、
`eta(s=0)`、总 `RHS(s=0)` 的全局 max-abs 上界依次为
`6.51e-11/9.98e-13/7.06e-11 MeV`，均通过 `1e-10 MeV` 门槛。每个体系的
主 driver 在 `s=0` 和 `ds=1e-4` RK4 后又分别与直接 C++ 路径逐元素为零。
这些真实 RHS 随后已进一步拆成八个既有 SR contraction 及
`mr_lambda2_one_body/zero_body`。把原始 m→J 输入投影误差与代数误差分离后，
同一 J 可表示输入上的逐命名全局 max-abs 为 `1.954e-14 MeV`、相对
Frobenius 为 `1.380e-15`；最坏索引已持久化。共同
`s=0.001,0.002,0.003` checkpoints 的最坏 H/RHS 为
`1.279e-13/2.665e-14 MeV`，四核 2832 条 EN 分母最坏
`1.421e-14 MeV`。当前 Gate 5 已增加相关参考固定-s Python/C++ 验收器，
本地 He4 Nrefmax=2 `s=0.01` 烟雾流的最终 H/eta/RHS/vacuum 全局最坏
`4.24e-11 MeV`。point7 上统一 `rtol=atol=1e-9`、`s=100` 的正式比较现已
完成：`Be8/C12 Nrefmax=0` 与 `He4/O16 Nrefmax=2` 的终点
`H/eta/RHS/vacuum H` 全局 max-abs 分别为
`4.40e-9/6.88e-8/1.87e-8/1.04e-8 MeV`，机器记录见
`docs/MR-IMSRG-Jscheme-full-flow.json`。流后的 float64 J64 已按显式轨道
量子数和 `(a,b,c,d,J)` 记录映射到独立 `shell-model-obs` channel 顺序；
四核 NCSM 均完成读回。再转成下游 float32 `no2bpack` 后，基态相对 J64
的最大差为 `1.225e-6 MeV`（`0.001225 keV`）。因此 Gate 5 的相关参考
完整流、真空物化和谱闭环已经通过。SR 逐项证据也已汇总为
`docs/MR-IMSRG-SR-degeneration.json`：He4/O16 的输入、分母、eta、十个
命名 RHS 项、总 RHS、Euler/RK4 和 `s=100` 完整流全部通过，代数层最坏
为 `4.26e-14 MeV`（分母）和 `1.15e-14 MeV`（命名 RHS）。ODE 容差从
`1e-9` 收紧到 `1e-10` 后，四核各三条最低 NCSM 能级的全局最坏变化为
`1.683e-11 MeV = 1.683e-8 keV`，远低于 `1 keV`。因此 Gate 5 全部通过；
逐核 Slurm job、两套谱和差值已追加到
`docs/MR-IMSRG-Jscheme-full-flow.json`。

Gate 6 的首个可复现基准也已完成。真实 `He4 Nrefmax=2`、单线程、三次
中位数下，C++ J-scheme MR RHS 为 `6.88e-3 s`，Python dense m-scheme 为
`1.766 s`，快 `256.8` 倍；两个输入算符加参考密度的主数值存储分别为
`41,872` 与 `61,478,400 bytes`，相差 `1468` 倍。只构造 ModelSpace 并
解析 channel 尺寸的外推表显示，`emax=14` 一个 J-scheme 0/1/2B 算符约
`3.94 GB`，相应稠密 m-scheme 二体张量约 `438 TB`。完整命令和数据见
`prototype/mrimsrg/benchmark_jscheme.py` 与
`docs/MR-IMSRG-Jscheme-performance.json`；这证明存储结构能外推，不等于
已经完成 `emax=14` 的物理流。

生产 `lambda2` 的 ordered particle-hole 收缩原先会同时保留所有
J/parity/Tz 块的 `X/Y/Lambda` scratch；现改为逐块构造、BLAS 收缩后立即
释放。解析峰值估算在 `emax=14` 由 `48.46 GB` 降为 `2.87 GB`，而真实
He4(ref2) RHS 与 Python oracle 的 max-abs 仍为 `3.11e-12 MeV`。该数字只
计算这项主要 cross-block scratch，不冒充整个进程的峰值 RSS。

Gate 7 的 reference embedding 已实现为生产 `MRReference` 操作，而不是
重建 m-scheme RDM。它按 `(n,l,j2,tz2)` 映射旧轨道，保持旧占据、NAT 和
normalized-pair `lambda2` 逐元素不变；新增轨道严格取零占据、零 cumulant
与单位 NAT 块。`MRReference::WriteBinary()` 写出全部目标 J-channel 并拒绝
覆盖文件，随后由新 ModelSpace 独立读回。真实 `He4 Nrefmax=2` 已从
emax2 嵌入本机 NNLOopt `emax=4,e2max=8`，目标 interaction SHA-256 为
`d3dff5faa2a58d8c234914170caffa3649d0cd3805a1344818f4f4d3c37fd19e`；
读回的 cumulant 收缩误差 `1.85e-15`、Hermiticity 误差为零。这里保留原始
RDM/波函数校验和以标识固定参考态来源，不把 embedding 记作较大参考空间的
新 NCSM 解。

该 emax4 实测也已完成。普通 minipack 先通过不构造 m-scheme 张量的
`mrimsrg_minipack_to_j64` 转为 A-dependent float64 J64，并逐矩阵元零误差
读回；生产 driver 随后对 30 个 J-orbit、34320 条 TBME 做一个固定
`ds=1e-4` RK4 步。单线程 wall time `1.37 s`、峰值 RSS `54,360 KiB`；
1/64 线程输出的 0B/1B 差为零，2B 最坏 `4.44e-16 MeV`。新 profiler
显示 MR `lambda2` 的 V/VI 收缩分别用时 `0.628/0.492 s`，是下一轮明确的
优化对象。机器记录见 `docs/MR-IMSRG-Jscheme-large-space.json`；这仍是
RHS/短流门禁，不冒充 emax4 完整收敛流。

VI 的首轮机械优化已完成：把与外部 `(one,two,J1)` 无关的 `(J2,w)`
有限迹预计算为 `ly_trace/lx_trace`，不改任何 channel、相位或系数。emax4
VI 从 `0.492 s` 降至 `0.022 s`，单线程短流从 `1.37 s` 降至 `0.88 s`；
优化前后 J64 逐元素为零，1/64 线程最坏仍为 `4.44e-16 MeV`。emax2 对
Python oracle 的 RHS 误差保持 `3.11e-12 MeV`，完整相关参考 CTest 用时
`350.07 s` 并通过。当前真实主瓶颈已收敛到 V 的 ordered Pandya 块。

V 的下一轮优化也保持公式和 contraction 次序不变：一次 Pandya 求和同时
构造 `X/Y/lambda2`，Six-J 改用 `imsrg++` 已预计算的 dense cache，并用
现有双算符 TBME accessor 合并 `X/Y` 的 channel、ket 与 phase 查询。
emax4 的 V 从 `0.616 s` 降至 `0.332 s`，单线程短流由 `0.88 s` 降至
`0.602 s`；相对最初 `1.37 s` 快 `2.28` 倍。优化前后单线程 J64 bitwise
相同，1/64 线程两体最坏差仍为 `4.44e-16 MeV`；emax2 Python oracle 为
`3.11e-12 MeV`，完整相关参考 CTest 用时 `351.93 s` 并通过。单块内部
OpenMP 在 emax4 未产生实测加速且增加线程驻留，试验改动已撤回。Gate 7
的 profiler 优化仍须完成 point7 完整流和 NCSM 回归后才最终勾选。

该远端回归现已完成。point7 用生成器提交 He4/Be8/C12/O16 jobs
`100377/100379/100381/100383`，固定旧 Python `s=100, rtol=atol=1e-10`
轨迹并以优化后的 C++ 重新积分。四核 C++/Python 全流最坏差依次为
`2.212e-8/3.898e-8/8.870e-10/1.071e-9 MeV`，和优化前一致；新旧流后
真空 Hamiltonian 的全局 max-abs 为 `7.44e-13 MeV`。下载新 J64 后由
现有 `simpleFCI` 按 `Nmax=8/0/0/2` 重新对角化，12 条低能谱相对优化前
全局最坏变化 `2.13e-13 MeV`。因此 Gate 7 的 profiler 优化子项通过；
emax6 仍必须等可验证 interaction，不能由本结果外推冒充实跑。

当前 MR-P2 已继续推进到标准耦合 IV/VI 中间量的显式内存审计。原实现即使
`lambda2` 只支撑固定 emax2 参考空间，也会为每个 channel 保存完整
`LY=lambda2*Y` 和 `LX=lambda2*X`；commit `b0c832a8` 改为只保存 cumulant
精确 pair 支撑对应的活跃行，并通过一个保持原有交换相位、channel 选择和
normalized-pair `sqrt(2)` 因子的 accessor 为 VI 提供矩阵元。IV 所需的
`XLY/YLX` 仍完整形成；全空间活跃参考仍走原 `lambda*Y/X` 的 BLAS 分支，
因此没有加入阈值、低秩近似或新的 MR 截断。profiler 现在直接报告本组
standard products 的显式存储，以及相对完整 `LY/LX` 避免的字节。

对同一个 He4 Nrefmax=2 嵌入参考，emax4/6/8 分别避免
`907/14067/108943 KiB` 的 dense `LY/LX`，新 standard-product 存储分别为
`1198/14949/110935 KiB`。三个空间一个 `ds=1e-4` RK4 步的完整流后 J64
都与提交前冻结文件逐字节相同；emax6 的 1/16 线程输出 SHA-256 也相同。
emax8 当前 MR addon 为 `1.91721 s`，但进程峰值 RSS 仍约 `1.35 GiB`，
说明完整 0/1/2B 算符和现有 SR Pandya scratch 主导总峰值，不能把消除的
约 106 MiB MR 临时分配直接等同于 RSS 降幅。Release 下
`MRReference/MRDriver/MRSRDriver/MRDenominators` 与完整四体系
`MRCorrelatedDriver`（`354.88 s`）全部通过；ASan+UBSan 下生产 MR driver、
SR 退化、分母以及 emax4 稀疏活跃行短流也通过。下一轮只继续优化 MR
IV/VI 的精确 channel 遍历和临时量；共享 `comm222_phss` 虽已是 emax8
最大单项，但任何改动都必须另带完整 SR/VS 逐元素与性能门禁。

随后 commit `21760bbf` 只合并 MR 中成对张量的重复查询。IV 的
`XLY/YLX` 和 VI trace 的 `LY/LX` 共享一次 channel、ket、交换相位及
normalized-pair 归一化查询；VI 主和直接复用 `imsrg++` 已有的
`TwoBodyME::GetTBME_J_twoOps()` 同时读取 `X/Y`，并在调用前按精确角动量
三角条件排除零项。外部 `one/two` 只遍历现有 scalar one-body channel，
其顺序仍是原轨道序；没有重排任何有效求和，也没有改变公式或密度截断。
emax4 的 IV/VI/addon 从 `0.00860/0.01415/0.06659 s` 降至
`0.00404/0.00562/0.05402 s`，emax6 从
`0.07247/0.15042/0.46645 s` 降至
`0.02706/0.04122/0.25384 s`，emax8 从
`0.25080/0.77368/1.91721 s` 降至
`0.11566/0.19998/1.12244 s`；相邻 emax8 短流 wall 实测为
`13.86 -> 11.36 s`，峰值仍约 `1.35 GiB`。

emax4/6/8 的流后 J64 相对 `b0c832a8` 冻结输出全部 bitwise 相同，emax6
的 1/16 线程文件 SHA-256 相同。Release 的随机全活跃参考、生产 MR/SR
driver、分母和完整四体系 Python m-scheme oracle（`355.27 s`）通过；
ASan+UBSan 下生产 MR/SR driver 和 emax4 稀疏嵌入短流通过。下一项集中在
MR V 的 ordered-cross-block build：先检查 active-active Pandya 子块是否
在 row/column 两次构造中被重复计算，只有能保持全活跃分支和逐位输出时才
允许消除重复；不转向有限-s NCSM 行为研究。

实测否定了两个表面上合理但不值得保留的方案：复制 row/column 交叠的
active-active Pandya 方块对 emax8 V build 没有可测收益；块内 OpenMP 也未
缩短 16 线程 wall，反而增加约 30 MiB 峰值。两段试验代码均已撤回。随后
定位到真正的结构浪费：partial-active V 原先形成两个完整 `d*d` 的
`X*lambda*Y/Y*lambda*X`，而一体 partial trace 只读取 ordered pairs 共享
spectator `t` 的元素。commit `287adea4` 保留全活跃参考的原 BLAS 分支；
嵌入大空间分支只保留 `X*lambda` 与 `Y*lambda` 的 `d*a` 左因子，并对
实际读取的元素按相同活跃指标顺序做点积，不形成未使用的 dense 输出。

emax4/6/8 的单个最大 block 分别避免 `441/4390/26244 KiB` dense scratch。
emax6 的 V/BLAS/addon 从 `0.16431/0.04214/0.25384 s` 降至
`0.13445/0.00283/0.22717 s`；emax8 从
`0.62761/0.25488/1.12244 s` 降至
`0.34349/0.00590/0.86220 s`。emax8 单线程 wall 为
`11.36 -> 10.95 s`，16 线程为 `7.27 -> 6.95 s`；总 RSS 仍由完整算符和
共享 SR scratch 的另一时刻主导。emax4/6/8 J64 与优化前 bitwise 相同，
emax8 1/16 线程 SHA-256 相同；Release 完整四体系 oracle 用时
`352.75 s` 并通过，ASan+UBSan 的 MR/SR driver 与 emax4 partial-active
短流通过。下一项把同一“只形成实际读取的迹元素”原则用于标准耦合 IV 的
`XLY/YLX`，它们仍是当前 standard products 中最大的 MR dense 输出。

commit `460a009d` 已把该原则落实到 IV。`StandardProducts` 对空 cumulant
channel 不再分配两个全零 `d*d` 矩阵；partial-active channel 只保存
`X/Y(:,A)` 与 `lambda*Y/X(A,:)`，并对实际进入 IV partial trace 的坐标
点积；全活跃 channel 继续执行原完整 BLAS。完整矩阵、skinny IV 和 VI
活跃行现在共用同一个 pair 坐标函数，统一 channel 检查、交换相位和相同
轨道的 `sqrt(2)` normalized-pair 因子。

emax4/6/8 分别避免 `1050/14495/109900 KiB` IV dense 输出；全部 channel
同时驻留的 standard products 数值存储从 `1198/14949/110935 KiB` 降至
`291/883/1993 KiB`。emax8 setup/IV/addon 从
`0.13488/0.12549/0.86220 s` 降至
`0.04524/0.04068/0.67228 s`，单/16 线程 wall 分别为
`10.95 -> 10.80 s` 与 `6.95 -> 6.85 s`。总进程峰值仍约 `1.35 GiB`，
因为它出现在完整 Hamiltonian/SR scratch 占用阶段；这里报告的是已明确
消除的 MR 中间量，不能混写为 RSS 降幅。emax4/6/8 J64 继续 bitwise
相同，emax8 1/16 线程 SHA-256 相同；Release 完整四体系 oracle
（`361.73 s`）及 ASan+UBSan MR/SR driver、emax4 空/full/partial channel
短流均通过。当前 emax8 MR addon 只占单线程 profiler real 的 `6.3%`。
point7 没有可追溯的 canonical emax10 输入，因此没有临时裁剪或虚构一个；
容量门直接提升到现存 canonical NNLOopt `emax=12,e2max=24`，比原计划更强。

### 6.1 emax12 实测容量门

根目录生成器先增加独立、单参数 `--mr-prepare-jscheme` 模式，冻结 minipack、
源参考、转换器、`pyIMSRG.so`、嵌入器、环境和 git commit 的 SHA-256。
第一次计算节点试跑发现 point7 没有 `/usr/bin/time`，commit `6b786e4f`
改为与既有 NCSM runner 相同的 `/proc/<pid>/status` VmHWM 采样；第二次试跑
发现 embedder 只为三个格式常量意外加载 `export_jref -> sympy`，commit
`6975c2e6` 把 jref 格式常量和 `pyIMSRG` loader 抽成无符号代数依赖的小模块。
这两个失败均发生在输入准备阶段，没有产生可用物理结果，也没有靠手写脚本
绕过生成器。

修复后的 point7 job `100531` 为 `COMPLETED 0:0`。canonical minipack
SHA-256 为 `b735e7d7...e1c3b1`；A=4 转换得到 182 个 J-orbit、156 个
J-channel、73842036 条 TBME，J64 round-trip 最坏误差严格为 0。裸 J64
大小 `2067713204 bytes`、SHA-256 `0f11452d...c35553`；嵌入 jref 大小
`1181746556 bytes`、SHA-256 `f273a851...67b07`，Hermiticity/收缩误差为
`0/1.8492e-15`。转换和嵌入分别用 `42/13 s`，采样峰值
`8087548/4298452 KiB`。

第一次 flow 调度测试 job `100533` 被 Slurm 放到只有约 31 GiB 实际空闲
内存的 `master`；它在相互作用读取前被主动取消，不能算容量结果。重新由
生成器指定当时约 617 GiB 空闲的 node8，job `100535` 用 64 threads、
`flow_RK4`、`ds=smax=1e-4` 完成 4 次完整 RHS，并同时写出 J64 和
no2bpack。进程 wall 为 `185 s`，`/proc` 与内部 profiler 一致给出峰值
`26264640 KiB = 25649.0625 MiB`；Slurm 连同输入/输出哈希用时 `227 s`。
初末 MR 零体项为 `-17.0134907986/-17.017145544706 MeV`，组合
`||eta||` 比为 `0.9997883895`。该比值只确认一步流实际执行，不是脱耦
判据。反正规序后参考对角元与 MR 零体项差为 0；流后 J64/no2bpack 的
SHA-256 分别为 `3f8a9e84...13cfc` 和 `79b43f95...dd585`。

四次 RHS 合计的 MR `lambda2` addon 为 `8.69811 s`，占 profiler real
`183.76964 s` 的 `4.73%`。其中 IV/V/VI/setup 为
`0.30983/2.54807/4.13056/0.92176 s`，V build 只有 `1.13381 s`
（`0.6%`）。活跃支撑实现让 standard products 实存峰值保持
`6494 KiB`，并分别避免 `2305898/381307/2302859 KiB` 的 IV/V/LY-LX
dense 临时量。这个实测否定了当前新增 V recoupling plan/cache 的必要性：
MR 路径已不是 emax12 的容量或 wall 主瓶颈，继续复杂化不会有成比例收益。
本节只验收生产 C++ J-scheme MR 实现和物化能力，不启动 finite-s NCSM
研究，也不把一步结果称为收敛 Hamiltonian。完整机器记录见
`MR-IMSRG-Jscheme-large-space.json` schema v14。

### 6.2 sanitizer Python 异常门闭环

此前完整 sanitizer 数值路径已经通过，但 `MRReference.WriteBinary()` 的
“拒绝覆盖已有文件”测试一抛 C++ 异常就由 ASan 在
`__interceptor___cxa_throw` 内部终止。最小脚本证明 C++ 语义本身正确：
第一次写成功，第二次确实进入拒绝分支。根因是 `/usr/bin/python3` 是纯 C
host，ASan 初始化时 libstdc++ 尚未加载，所以 `real___cxa_throw` 为空；
这不是 MRReference 越界或 pybind 未转换异常。

commit `04202213` 在 CMake 的 MR Python tests 上自动识别 address-sanitized
build，保持 libasan 为 `LD_PRELOAD` 第一项并同时预加载编译器对应的
libstdc++。每个测试还把 `PYTHONPATH/LD_LIBRARY_PATH` 固定到自己的 binary
tree，防止 sanitizer build 意外导入普通 Release `pyIMSRG.so`。由于 host
CPython 未 instrument 且退出时保留全局分配，Python-hosted tests 设置
`detect_leaks=0`；这项限制只针对该 CTest 环境，AddressSanitizer 与 UBSan
仍以 `halt_on_error=1`、ASan `abort_on_error=1` 运行，任何地址或 UB 命中
都会立即失败。

修复后 ASan+UBSan `MRReference`（包括第二次 `WriteBinary` 抛出并由 Python
捕获）以 `15.51 s` 通过；`MRDriver/MRSRDriver/MRDenominators/
MRJobGenerator/MRTailDiagnostic/MRDownstreamWindow` 六项合计 `46.85 s`
全通过；最重的四体系 `MRCorrelatedDriver` 以 `465.57 s` 通过。普通 build
对应七项测试合计 `41.40 s` 通过。至此 sanitizer 环境不再是 MR-P0 的
未决项。

### 6.3 emax2 严格脱耦门禁闭环

emax4 的 `s=1000` 流只把 `Rgen/Rgen(0)` 降到 `1.31e-3`，因此此前生产
C++ J-scheme 路径仍缺一份真正达到项目 `1e-6` 相对门槛的证据。这里不再
扩展 finite-s IM-NCSM 研究，而是选择仍能由 Python m-scheme oracle 独立
检查、同时具有非零 `lambda2` 的 `He4 Nrefmax=2, emax=2/e2max=4`，专门
关闭实现验收缺口。输入固定为 NNLOopt、`hw=20 MeV`、White-NCSM/EN 分母、
`Delta e != 0` 放松掩码和 `lambda3=0`；bare J64/reference SHA-256 分别为
`55a3c161...44e2` 与 `1c9934aa...8898`。

初始独立重算为
`Rgen(0)=0.7035911887073647`，所以相对 `1e-6` 对应绝对上界
`7.0359118871e-7`。作业停止值预先取更严格的 `7.0e-7`。point7 job
`100540` 先以 `rtol=atol=1e-10` 流到 `s=1000`，得到
`Rgen/Rgen(0)=1.13293852e-3`；lossless J64 SHA-256 为
`7b02f677...44df0`。随后 jobs `100548/100550` 从同一个 J64 Hamiltonian
checkpoint 继续，除 ODE 容差 `1e-10/1e-9` 外输入、reference、executable
和物理参数完全相同。两条单线程作业分别用 `8037/8099 s`、峰值
`20348/22504 KiB`，均正常 `exit 0`。28 线程重复轨迹数值一致但在这个
极小空间更慢，单线程反超后被取消；它不计作终点证据，也说明 emax2 不应
默认用大 OpenMP 线程数。

终点 J64 由 `diagnose_white_ncsm_tail.py` 重新读取并逐通道构造生产
White-NCSM generator，而不是采用 flowfile 的九位小数：

- `rtol=1e-10`：`Rgen=6.99999986354e-7`，相对初值
  `9.94895896352e-7`；
- `rtol=1e-9`：`Rgen=6.99999983706e-7`，相对初值
  `9.94895892588e-7`；
- 两条的 `eta1≈1.26e-11`，终点最小绝对分母约 `18.929 MeV`，不存在
  denominator cutoff 伪停止；
- 物理流约在累计 `s≈4.0774e5` 过线；其后停止判据把 generator 置零，
  adaptive solver 快速前进到名义输出点而 Hamiltonian 不再变化；
- White-NCSM 未除分母的 numerator norm 相对初值为 `3.14690e-7`。

含 `lambda2` 的 `D-D†` 诊断从 `19.1285` 降到 `10.8775 MeV`，也随
机器记录保存，但它不是 Vobig White-NCSM 生成元的固定点定义，不能事后
替换预先冻结的 `Rgen` 门槛。正式结论只按上述 `Rgen/Rgen(0)`，两条均
严格通过。

数值稳定、重启与物化门禁也同时关闭：`1e-10` 对 `1e-9` 的 lossless J64
最坏 0B/1B/2B 差为 `5.80e-12/3.19e-12/4.28e-10 MeV`，`Nmax=2` 三态
最大差 `2.38e-8 keV`。生产 driver 从严格终点 J64 重启，停止判据使其执行
0 次 commutator；输出相对输入的 0B/1B/2B 最坏差为
`1.60e-14/8.88e-15/7.11e-15 MeV`，反正规序零体核对严格为零。同一个
NCSM reader 对 J64/no2bpack 的 59 维三态读回最大格式差
`9.35e-4 keV`，低于 `1 keV` 门槛。最后七项
`MRReference/MRDriver/MRSRDriver/MRCorrelatedDriver/MRDenominators/
MRJobGenerator/MRTailDiagnostic` 以 `417.82 s` 全通过，其中完整四体系
m-scheme oracle 为 `373.60 s`。

因此生产 C++ J-scheme MR-IMSRG 已获得一份非平凡相关参考态的严格脱耦、
十倍 ODE 容差、J64 restart、真空物化、SR 退化和 m-scheme oracle 的闭环
证据。该结论验收的是 MR 实现本身，不把此长流 Hamiltonian 宣称为下游
IM-NCSM 生产结果。完整机器记录见
`MR-IMSRG-Jscheme-large-space.json` schema v15。

### 6.4 完成审计与整库回归

严格门关闭后又从生产入口向下做了一次源码级审计。显式 MR 参数只在
`imsrg++.cc` 中启用 `MRReference`；`IMSRGSolver` 无参考态时仍调用原
`Commutator::Commutator`，有参考态时先执行相同 SR commutator，再由
`MRCommutator` 加入 `lambda2` 的 1B/0B 收缩。`lambda2=0` 分支在任何
MR 浮点加减之前直接返回 SR 结果。White-NCSM 继续使用现有
`Generator`/Epstein--Nesbet 分母和球形 HO `Delta e != 0` 掩码；输出先在
自然基反 MR 正规序，再转回 HO 并由同一真空算符写 J64/no2bpack。因此
生产流中不存在测试专用核分支、Python 后端或完整 m-scheme tensor。

审计同时发现 `src/CMakeLists.txt` 中六个重复上游测试把脚本路径写成
`../test/...`；从 `build/src` 执行时会错误解析到不存在的 `build/test`。
把它们改为 `${CMAKE_SOURCE_DIR}/test/...` 后，从当前 checkout 重新执行：

```bash
source ./sourceme.sh
cmake -S . -B build
ctest --test-dir build --output-on-failure
```

结果为 `20/20` 通过、总用时 `658.08 s`。其中既有上游 SR/通用测试、真实
He4/O16 单参考生产退化、四核相关参考 m-scheme oracle（`505.11 s`）、
分母、作业生成器、严格尾部诊断和下游格式窗均在同一次 CTest 中通过。
结合 P0--P4 的 QCombo/显式 Fock-space、十倍 ODE 容差、严格残差、NCSM
读回以及 emax12 容量证据，当前声明的 `Jref=0, lambda3=0, NN/NO2B,
direct-flow` 生产 C++ J-scheme MR-IMSRG 范围没有未关闭的实现或验收项。
这不扩大范围到显式 3N、`lambda3`、非标量参考、Magnus 或有限-s IM-NCSM
效果研究。

## 7. 错误定位原则

- SR 极限失配：先检查是否真正复用了现有 contraction、occupation 和
  TBME normalization，不通过加核特例修正。
- MR 失配而 SR 通过：按命名 `lambda2` 收缩逐项查相位、6j、pair
  normalization 和 J-channel 选择。
- 只有流后失配：用共同 Euler/RK4 隔离 RHS 与 ODE driver，再检查对称化和
  停止条件。
- 能量接近但矩阵元不同：仍判为失败，禁止以截断误差或基底差异
  做无证据解释。

## 8. 当前阶段：MR-Magnus(2)

### 8.1 文献冻结与数学定义

本阶段以 Morris--Parzuchowski--Bogner 2015 Eqs. (28),(30)、Vobig 2020
Secs. 4.5--4.6/App. B.2 和 Mongelli 2022 Secs. 5.4--5.6 为直接依据：

```text
dOmega/ds = sum_k B_k/k! ad_Omega^k(eta),  Omega(0)=0
H(s)       = sum_k 1/k! ad_Omega^k(H(0))
ad_Omega^0(X)=X,  ad_Omega^k(X)=[Omega,ad_Omega^(k-1)(X)]
```

`Omega` 与 `eta` 为 anti-Hermitian，`H(s)` 为 Hermitian；所有对象和每个
中间嵌套对易子均按 MR 正规序截到 0B/1B/2B，`lambda3=0`。Mongelli
Sec. 5.6 明确把 Magnus commutator 写成普通 MR commutator，因此不能给
BCH 另造近似收缩。当前不处理 rank-J 非零张量算符。

原 `method=magnus` 不直接积分上述 Bernoulli ODE，而用一阶指数 Euler：每步
以 `BCH_Product(ds*eta,Omega_old)` 合成有限变换。它适合检查 BCH-product
代数，但真实 He4 有限流表明默认 `domega=0.1` 没有达到 keV 稳定性。因此
它保留为 oracle，生产入口改用论文中的完整 Bernoulli ODE：已有
`method=magnus_adaptive` 现在最多逐阶计算到 `k=8`，并按 Vobig
Eqs. (4.6.7),(4.6.8) 在非零 Bernoulli 阶用相对 `1e-2` 项阈值和嵌套范数
单调性停止；每层 scalar commutator 均接到 MR dispatcher，并以原生
`ode_tolerance/dsmax/omega_norm_max` 做
Dormand--Prince RK45 和分段。每个 RHS 仍由
`exp(Omega) Hstart exp(-Omega)` 的 MR-BCH 重建 Hamiltonian。direct flow
保留作独立交叉验证。

#### 出处与本仓库决策的边界

上述数学与数值选择的出处如下；“论文公式”和“本仓库根据验收做出的
工程决策”不得混为一谈。

1. T. D. Morris, N. M. Parzuchowski, and S. K. Bogner,
   *Magnus expansion and in-medium similarity renormalization group*,
   Phys. Rev. C **92**, 034331 (2015),
   [doi:10.1103/PhysRevC.92.034331](https://doi.org/10.1103/PhysRevC.92.034331),
   [arXiv:1507.06725](https://arxiv.org/abs/1507.06725)。Eq. (28) 给出完整
   Bernoulli 形式的 `dOmega/ds`，Eq. (30) 给出用 `Omega` 变换
   Hamiltonian 的 BCH 级数；论文紧接 Eq. (29) 后明确说实际计算通过
   数值积分 Eq. (28) 构造 `Omega`。本地原文为
   `refs/papers/2015_Morris_Magnus_IMSRG_arXiv1507.06725.pdf`。
2. K. Vobig, *Electromagnetic Observables and Open-Shell Nuclei from the
   In-Medium No-Core Shell Model*, TU Darmstadt dissertation (published 2020),
   [URN:nbn:de:tuda-tuprints-113758](https://nbn-resolving.org/urn:nbn:de:tuda-tuprints-113758)。
   Eq. (4.5.7) 是完整 Magnus ODE，Eq. (4.5.14) 是 BCH 观测量变换；
   Eqs. (4.6.1)--(4.6.4) 定义 Magnus(2) 中对 `Omega`、导数和每个
   嵌套对易子的 NO2B 截断，Eqs. (4.6.7),(4.6.8) 给出相对
   `1e-2` 项阈值与嵌套对易子范数递减条件。Appendix B.2 明确
   说明他们在每个 `Omega` ODE 积分步用 BCH 计算 `H(s)`，并用
   GSL 的自动步长 RKF45。本地原文为
   `refs/theses/2020_Vobig_IMNCSM_PhD.pdf`。
3. T. Mongelli, *The In-Medium No-Core Shell Model as Comprehensive Ab-Initio
   Tool*, TU Darmstadt dissertation (2022),
   [URN:nbn:de:tuda-tuprints-216719](https://nbn-resolving.org/urn:nbn:de:tuda-tuprints-216719)。
   Eqs. (5.92),(5.94),(5.96),(5.97) 独立复述了 Bernoulli ODE、BCH 及
   Magnus(2) 截断；Sec. 5.6 明确用同一套 MR-IMSRG commutator
   equations 计算 Magnus/BCH 中的嵌套对易子。本地原文为
   `refs/theses/2022_Mongelli_IMNCSM_PhD.pdf`。
4. `method=magnus` 的有限步算法不是从上述论文反推出来的，而是
   现有 `imsrg++` 代码 `IMSRGSolver::Solve_magnus_euler()` 的直接语义：
   `Eta *= ds` 后调用 `BCH_Product(Eta,Omega)`。Vobig Sec. 4.6 和
   Appendix B.2 仅说固定步长 Euler 在形式上可以使用，但他们的实际
   实现选用自动步长 RKF45。

因此，将 `method=magnus_adaptive` 作为本项目的生产入口是“按上述文献
实现完整 ODE，再由 keV 稳定性验收决定”，不是论文规定的命令行
名称。本仓库用 Boost odeint 的 Dormand--Prince 5(4)，而 Vobig 用
GSL RKF45；两者都是带局部误差控制的显式嵌入式 5(4) Runge--Kutta，
但 Butcher tableau 不同，不宣称数值步进逐项复刻 Vobig。将原
`method=magnus` 保留为 Euler/BCH-product oracle 同样是本仓库根据
He4 `76.8 keV` 步长敏感性测试做出的工程决策，不是文献结论。

### 8.2 源码最小改动边界

1. `BCH_Transform`/`BCH_Product` 增加显式、非全局的 MR reference 上下文；
   无上下文时仍调用原 `Commutator::Commutator`，有上下文时调用
   `MRCommutator::Commutator`。
2. `Solve_magnus_euler`、完整 Bernoulli `MagnusDerivative`、
   `Solve_ode_magnus`、`NewOmega`、`GatherOmega`、中途 Hamiltonian 重建和
   最终 Hamiltonian `Transform` 全部传递同一上下文。遗漏任一调用点都会
   让第二阶及以上嵌套项错误退回 SR，因此用调用点测试锁定。
3. 第一版继续拒绝 goose-tank/factorized IMSRG(3) corrections、Brueckner
   BCH、一般附加 flowing operators 与 MR+IMSRG(3)；这些分支没有
   `lambda2` 公式验收，不能静默走 SR。
4. 不新增 MR 专用 ODE 参数，不改 White-NCSM、`Delta e!=0` 掩码、自然基、
   分段 `Omega` 或现有真空物化语义；生产自适应路径直接使用 imsrg++ 的
   `ds_0/dsmax/ode_tolerance/omega_norm_max`。

### 8.3 验收梯子

1. Python m-scheme oracle 保存固定 `Omega,O` 的每阶 BCH 以及固定
   `dOmega,Omega` 的每阶 BCH-product；随机 scalar J-scheme 展开后
   0B/1B/2B max-abs `<=1e-10 MeV`。
2. `Omega=0` 时首步严格为 `dOmega=ds*eta`；每步保持 `Omega/eta`
   anti-Hermitian、BCH 后 Hamiltonian Hermitian。
3. `He4/O16 Nrefmax=0` 的同一 MR driver 对每段 `Omega`、逐阶 BCH、
   `eta` 和最终 Hamiltonian 逐元素退化到原 SR `method=magnus` 与
   `method=magnus_adaptive`。
4. `Be8/C12 Nrefmax=0` 与 `He4/O16 Nrefmax=2` 逐项复现 Python oracle；
   短流与 direct flow 在共同小步长极限下收敛一致。IMSRG(2) 截断使两种
   长流路径不保证 bitwise 相同，不能把 direct flow 当成唯一真值。
5. 收紧步长/BCH 阈值十倍，最终 Hamiltonian 和下游 NCSM 低能谱变化
   `<1 keV`；J64/no2bpack 物化读回继续通过。最后运行完整 CTest、
   sanitizer 和 emax4/6 性能门。

### 8.4 首轮实现证据

commit `382eb3e2` 已把现有 scalar `BCH_Transform`/`BCH_Product`、分段
`Omega` 和最终 `Transform` 接到显式 MR dispatcher；SR 两参数 API 和默认
参数不变。独立 emax1、非零 `lambda2` 的随机 J-scheme 输入展开到 m-scheme
后，40 阶 MR-BCH 的 0B/1B/2B 最坏差为 `6.66e-16`，生产 BCH-product
最坏差为 `3.04e-18`。同一测试的 `lambda2=0` MR/SR transform 差为
`6.94e-18`，product 为逐位零差。

真实参考门在 `He4/Be8/C12/O16` 的
`Nrefmax=2/0/0/2` 上用 `ds=1e-4` 比较首步 `Omega` 与零阈值完整 MR-BCH：
四核 `Omega` 最坏差依次为
`1.694e-21/2.118e-21/1.271e-21/1.694e-21`，Hamiltonian 最坏差依次为
`1.776e-15/1.776e-15/1.776e-15/3.553e-15 MeV`。恢复原
`bch_transform_threshold=1e-9` 和 `bch_product_threshold=1e-4` 后，同一
生产 driver 的 J64 物化 0B/1B/2B 均逐位复现默认 C++ 路径。完整
`MRCorrelatedDriver` 用时 `409.10 s`。

`He4/O16 Nrefmax=0` 还在 `s=0,1e-4` 对 `flow_RK4` 和 `magnus` 分别运行
MR/SR executable；八组 J64 的 0B/1B/2B 差全部为零。这里验收的是
MR-Magnus 代数与生产入口，尚未替代后续步长/BCH 阈值稳定性、direct-flow
共同极限、完整流和 NCSM 谱门。

随后在同一个非零 `lambda2` 的 emax1 参考态上构造有明确单粒子能隙的
Hamiltonian，排除 White 分母 cutoff 的非光滑影响。生产 RK4 direct flow
与零阈值生产 Magnus 在共同区间
`ds=1e-2,5e-3,2.5e-3` 的 Hamiltonian 差依次为
`1.338e-4,3.368e-5,8.447e-6`，相邻比值为 `0.252,0.251`。这与首阶
Magnus 步相对 RK4 的 `O(ds^2)` 局部差一致，验证两条实现趋向同一个连续
MR 流；此前无隙随机 Hamiltonian 会触发 denominator cutoff，不能用作
ODE 阶数测试。

生产作业生成器现固定调用 `method=magnus_adaptive`，不增加 MR 专用 ODE 输入；
它启用原 driver 的 `write_omega=true`，使用作业私有 scratch，并逐个校验和
复制出的 `<intfile>_Omega_*` 文件。manifest schema v3 明确把这些文件标为
`segment` 变换及其累计流起止点，避免将重启段的 `Omega` 误认为从裸
Hamiltonian 开始的单一总变换。`UnitTestMRJobGenerator.py` 锁定脚本和
metadata，`UnitTestMRSRDriver.py` 已由实际 MR/SR Magnus executable 验证
`Omega` 文件存在且非空。

随后发现原一阶 Magnus Euler 在 `He4 Nrefmax=2, s=1` 上把 `domega` 从
`0.1` 收紧到 `0.01` 时，NCSM 三条低能级最多变化 `76.8 keV`，未通过数值
门，而 packed-format 误差仅 `0.90 eV`。这排除了输出格式问题。按照 Vobig
Eq. (4.6.4) 和 App. B.2，现已补齐完整 Bernoulli 导数与 RK45；随机相关
J-scheme 导数展开到 m-scheme 的最坏差为 `2.17e-19`。同一 He4 点用默认
`dsmax=0.5, ode_tolerance=1e-6` 和十倍收紧
`dsmax=0.05, ode_tolerance=1e-7`（BCH 阈值也收紧十倍）比较，Hamiltonian
最大矩阵元变化 `0.2708 keV`，三条 `Nmax=8` 能级最大变化 `0.4668 keV`；
J64/no2bpack 差最多 `1.09 eV`。因此生产生成器最终固定为
`method=magnus_adaptive`，Euler 仅保留为代数与收敛 oracle。

完整 Release CTest 的 `21/21` 项以 `411.32 s` 通过，其中四核相关参考
oracle 为 `411.32 s`；Bernoulli 停止条件改动后又单独通过
`MRMagnus/MRSRDriver/MRJobGenerator`，`He4/O16 Nrefmax=0` 的 adaptive
MR/SR driver 在 `s=0,1e-4` 的所有 J64 rank 均为零差。独立 ASan+UBSan
构建下这三项以 `10.03 s` 通过。真实相关 He4 adaptive driver 与同参数
library 路径的 0B/1B/2B 最坏差为 `3.55e-15 MeV`。

Vobig Eqs. (4.6.7),(4.6.8) 的 Bernoulli 相对项停止门将 emax4/6
`s=1e-4` 单线程流的 scalar commutator 数从固定 `k=8` 的 `65` 降至 `23`。
emax4 wall 从 `1.63` 降至 `0.64 s`，峰值约 `61.5 MiB`；emax6 从
`14.31` 降至 `5.59 s`，峰值约 `390.1 MiB`。两空间停止版与固定做到
`k=8` 的 J64 输出在 0B/1B/2B 上均逐位相同。这是计算能力/实现门，不把
`s=1e-4` 称作已脱耦的物理 Hamiltonian。

### 8.5 Magnus 级数失败与拒步

本阶段的生产收尾改动将 Vobig Eq. (4.6.8) 从只用于提前
停止的条件升级为真正的失败信号。对非零生产阈值，嵌套对易子
范数不再严格递减，或到 `k=8` 仍未达到 Eq. (4.6.7) 的相对项
门，都抛出带失败阶数和相邻范数的 `MagnusSeriesError`；不再将最后
一个部分和静默当作导数。`relative_threshold=0` 仍固定计算到
`k=8` 且不抛收敛异常，保留独立代数 oracle。

`Solve_ode_magnus` 在任一 Dormand--Prince stage 遇到该异常时，显式把
`state/current_s` 恢复到 `Omega.back()/s`，清除 FSAL 导数缓存，把
trial step 减半后重试。若失败就在已接受起点发生，或同一
Omega 段累计八次真实缩步仍不能安全前进，则先用 MR-BCH
物化该已接受段，再以局部 `Omega=0` 继续。普通
`omega_norm_max` 分段后也现在清除 FSAL 缓存，避免把上一个局部
坐标的末导数误用为新段首导数。

固定随机非零 `lambda2` 模型已分别强到：

- 非递减级数拒绝，以及不可能的 `1e-30` 项阈值在 `k=8`
  仍未收敛时拒绝；
- 起点 `Omega` 失败时自动分段，其最终 Hamiltonian 与先显式
  物化同一段再流的 oracle 差小于 `1e-13 MeV`；
- 近边界段在 RK45 中间 stage 于同一段累计触发八次
  拒步，然后有界
  分段恢复并到达目标 `s`，没有耗尽浮点步长。

真实 `He4 Nrefmax=2, s=1` 默认/十倍收紧流均无需级数拒步；
二者 Hamiltonian 最大矩阵元差 `0.2543 keV`，`Nmax=8` 三条
低能级最大差 `0.4807 keV`，driver/library 最大差
`5.33e-15 MeV`，J64/no2bpack 最大差 `0.907 eV`。上述强制恢复、
SR driver 与作业生成回归还在独立 Debug
`-fsanitize=address,undefined` 构建下全部通过；完整 Release
CTest `21/21` 以 `499.20 s` 通过，其中相关参考态验收项
`MRCorrelatedDriver` 为 `405.61 s`。

生产作业仍不要求 MR 专用 ODE 输入：默认 profile 完全继承
`IMSRGSolver` 的 `ds_0/dsmax/ode_tolerance`。为执行独立的十倍收紧
验收，单点生成器另提供 `--mr-tight-validation`布尔 profile；它只在
复算脚本中写入既有 `ode_tolerance=1e-7`，不改步长参数，并将
唯一 override 和 profile 写入 metadata。这不把 ODE 参数重新变成
默认生产必填项。tight 复算必须以 default 流实际停止的 `s`
为固定 `target_s`，并设 `eta_criterion<=1e-20` 禁用残差提前停止；
否则两个容差会在略有不同的 `s` 越过门槛，其 Hamiltonian 差混入
真实流变化，不再是纯 ODE 误差。生成器对该条件作强制检查。

### 8.6 emax2 四核严格 Magnus 门禁

point7 唯一生产 checkout 已回到
`/tns/mengziyan/mr-imsrg`；原有被 Git 忽略的结果保留，旧的 J-scheme
direct-flow 结果通过 `result/mr-jscheme-flow-legacy-e604e181` 链接统一
访问。生产 executable SHA-256 为
`a70a28582d232d675c41a8db26a327b6d29bde8e5cf0317695d7d26745218de5`。
配套 `libIMSRG.so` SHA-256 为
`c7eb668df585e96c29b0c7e659288ef16cbeecefcf358d9d390399488261927b`；
最终聚合器要求四核全部逐字匹配这两个冻结哈希，而不只要求彼此相同。
该 executable 构建于 `b15ac899`；从该提交到当前验收工具提交，
`BCH.cc`、`IMSRGSolver.cc`、`MRCommutator.cc`、`MRReference.cc`、
`Generator.cc`、`imsrg++.cc`、`Operator.cc` 和 `TwoBodyME.cc` 的
`git diff --quiet` 返回零。`src/` 后续变化只有 CMake 注册和独立
`imsrg_operator_validate`，没有生产物理代码漂移；因此 default/tight
继续使用同一已冻结 executable/library，不为版本字符串重链后重跑长流。
He4 的 bare J64/JREF 分别复现已有严格门禁 SHA-256
`55a3c161...44e2`/`1c9934aa...0898`；四核均从同一固定 NNLOopt
minipack 按各自 `A` 生成 bare J64。

`Be8 Nrefmax=0` 和 `C12 Nrefmax=0` 已首先完成。default jobs
`100602/100604` 在 `s=330.08347/226.53612` 停止，最终
`||eta||/||eta(0)||=9.78935e-7/9.71606e-7`。同终点 tight jobs
`100618/100616` 禁用提前停止并将 ODE 容差从 `1e-6` 收紧为
`1e-7`。NCSM `Nmax=0` 三态最大变化分别为
`0.648406/0.769901 keV`，通过 `<1 keV` 门禁。完整 J64
default/tight 矩阵元最大差为：

- Be8：0B `0.58385 keV`、1B `0.55189 keV`、2B `0.28986 keV`；
- C12：0B `1.68784 keV`、1B `0.83744 keV`、2B `0.27109 keV`。

C12 的零体矩阵元差大于 `1 keV`，但后 NCSM 谱差为
`0.770 keV`；这两个量已分开记录，不用零体项代替谱验收。
J64/no2bpack 三态最大差为 `1.145/3.305 eV`，符合
float32 packed 路径的 `5 eV` 四核读回窗；高精度物理比较始终用
lossless J64。两核所有级数拒步计数为零。

新增 `prototype/mrimsrg/summarize_magnus_gate.py` 统一读取 flow、
metadata、resource usage、拒步日志与 `Omega` 段，比较完整
J64 0B/1B/2B，并分别对角化 J64/no2bpack 后输出 JSON。
`MRMagnusGate` 回归固定其 flow/resource/NCSM 解析，并把 default/tight
均存在至少一个非空 `Omega` 段列为硬门禁。新增
metadata 门同时要求两条流均为 `magnus_adaptive`、default 无 ODE override、
tight 只含 `ode_tolerance=1e-7`、`eta_criterion<=1e-20`，且 tight 的配置
终点逐字对应 default 的实际停止点；因此不能用两个偶然在同一时刻输出、
但实际 profile 不同的目录冒充十倍收紧验收。单核 schema 因此升级为
`mrimsrg_magnus_gate_v2`；最终聚合器要求完整 13 项门且全部为真，拒绝
缺少新门的旧 JSON，即使旧文件的顶层 `passed` 曾为真。新增
`imsrg_operator_validate` 通过同一生产 `Operator::ReadBinary` 读回每段，
要求 scalar rank/parity、anti-Hermitian metadata、有限矩阵元、零 0B，
以及完整 1B/各 J-channel 2B 的数值 anti-Hermiticity 均在 `1e-10`
以内；人为写入 `1e-3` 的 1B 对角元会被负测试拒绝。He4/O16
default jobs `100600/100606` 正在运行，完成后执行同一门禁。Be8/C12
已在 point7 计算节点用远程原始产物生成
`result/mr-jscheme-magnus/gates/Be8.json` 与 `C12.json`，两份
`passed=true`；验证器由独立的 `prototype/mrimsrg/build-current`
构建，不重链生产 executable。加入 `MRMagnusGate` 后的本机 Release
完整 CTest 为 `22/22`。再加入落盘 Omega 数值 checker 后重跑仍为
`22/22`，总用时 `491.61 s`。point7 原始产物中 Be8 default/tight
分别有 `20/21` 段、C12 有 `12/13` 段；所有段的 1B/2B 数值
anti-Hermiticity 违例均为零，0B 最大绝对值仅
`5.063e-22`，完整 gate 保持 `passed=true`。checker 另在独立
ASan+UBSan Debug 构建中读取 Be8/C12 的真实 production `Omega_0`；
两次均无 sanitizer 诊断，1B/2B 违例为零，0B 分别为
`1.212e-23/4.019e-23`。

He4/O16 长流的 scalar commutator 计数超过六位后，原生
`WriteFlowStatus` 的七字符 `Ncomm` 紧接前一个 16 字符 `Eta_3` 字段；
因此纯 `split()` 会把 `0.0000000001184939` 误当成一个浮点数并把计数读成
零。gate reader 现按源码的 `5,12,9x16` 固定列宽解析到 `Eta_3`，随后
读取不限于 `setw(7)` 的完整十进制计数；真实 He4 七位行和合成八位溢出行
均进入回归，既不改动流也不丢失计数。
单核 JSON 同时记录该 reader 文件本身的 SHA-256；最终聚合要求四核一致，
避免同一 v2 schema 下混入修复前后的两种解析实现。

最终四核完成判据由 `prototype/mrimsrg/aggregate_magnus_gates.py`
聚合，不从日志手工摘录。它要求恰有 `He4/Be8/C12/O16`，对应
`Nrefmax=2/0/0/2`、NCSM `Nmax=8/0/0/2` 和各三态，并强制固定
NNLOopt SHA-256、生产 executable、`libIMSRG.so`、下游 NCSM/Omega
validator 与正式阈值逐字一致；四核 JREF 及各自 A-dependent bare J64
也按已冻结 SHA-256 逐核检查。任一单核
`passed=false` 或 provenance 不一致即拒绝。输出同时汇总四核最大
残差比、tight/default 谱差、packing 差、Omega anti-Hermiticity、
拒步计数、作业号和产物哈希。
