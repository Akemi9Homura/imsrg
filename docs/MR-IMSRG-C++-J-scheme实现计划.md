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
- 首先复用现有直接 flow 和 ODE driver；Magnus 是后续独立门禁。

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

## 6. 实现顺序与 commit 边界

1. 文献/C++ 调研表与 J-coupled density 约定；
2. MR reference 读入及 J/m 密度闭环；
3. 单个命名 `lambda2` contraction 及其随机 oracle；
4. 完整 MR commutator 与 SR 退化；
5. White-NCSM generator 与 `Delta e` 掩码；
6. solver/driver 集成和真实体系 `s=0`/short-flow 验收；
7. 完整流、物化、NCSM 与性能验收。

每个边界在自身测试通过后立即独立 commit。任何一个未通过 Gate 1--4
的中间实现都不能用于新物理结果。

当前进度（2026-08-25）：步骤 1 完成；步骤 2 的显式 `MRReference`、
J-coupled 恒等式检查、正规序往返及 `mrimsrg_jref_v1` 文件 reader 已完成；步骤 3--4
的慢速参考、块矩阵生产收缩、完整 0B/1B/2B MR commutator 及精确 SR
退化门禁已完成。步骤 5 的 White-NCSM 方向分子、MR EN 分母和
`Delta e != 0` 掩码也已完成并对 Python m-scheme 验证。当前不得跳到
完整流出数；步骤 6 的 solver 显式 MR context 和直接 RK4/ODE RHS 接线
已完成，剩余边界是 `imsrg++` driver 参数、自然基 Hamiltonian 变换及真实
体系 `s=0`/short-flow 验收。

## 7. 错误定位原则

- SR 极限失配：先检查是否真正复用了现有 contraction、occupation 和
  TBME normalization，不通过加核特例修正。
- MR 失配而 SR 通过：按命名 `lambda2` 收缩逐项查相位、6j、pair
  normalization 和 J-channel 选择。
- 只有流后失配：用共同 Euler/RK4 隔离 RHS 与 ODE driver，再检查对称化和
  停止条件。
- 能量接近但矩阵元不同：仍判为失败，禁止以截断误差或基底差异
  做无证据解释。
