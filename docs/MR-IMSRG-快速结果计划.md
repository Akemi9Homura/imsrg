# NNLOopt emax2 MR-IMSRG 快速结果计划

## 1. 交付目标

本计划只追求最短的可信闭环，不追求高效、通用或生产级 MR-IMSRG：

```text
Nmax 截断 NCSM 参考态
  -> gamma1/gamma2/lambda2
  -> m-scheme MR-IMSRG(2)
  -> 真空正规序 0B/1B/2B Hamiltonian
  -> NCSM 确定性读回与对角化
  -> 可选：同一文件交给 FCIQMC
```

完成定义是为 `He4`、`Be8`、`C12`、`O16` 各产生一份 A-dependent 演化 Hamiltonian，并由 NCSM 成功读回；至少 `Be8`、`C12` 使用非零 `lambda2` 的多参考流。FCIQMC 不是第一批结果的阻塞条件。

## 2. 冻结的物理设置

| 项目 | 固定值 |
|---|---|
| 相互作用 | NNLOopt；本机命名 `N2LO_opt` |
| 核力 | `/home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack` |
| SHA-256 | `76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84` |
| HO 频率 | `hw=20 MeV` |
| 单粒子截断 | `emax=2` |
| 二体截断 | `e2max=4` |
| 相互作用 rank | NN-only；无显式 3N |
| 单粒子基 | HO |
| 质心项 | `BetaCM=0` |
| 其他约定 | 无库伦、无额外核子质量修正 |
| MR 密度截断 | 保留 `lambda2`，设 `lambda3=0` |
| 流动算符截断 | 参考态正规序 0B/1B/2B |
| 流形式 | 直接积分，不先做 Magnus |

`minipack` reader 根据目标 `A` 加入内禀动能，因此必须分别为四个核读取并演化，不能把一个核的物化 Hamiltonian 用于另一个核。该普通 NN `minipack` 也不能当作 `fmt2=no2bpack` 读取。

## 3. 两批参考态

### 第一批：最快四核结果

四个核统一使用 `Nrefmax=0`，降低参考态生成和一般非对角一体密度带来的复杂度。

| 核 | 第一批角色 | 说明 |
|---|---|---|
| `He4` | SR 极限与 I/O 基准 | 闭壳 `Nrefmax=0` 是单 Slater；要求复现已知 emax2 FCI 基态。 |
| `Be8` | 非平凡 MR 主测试 | `J=0+` 的 `Nrefmax=0` 子空间是多行列式；必须得到非零 `lambda2`。 |
| `C12` | 非平凡 MR 文献测试 | 与原始 IM-NCSM 的标准开壳案例一致；必须得到非零 `lambda2`。 |
| `O16` | 较重 SR 极限对照 | 闭壳 `Nrefmax=0`，用于观察同一实现的质量数扩展。 |

这里必须在输出中把 `He4/O16` 标记为 SR-control，不能把它们的第一批结果描述成非平凡 MR。

### 第二批：闭壳核相关参考态

第一批闭环通过后，只增加：

- `He4, Nrefmax=2`；
- `O16, Nrefmax=2`。

这两点用于验证闭壳体系采用相关 NCSM 参考态后的真正 MR 路线。第二批不得反过来阻塞第一批四核物化结果。

## 4. 最小软件结构

优先新增独立原型目录，例如：

```text
prototype/mrimsrg/
  README.md
  interaction_io.py
  reference_io.py
  densities.py
  normal_order.py
  commutator.py
  generator.py
  flow.py
  export.py
  run.py
  tests/
```

语言优先 Python/NumPy/SciPy；若现有 C++ reader 的复用能明显缩短工作量，可以提供一个只负责格式转换的小工具。第一版不要为了统一架构而把 MR 数据结构塞进现有 J-coupled `Operator`/`Commutator` 主路径。

内部张量可以低效，但必须利用粒子种类、`M`、宇称等明显的零块，避免把 `emax=2` 的所有六指标循环写成无条件稠密扫描。性能目标只是四个指定体系能完成，不要求外推到更大 `emax`。

## 5. 里程碑与退出条件

### M0：冻结普通 Hamiltonian 和基准

任务：

1. 读取固定 `minipack` 的 header，验证 `hw=20, emax=2` 和校验和。
2. 固定轨道表、质子/中子和 m-scheme bit/index 顺序。
3. 对每个目标质量数按相同约定构造内禀 Hamiltonian。
4. 用现有确定性 NCSM/FCI 程序记录未演化谱和基维数。

退出条件：

- `He4` 全 emax2 FCI 基态复现 `-20.33883250 MeV`；
- 同一普通 Hamiltonian 经读出/写回后谱不变，误差不超过 `1e-8 MeV`；
- 四个核的基准命令、轨道顺序、输出和 SHA-256 写入一个小型 manifest。

### M1：NCSM 参考态和密度接口

任务：

1. 读取最低 `J=0+` NCSM 参考态或直接读取其 `gamma1/gamma2`。
2. 构造

   \[
   \lambda^{(1)}=\gamma^{(1)},\qquad
   \lambda^{(2)}=\gamma^{(2)}-\mathcal A\{\gamma^{(1)}\gamma^{(1)}\}.
   \]

3. 把所有轨道、相位、RDM 指标约定写进文件 header/metadata。

退出条件：

- `Tr(gamma1)=A`；
- Hermiticity、反对称性和一/二体收缩关系在 `1e-10` 相对精度内成立；
- `He4/O16 Nrefmax=0` 给出 `lambda2=0`；
- `Be8/C12 Nrefmax=0` 给出稳定、明确非零的 `||lambda2||`；
- 直接由参考波函数计算的能量等于 RDM 收缩能量。

### M2：MR 正规序往返

任务：

1. 将普通

   \[
   H=E_0+t+V
   \]

   转成参考态正规序的 `E,f,Gamma`。
2. 从 `E,f,Gamma` 转回普通真空正规序 `E0,t,V`。
3. 保留零体常数，禁止下游静默丢弃。

退出条件：

- `E=<Psi_ref|H|Psi_ref>`；
- 两次转换后所有 0B/1B/2B 元素相对误差不超过 `1e-10`；
- 转回文件由 NCSM 读入后复现 M0 谱。

### M3：MR-IMSRG(2) commutator

任务：

1. 实现保留线性 `lambda2`、设置 `lambda3=0` 的 0B/1B/2B 对易子。
2. 用 `refs/qcombo/examples/MR_IMSRG2.ipynb` 保存的符号输出核对指标、系数和交换项。
3. 用小随机反对称张量和显式 Fock-space 矩阵对易子做数值比较。

退出条件：

- 单 Slater 极限逐项退化为仓库现有 SR-IMSRG(2)；
- 随机小模型所有可比较张量元素误差不超过 `1e-10`；
- 对易子保持预期 Hermitian/anti-Hermitian 类型和二体反对称性。

### M4：单一生成元和直接流

任务：

1. 只实现 Vobig 论文 6.5.28--6.5.34 及其命名表中实际使用的
   `White-NCSM`：生成元分子保留 `D1/D2` 的领头占据数项，不含
   `lambda^(2,3,...)`；Epstein--Nesbet 分母同样舍去文献标为
   `O(lambda2)` 的修正。这是生成元截断，不是对易子截断；
   MR-IMSRG(2) 对易子仍保留线性 `lambda2`。
2. 应用 IM-NCSM 掩码：

   \[
   e(p)\ne e(q),
   \]

   \[
   e(p_1)+e(p_2)\ne e(q_1)+e(q_2).
   \]

   对 `gamma1` 非对角的相关参考态，按 Vobig Sec. 6.5.3 先转到
   球形自然轨道基，在该基中用自然轨道输出槽位继承的
   `e=2n+l` 标签施加掩码。不得把已经用自然基 EN 分母加权的
   生成元先变回 HO 基再投影；两个操作不对易，会改变已发表的
   `White-NCSM` 生成元并产生伪残差长尾。流结束后再把
   Hamiltonian 变回原 HO 轨道顺序供 NCSM 读回。

3. 直接积分 `dH/ds=[eta,H]_(0,1,2B)`，记录每步残差和对称性误差。
   和现有 `IMSRGSolver` 的 `Eta.Norm()` 停止条件一致，正式停止量
   `Rgen` 是 Vobig 式 (6.5.28)--(6.5.29) 中完整的、带 EN 分母并反对称化后的
   掩码生成元范数。另外分别记录不带分母、不含不可约密度的
   `Rnum` 和含线性 `lambda2` 的严格 `Rstrict=D-D^dagger`，
   用于量化生成元截断差异。对分数占据，分母加权与反对称化不对易，
   所以 `Rnum` 不是这个生成元的正式固定点条件。三种量不得混称。

退出条件：

- `H` 始终 Hermitian，`eta` 始终 anti-Hermitian；
- 目标掩码内的 White-NCSM 生成元范数 `Rgen` 相对初值至少下降
  `1e-6`，并同时报告 `Rnum` 和严格 `Rstrict` 的变化；
- 同 HO 量子数的参考空间内部耦合没有被错误清零；
- ODE 容差缩小十倍后，后对角化能量变化小于 `1 keV`；
- 找到稳定流参数窗口；不要求盲目积分到无限大。

### M5：物化输出和四核结果

任务：

1. 将流后的 `E,f,Gamma` 转回普通 0B/1B/2B。
2. 输出下游可读格式和完整 metadata。
3. 先用 NCSM 读回，再按 `He4 -> Be8 -> C12 -> O16` 运行。

每个核至少保存三个点：

- `s=0`；
- 稳定平台中的中间点；
- 选定最终点。

退出条件：

- `s=0` 文件严格复现 M0；
- 后 NCSM 能读取所有四份演化 Hamiltonian；
- 报告基态和最低若干 `J=0+`/可用低激发态，不把 `E(s)` 当最终能量；
- 报告原始与演化 Hamiltonian 的 NCSM 收敛、解耦残差和本征值漂移；
- 同一文件可选地由 FCIQMC 读入，但不因 FCIQMC 调参延迟本里程碑。

### M6：闭壳相关参考态补充

只对 `He4/O16` 改用 `Nrefmax=2`，重复 M1–M5。若 `gamma1` 在 HO 基非对角，支持一般一体密度，并只在流方程内部按 M4 做临时自然基协变求值；不要为了这一补充先开发自然轨道 NCSM 变分优化或独立生产链。

退出条件：

- 两个核均有非零 `lambda2`；
- 与各自 `Nrefmax=0` 控制结果并列报告；
- 量化参考态改变、IMSRG(2) 截断和后 NCSM 结果的变化。

## 6. 四核结果表模板

每个计算最终至少填充：

| 核 | `Nrefmax` | `||lambda2||` | `s` | `Rgen/Rgen0` | `Rnum/Rnum0` | `Rstrict/Rstrict0` | `E0B(s)` | 后 NCSM 基态 | 原始有限空间基态 | 漂移 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| He4 | 0 / 2 | | | | | | | | `-20.33883250`（全 emax2） | |
| Be8 | 0 | | | | | | | | | |
| C12 | 0 | | | | | | | | | |
| O16 | 0 / 2 | | | | | | | | | |

另附：Hamiltonian 文件名、核力 SHA-256、参考态/RDM SHA-256、ODE 容差、生成元和 `lambda3=0` 标记。

## 7. 明确不做的工作

第一批结果出来前不做：

- 高效 J-scheme MR-IMSRG；
- 完整 `lambda3` 或 3-RDM 接口；
- 显式 3N 和 MR-NO2B；
- Magnus 与观测算符一致演化；
- natural-orbital NCSM 变分优化或独立生产工作流（不包括 M4 必需的内部临时自然基求值）；
- 奇核或 `Jref != 0`；
- MPI/OpenMP 优化和大型参数扫描；
- 修改生产 `gen_job.py` 或提交集群大作业；
- 为 FCIQMC walker/initiator 参数做优化。

任何新增任务若不能直接缩短 M0–M5 的完成路径，应移到第一批四核结果之后。

## 8. 风险与快速处理

| 风险 | 快速处理 |
|---|---|
| `minipack` 与本仓库 reader 不兼容 | 复用已验证 reader 做一次性转换，不把普通 `minipack` 冒充 `no2bpack`。 |
| NCSM 波函数 reader 相位不一致 | 用能量期望值、RDM 收缩和小空间 overlap 三重核对。 |
| 一般 `gamma1` 增加公式复杂度 | 第一批统一 `Nrefmax=0`；闭壳 `Nrefmax=2` 留到 M6。 |
| MR 对易子指标错误 | QCombo 输出 + 显式 Fock-space 小模型双重 oracle。 |
| 流在大 `s` 漂移 | 在后 NCSM 能量稳定平台选终点，并报告平台宽度；不强求 `s->infinity`。 |
| 下游格式丢零体项 | 在 header 中强制保存 `E0`，读回测试单独核对总能量。 |
| O16 低效原型过慢 | 保持 `emax=2`，利用量子数零块；只做单个生成元和少量流参数 checkpoint。 |
