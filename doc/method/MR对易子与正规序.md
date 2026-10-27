# MR 对易子与正规序

代码：`src/MRReference.{hh,cc}`（正规序、λ2 收缩）、`src/MRCommutator.{hh,cc}`（λ2 项）、`src/BCH.cc`（MR 分派）。
公式来源与逐项系数见 [公式与代码映射](公式与代码映射.md) 第 3 节；这里只说结构、为什么这样写、怎么调用。

## 1. 设计：SR 对易子 + λ2 增量

在自然轨道基下，γ1 对角（占据 n_a），MR-IMSRG(2) 的 [X, Y] 等于“分数占据的 SR IMSRG(2) 对易子”再加上含 λ2 的项。
所以 `MRCommutator::Commutator(X, Y, ref)` 的结构是：

```
result = Commutator::Commutator(X, Y)        // imsrg++ 原有 SR 对易子，占据取 n_a
if ‖λ2‖ == 0: return result                   // SR 退化：逐位等于 SR
标量：result.OneBody  += IV + V + VI          // comm221_lambda2
      result.ZeroBody += ¼ Σ (2J+1) tr(C2^J λ2^J)，C2 为完整的两体结果
      result.ZeroBody += comm220_lambda3        // 参考态带 λ3 时
张量：result.OneBody  += comm221_lambda2_tensor
```

理由：

- SR 部分的大量对易子（含 imsrg++ 的全部优化、并行、张量项）直接复用，不复制公式；MR 只维护 λ2 增量；
- λ2 = 0 且无 λ3 时直接返回 SR 结果，不做“加零”运算，SR 退化是逐位相同而不是近似相同（`MRSRDriver` 依赖这一点）；
- 两体输出不含 λ2（λ3 = 0 时 Hergert 式 (51) 无 λ2 项），所以 MR 只改一体、零体。

λ3：参考态文件带 λ3 时加零体项 ¼ Σ (X^{ab}_{cd} Y^{kl}_{am} − Y^{ab}_{cd} X^{kl}_{am}) λ^{bkl}_{cdm}（Hergert 式 49），
加上后标量 [X, Y] 的零体等于 ⟨Ψ|[X,Y]|Ψ⟩（λ2 二次项与 λ4 项相消）；一体、两体不含 λ3。J 耦合形式与实现见 三体力与λ3.md。
生成元分子不含 λ（见 生成元.md）。

## 2. 标量的 λ2 项

一体（Gebrerufael 式 (4.89b) 的 IV、V、VI，对应 Hergert 式 (50) 的四个 λ2 指标式）：

| 项 | 拓扑 | 实现 |
|---|---|---|
| IV | 直接通道 X λ Y，旁观轨道收缩 | 每个二体通道的 X^J λ^J Y^J 块乘积，再按旁观轨道偏迹 |
| V | 交叉耦合 X̄ λ̄ Ȳ | X、λ、Y 分别 Pandya 变换，在交叉耦合通道做三矩阵乘积 |
| VI | λY、λX 后与 X、Y 收缩 | 普通通道先形成 λY、λX，再直接偏迹 |

零体：在完整的两体对易子（1B–2B 与两种 2B–2B 拓扑都已累加）上做 `ContractLambda2`。只收缩 2B–2B 子项是错的（Hergert 式 (49)）。

实现分两版：`comm221_lambda2_reference`（慢速逐轨道循环，正确性基准）和 `comm221_lambda2`（生产，块矩阵乘法、OpenMP）。
两者对 m-scheme 随机镜像的误差为 1e-12 量级。归一化对存储下的系数（IV 0.5、V 0.5、VI −1）与论文写法不同，原因见映射文档。

文献排印：Mongelli 式 5.183 印作 (−2) 且无 Ĵ²，应为 +2Ĵ²；式 5.160 也漏了 Ĵ²。代码按 Gebrerufael 写并已验证，不要按 Mongelli 改。

## 3. 张量的 λ2 项

用途：Ω（标量）与秩 L > 0、不改变同位旋的张量 O（E2、M1、E1 等）的 [Ω, O]，算符演化需要。λ2 只进一体项（零体项只对标量存在）。

没有照抄论文的 J 耦合张量公式（Vobig 第 5 章），而是从 m-scheme 式（Hergert 式 (50)）自行重排，每一块对暴力 m-scheme 求和逐项核过：

Z = PT(XλY − YλX) + C(X, PT(Yλ − λY)) − D(Y, PT(Xλ − λX)) + R(X, Y) − R(Y, X)

| 记号 | 含义 | 来源 |
|---|---|---|
| PT_ij | Σ_c W_{ci,cj}，系数 Ĵ1Ĵ2 {J1 J2 L; j_j j_i j_c} (−1)^{j_j+j_c+J1+L} | imsrg++ SR 221 张量项的部分迹 |
| C | −Σ_ab N_ab Σ_J (2J+1) {j_i j_j L; j_a j_b J} X^J_{(ib),(aj)} | comm121st 的 2B×1B 收缩（N 为张量一体量） |
| D | Σ_ab S_ab Σ (−1)^{j_j+j_a+J1+L} Ĵ1Ĵ2 {J1 J2 L; j_j j_i j_a} T^{J1J2}_{(bi),(aj)} | 同上（S 为标量一体量） |
| R | 环形项：R = −Σ_t Σ_{JJ'} ĴĴ' {J J' L; j_j j_i j_t} (−1)^{j_i+j_t+J'+L} [x̄ λ̄ ȳ]^{JJ'}_{(it),(jt)} | 交叉耦合 |

交叉耦合：标量 X̄^J = Σ_{J3} (2J3+1) {j_a j_b J; j_c j_d J3} X^{J3}_{(ad),(cb)}；
张量 X̄^{JJ'} = Σ_{J3J4} ĴĴ'Ĵ3Ĵ4 (−1)^{j_b+j_d+J'+J4} {j_a j_d J3; j_b j_c J4; J J' L} X^{J3J4}_{(ad),(cb)}。

约定：imsrg++ 的 Edmonds Wigner–Eckart、`GetTBME_J` 的非归一化矩阵元、标量按普通矩阵元存。

实现要点（为什么快，详见 性能与并行.md）：

- 矩阵元经 `PairTable`/`BlockAccess` 直接按块指针读，不走 `GetTBME_J` 的 map；6j/9j 按 j 组合制表、跨调用保留；
- λ2 只在“活跃轨道”（出现在某个非零 λ2 元里的轨道，即 Nmax_ref 盒子内部分占据的轨道）上非零，含 λ 的乘积只跑活跃 ket；
- O 厄米、Ω 反厄米时，PT(OλΩ)、R(O, Ω) 由 PT(ΩλO)、R(Ω, O) 的厄米共轭得到：P_ij + h_Z (−1)^{j_i−j_j} P_ji；
- 交叉耦合矩阵按 (J, 宇称, Tz 差) 分块；P、W2 等中间量在偏迹任务内逐块生成、用完即丢（`std::function` 块生成器）。

## 4. MR 正规序

两体算符的 MR 正规序只需 γ1、λ2；3N 的 MR-NO2B（含 λ3 零体项）在 `MRThreeBody` 里，见 三体力与λ3.md。


`MRReference::NormalOrder(vacuum_op)`：输入必须已在自然轨道基（驱动先做 `TransformOneAndTwoBody(NaturalOrbitTransformation)`）。

- 先做分数占据的 SR 正规序 `DoNormalOrdering2(+1, all_orbits)`：E += Σ n_a h_aa + ½ Σ n_a n_b V_abab，f_pq = h_pq + Σ_a n_a V_paqa，Γ = V；
- 标量再加 E += ¼ Σ λ2 V（`ContractLambda2`）；
- 张量（L > 0）不加 λ2 项：零体部分只对标量存在，一体部分在自然基下与 SR 相同。

`UndoNormalOrder` 是严格逆（`DoNormalOrdering2(−1)` 并减去 λ2 收缩）。导出时：撤销 → 用自然轨道变换的转置转回 HO → 写文件。
驱动打印的 `IM-FCIQMC export check` 就是“撤销再正规序”的往返差。

为什么在自然基里做：γ1 对角后 SR 收缩只需占据数 n_a，可以原样复用 imsrg++ 的正规序与对易子；
否则所有收缩都要带非对角 γ1，相当于重写全部对易子。

## 5. BCH/Magnus 中的 MR 处理与防护

`BCH::BCH_Transform`、`BCH_Product`、`MagnusDerivative` 都多了一个 `const MRReference*` 参数；非空时嵌套对易子走 `MRCommutator::Commutator`。
`IMSRGSolver::SetMRReference(ref)` 设好后，求解器的所有对易子、BCH 都自动走 MR。

防护（MR 下报错而不是静默算错）：

| 情形 | 行为 | 原因 |
|---|---|---|
| 三体算符或全局 `use_imsrg3` | `MRCommutator::Commutator` 抛异常 | MR 增量只有 MR-IMSRG(2) 版本，SR 会返回缺 MR 项的三体部分（3N 输入在流前就折成两体） |
| BCH 在 40 层嵌套内不收敛 | 抛异常（阈值为 0 的固定阶模式除外） | SR 只警告；MR 下发散的变换不应写出 |
| goose-tank、factorized、Brueckner BCH | 抛异常 | 没有 MR 版本 |
| `MagnusDerivative` 级数不收敛 | 抛 `MagnusSeriesError`，ODE 拒步 | 见 驱动参数与输出.md 1.2 |
| 求解器 `FlowingOps` 多于一个 | 抛异常 | MR 下算符只经 Ω 变换，不随流积分 |

驱动层的启动检查已经拒绝这些配置，库函数层的检查只拦截 Python 脚本或其他代码绕过驱动的情形。

## 6. Python 接口

pyIMSRG 暴露了 `MRReference`（`ReadOccupationMap`、`ReadBinary`、`WriteBinary`、`EmbedInModelSpace`、`NormalOrder`、`UndoNormalOrder`、
`Validate`、`Lambda2`、`has_lambda3`、`Lambda3Blocks`、`reference_energy`）和函数 `MR_Commutator`、`MR_comm221_lambda2`、
`MR_comm221_lambda2_reference`（返回 IV、V、VI 分量）、`MR_comm220_lambda3`、`MR_NormalOrderThreeBody`（返回 (算符, 报告)），测试与诊断脚本都用它：

```python
ms = pyIMSRG.ModelSpace(2, "He4", "He4")
ms.SetReferenceOcc(pyIMSRG.MRReference.ReadOccupationMap(ms, "He4.jref"))
ref = pyIMSRG.MRReference.ReadBinary(ms, "He4.jref", 1e-10)
Z = pyIMSRG.MR_Commutator(X, Y, ref)     # 标量或张量 Y 都走这个入口
```
