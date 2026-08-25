# MR-IMSRG 当前 TODO

## P0：生产 MR 路径严格退化到 `imsrg++` SR-IMSRG(2)

这是当前最高优先级门禁。比较必须运行现有生产 `prototype/mrimsrg/`
路径；不得另写一套只为通过测试的 SR 生成元、对易子或 flow。只允许增加
不改变生产物理定义的导出、诊断和逐项比较接口。

固定输入：

- `N2LO_opt`, `hw=20 MeV`, `emax=2`, `e2max=4`, NN-only；
- 固定 minipack SHA-256
  `76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84`；
- `BetaCM=0`，无 Coulomb、无额外质量修正；
- HO 基，A-dependent 内禀动能；
- `He4/O16, Nrefmax=0` 闭壳单 Slater 参考；
- `lambda2=lambda3=0`；
- `imsrg++` 的 `white` 生成元和 `Epstein_Nesbet` 分母。

### A. 冻结完全相同的输入

- [ ] 记录本仓库和 oracle checkout 的 commit、编译器及线性代数库版本。
- [x] 从同一个 double-precision J-coupled Hamiltonian 构造两边的初始算符，
      禁止用 float32 `no2bpack` 做代数级比较。
- [x] 核对 J-orbit/m-orbit 表、质子中子约定、占据数、core/particle 划分。
- [x] 核对 `He4/O16` 的 `lambda2` 最大绝对值为数值零。
- [x] 证明闭壳参考的所有非零 `ph/pphh` 分子均满足 `Delta e != 0`，因此生产
      IM-NCSM 掩码与 `imsrg++` SR 非对角块在这两个核上完全相同。

### B. `s=0` 逐项退化

- [x] 比较正规序后的 `E`、`f`、`Gamma`，并分别报告 max-abs、Frobenius
      和最坏矩阵元索引。
- [x] 比较每个被选择的 1B/2B Epstein--Nesbet 分母。
- [x] 检查两边是否触发 `1e-6 MeV` 小分母 cutoff；若触发，生产实现必须与
      当前 `src/Generator.cc` 的实际符号约定一致。
- [x] 比较 `eta(s=0)` 的 1B/2B 全部矩阵元。
- [x] 比较 `dH/ds=[eta,H]` 的 0B/1B/2B 全部矩阵元。
- [x] 对 commutator 的 0B/1B/2B 各收缩项增加可单独开关的诊断，定位首个
      不一致项，不能只比较总和。

代数级门槛：输入转换 max-abs `<=1e-12 MeV`；分母、`eta` 和 RHS 的
max-abs `<=1e-10 MeV`，且相对 Frobenius 误差 `<=1e-10`。接近零的张量只用
绝对门槛判断。

### C. 短流和完整流

- [x] 在完全相同的初始算符上比较一个 `ds=1e-4` Euler 固定步，并在步后
      重新比较 `H/eta/RHS`。
- [x] 用共同的固定步 RK4 比较 `s=0.001,0.002,0.003` 三个 checkpoint；
      这一步继续隔离 ODE
      driver 差异。
- [ ] 统一直接 flow、生成元、cutoff、停止条件和 ODE 容差后，比较完整流的
      `E/f/Gamma`、`eta` 范数及停止点。
- [ ] 容差缩小十倍后，最终 SR 能量变化 `<1 keV`。
- [ ] 转回普通 `E0+t+V` 后用 double-precision J-coupled readback 比较，
      再用 NCSM 对角化比较谱；float32 `no2bpack` 只做下游格式验收。

### D. 核顺序与完成定义

- [ ] `He4 Nrefmax=0` 全部 A--C 通过。
- [ ] `O16 Nrefmax=0` 全部 A--C 通过。
- [ ] 把命令、commit、输入摘要、误差表和最坏元素写入
      `docs/MR-IMSRG-验收结果.md`。
- [ ] 两核通过后，才恢复把生产 MR 流用于新的 `Nrefmax=2` 物理结果。

“最终能量接近”不能替代 B 阶段；只有 `E/f/Gamma -> denominator -> eta ->
RHS -> flow` 整条链逐层通过，才能称为退化到 `imsrg++`。
