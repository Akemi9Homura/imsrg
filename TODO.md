# MR-IMSRG 当前 TODO

## P0：生产版 C++ J-scheme MR-IMSRG(2)

当前唯一实现目标是把已验收的 MR 物理放入现有 `imsrg++`
生产 J-scheme 框架。不得为退化验收复制一套 SR commutator、生成元
或 ODE；必须复用现有 `Operator`/`TwoBodyME`/`ModelSpace`、SR
contractions、`Generator` 和 `IMSRGSolver`。Python m-scheme 只作相关参考
态的独立 oracle，不得成为生产 C++ 运行时后端。

固定输入继续使用 NNLOopt、`hw=20 MeV`、`emax=2`、`e2max=4`、NN-only、
`BetaCM=0`、HO 基和 A-dependent 内禀动能。生产范围限定为 `Jref=0`、
`lambda2` 保留、`lambda3=0`。

### A. 文献、约定与架构冻结

- [x] 逐篇阅读 `refs/` 中 Hergert 2016、Gebrerufael/Vobig/Mongelli 学位
      论文和 IM-NCSM 系列论文，建立 m-scheme 公式、J-scheme 张量约定与
      `lambda3=0` 截断表；见 `docs/MR-IMSRG-J-scheme公式与代码映射.md`。
- [ ] 重生并核对 QCombo `MR_IMSRG2.ipynb`，固定指标方向、反对称、
      组合系数和每个命名 contraction。
- [x] 逐文件审读 `Operator`/`TwoBodyME`/`TwoBodyChannel`/`ModelSpace`、
      `Commutator.cc`、`Generator.cc`、`IMSRGSolver.cc` 和现有测试，确定哪些
      SR contractions 在分数占据下已等于 MR 公式的公共部分；公共 2B/1B
      与 0B 项全部复用，新增范围冻结为 IV--VI 与 `C2-lambda2`。
- [x] 冻结 `Operator` TBME 归一化/相位、J-coupled `lambda2` 定义、
      natural-orbit 占据和 `Jref=0` 标量密度格式；用 J-to-m-to-J 随机
      张量闭环验证；独立 oracle 见 `prototype/mrimsrg/jcoupling.py`。

### B. 生产参考态与正规序接口

- [x] 新增最小 `Jref=0` MR reference 对象，保存每个 J-orbit 的自然
      占据与 J-coupled `lambda2`，不把相关密度隐藏成全局状态。
- [x] 读入时验证标量性、Hermiticity、pair antisymmetry、`Tr(gamma1)=A`
      和 `gamma2 -> (A-1)gamma1` 收缩，拒绝 `Jref!=0` 或非自然基密度。
- [x] 复用生产正规序/反变换路径，只增加 `lambda2` 必需的 MR 收缩；
      随机 J-scheme 算符真空往返误差 `<1e-10`。
- [x] 抽取并复用 HF/HFMBPT 的 normalized-pair J-scheme 0/1/2B 基变换，
      用真实 `He4 Nrefmax=2` 非单位自然轨道验证 HO→NAT→HO，并对独立
      m-scheme 四指标变换逐元素比较。
- [x] 将原型验收已有的 `mrimsrg_j64_v1` float64 J-coupled Hamiltonian
      格式纳入 `ReadWrite`，作为 driver 的无损输入/检查点；随机写读误差
      `<1e-13`，既有真实 He4 文件读入误差为零。下游 `no2bpack` 仍保持
      既有 float32 布局。
- [x] 对真实 NNLOopt `He4 Nrefmax=2` Hamiltonian 与 Python 比较自然基
      `E/f/Gamma`，max-abs 分别为 `2.95e-13/3.36e-11/9.95e-12 MeV`；
      `s=0` 最终真空反变换的 0/1/2B 无损输出相对直接 C++ 路径逐元素为零。

### C. J-scheme MR-IMSRG(2) commutator

- [x] 保留并直接调用现有 `comm110/220/111/121/221/122/222_pp_hh/222_ph`
      等 SR 项；禁止复制这些公共公式。
- [x] 只新增形式 MR-IMSRG(2) 在 `lambda3=0` 下剩余的线性
      `lambda2` contractions：1B 的 2B--2B 收缩和 0B 的
      `1/4 C2·lambda2`。
- [x] 每个新 J-scheme contraction 必须同时通过：QCombo/m-scheme 指标式、
      随机标量 J-to-m oracle、一个显式小 Fock-space 矩阵、Hermiticity/
      anti-Hermiticity 和 `lambda2=0` 极限。
- [x] 稀疏块实现必须保留 angular-momentum/parity/isospin 选择定则，不在
      生产 RHS 中展开完整 m-scheme 张量。

### D. White-NCSM 生成元与流

- [x] 在现有 `Generator` 中增加显式 MR/IM-NCSM 模式，复用
      Epstein--Nesbet monopole 分母、`1e-6 MeV` cutoff 和 anti-Hermitization。
- [x] 按 Vobig 6.5.28--6.5.34 与命名表使用舍去 `lambda[2,3,...]`
      的 White-NCSM 分子/分母，MR commutator 本身仍保留 `lambda2`。
- [x] 在球形自然基中实施 `Delta e != 0` 的 1B/2B 掩码，保留同
      HO 量子数参考空间内部耦合。
- [x] 复用 `IMSRGSolver` 直接 flow/ODE 和停止诊断；默认 SR/VS-IMSRG
      不得因 MR 参数或状态而改变。
- [x] 在 `imsrg++` 增加显式 `mr_reference_file` driver：重建分数占据
      ModelSpace、HO→NAT、MR 正规序、直接流、MR 反正规序、NAT→HO，
      并输出 float64 J64 与下游 float32 `no2bpack`；未给参数时不进入 MR。

### E. 单参考生产退化门禁

- [x] 对 `He4/O16, Nrefmax=0`，同一生产 MR 入口在 `lambda2=0`
      时直接运行现有 SR contractions；不允许根据测试核名或 fixture
      切换实现。真实 driver 在 `s=0` 与共同 `ds=1e-4` RK4 后的完整
      vacuum 0/1/2B 均与原生 SR driver 逐元素为零。
- [ ] 比较 `E/f/Gamma`、每个被选中的 EN 分母、`eta(s=0)`、每个
      命名 RHS contraction 与总 RHS，报告 max-abs/Frobenius/最坏元素。
- [ ] 代数级输入转换 max-abs `<=1e-12 MeV`，分母/eta/RHS max-abs
      `<=1e-10 MeV`、相对 Frobenius `<=1e-10`。

### F. 相关参考态对 Python m-scheme 的退化门禁

- [x] 先用随机 `Jref=0` 标量 Hamiltonian 和合法 `lambda2` 比较每个
      新 J-scheme MR contraction 展开后的完整 m-scheme 张量。
- [ ] 对 `Be8/C12, Nrefmax=0` 与 `He4/O16, Nrefmax=2` 使用同一份
      float64 Hamiltonian/RDM，比较 `E/f/Gamma`、生成元、每个 RHS
      contraction 和总 RHS。
- [ ] 比较 `ds=1e-4` Euler 一步与 `s=0.001,0.002,0.003` 共同固定步
      RK4 checkpoints，每点重新比较 `H/eta/RHS`。
- [ ] 统一生成元/cutoff/停止条件/容差后比较完整流、真空
      `E0+t+V` 与 NCSM 谱；float32 `no2bpack` 只做下游格式验收。

### G. 性能与完成定义

- [ ] J-scheme MR RHS 的主存储不随 m-substate 数量增长，主收缩使用块稀疏
      矩阵/现有 channel 索引，并给出相对 Python m-scheme 的时间和内存比较。
- [ ] `He4/O16 Nrefmax=0`、`Be8/C12 Nrefmax=0` 和 `He4/O16 Nrefmax=2`
      按 E--F 门禁通过。
- [ ] 四核完成流和 NCSM/no2bpack 读回通过，命令、commit、环境、
      误差表、最坏元素及性能写入新验收文档。

“最终能量接近”不能替代 E--F；只有
`E/f/Gamma -> denominator -> eta -> named RHS contractions -> flow -> vacuum H -> NCSM`
整条链逐层通过，才能称为 C++ J-scheme 生产 MR 实现完成。

## 已完成基准

- [x] Python m-scheme MR 原型的 RDM、正规序、QCombo/显式 Fock-space
      commutator、White-NCSM、真空物化与 NCSM/no2bpack 闭环。
- [x] Python 生产 MR 路径在 `He4/O16, Nrefmax=0` 下对当前 C++
      `imsrg++` 逐正规序、分母、eta、命名 RHS 收缩、Euler/RK4 和
      `s=100` 完整流严格退化；证据见 `docs/MR-IMSRG-验收结果.md`。
