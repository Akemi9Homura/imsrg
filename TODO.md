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
- [x] 比较 `E/f/Gamma`、每个被选中的 EN 分母、`eta(s=0)`、每个
      命名 RHS contraction 与总 RHS，报告 max-abs/Frobenius/最坏元素；
      完整机器表见 `docs/MR-IMSRG-SR-degeneration.json`。
- [x] 代数级输入转换、分母、eta 和 RHS 均通过：两核输入最坏
      `3.55e-15 MeV`，分母最坏 `4.26e-14 MeV`，命名 RHS 最坏
      `1.15e-14 MeV`；短流相对 Frobenius 最坏 `1.24e-15`。完整
      `s=100` 的矩阵元最坏为 `3.00e-8 MeV`，另按 ODE 门限报告。

### F. 相关参考态对 Python m-scheme 的退化门禁

- [x] 先用随机 `Jref=0` 标量 Hamiltonian 和合法 `lambda2` 比较每个
      新 J-scheme MR contraction 展开后的完整 m-scheme 张量。
- [x] 对 `Be8/C12, Nrefmax=0` 与 `He4/O16, Nrefmax=2` 使用同一份
      float64 Hamiltonian/RDM，比较 `E/f/Gamma`、生成元和总 RHS；四体系
      max-abs 上界分别为 `6.51e-11/9.98e-13/7.06e-11 MeV`。
- [x] 把上述四个真实点的总 RHS 继续拆为每个命名 SR contraction 与
      `mr_lambda2_one_body/zero_body`；逐项 J/m 对照的全局 max-abs 为
      `1.954e-14 MeV`（从同一 J 可表示输入展开；原始 m→J 输入误差另列）。
- [x] 为四体系各命名项生成机器可读误差表，补齐相对 Frobenius 与最坏
      J/m 元素索引；40 个命名项/体系项组合的全局相对 Frobenius
      `1.380e-15`，见 `docs/MR-IMSRG-Jscheme-contractions.json`。
- [x] 对四个相关参考测试点逐条比较 White-NCSM 生成元实际选中的 EN
      分母：每核 4 个 1B 与 704 个 2B 正/反向通道，四核 2832 条全部
      通过，全局 max-abs `1.421e-14 MeV`；完整表见
      `docs/MR-IMSRG-Jscheme-denominators.json`。
- [x] 比较共同固定步 RK4 checkpoints：四体系用 `ds=1e-3` 到
      `s=0.001,0.002,0.003`，每点从 driver 的 HO 真空 J64 输出读回、
      HO→NAT 并 MR 正规序后重新比较 `H/eta/RHS`；同 J 输入的全局最坏
      `H/RHS` max-abs 为 `1.279e-13/2.665e-14 MeV`。机器记录见
      `docs/MR-IMSRG-Jscheme-checkpoints.json`。四体系显式
      `H_1=H_0+10^-4 RHS(H_0)` Euler J/m 对照也通过，全局最坏
      `1.776e-15 MeV`；不使用历史上实际调用 RK4 的 `flow_euler` 名称冒充。
- [x] 增加相关参考固定-s 完整流验收器：按生产顺序执行
      HO→NAT→MR flow→反正规序→HO，并比较终点 `H/eta/RHS/vacuum H`；
      He4 Nrefmax=2 本地 `s=0.01` 烟雾流最坏 `4.24e-11 MeV`。
- [x] 统一生成元/cutoff/停止条件/容差后比较 `s=100` 完整流、真空
      `E0+t+V` 与 NCSM 谱；四体系全流 max-abs 为
      `4.40e-9--6.88e-8 MeV`。float64 J64 与 float32 `no2bpack`
      均由独立 NCSM reader 对角化，打包引起的基态差不超过
      `1.225e-6 MeV`；见 `docs/MR-IMSRG-Jscheme-full-flow.json`。

### G. 性能与完成定义

- [x] J-scheme MR RHS 的主存储不随 m-substate 数量增长，主收缩使用块稀疏
      矩阵/现有 channel 索引。真实 `He4 Nrefmax=2` 单线程 RHS 比 Python
      m-scheme 快 `256.8` 倍，输入主数组小 `1468` 倍；解析存储估算可到
      `emax=14` 且不分配稠密 m-scheme 张量。见
      `docs/MR-IMSRG-Jscheme-performance.json`。
- [x] `He4/O16 Nrefmax=0`、`Be8/C12 Nrefmax=0` 和 `He4/O16 Nrefmax=2`
      按 E--F 门禁通过。
- [x] 四核完成流和 NCSM/no2bpack 读回通过，命令、commit、环境、
      误差表、最坏元素及性能写入新验收文档。ODE 容差从 `1e-9` 收紧到
      `1e-10` 后，四核各三条最低 NCSM 能级的最坏变化为
      `1.683e-11 MeV = 1.683e-8 keV`，通过 `<1 keV` 门禁；见
      `docs/MR-IMSRG-Jscheme-full-flow.json`。

### H. 实际较大空间门禁

解析 channel 尺寸和 scratch 上界只证明数据结构没有退回稠密 m-scheme，
不能替代真实较大空间运行。完成 emax2 的 G 门禁后，按以下顺序继续：

- [x] 增加受验证的 reference embedding：把既有 `Nrefmax` 截断 NCSM
      波函数的自然轨道占据、变换和 `lambda2` 原样嵌入更大的单粒子空间；
      新增轨道必须是零占据、零 cumulant、自然变换单位块。读回后重新检查
      粒子数、质子数、Hermiticity、`lambda2` 收缩及低空间逐元素不变，
      并记录新 interaction SHA-256；不得把 embedding 宣称为重新求解过的
      较大 `Nrefmax` 参考态。`MRReference::EmbedInModelSpace/WriteBinary`
      与 `prototype/mrimsrg/embed_jref.py` 已完成；真实 He4 ref2 从 emax2
      嵌入 emax4 后收缩误差 `1.85e-15`、Hermiticity 误差为零，并通过
      完整 J-scheme 文件写读闭环。
- [x] 使用本机已存在的 NNLOopt `hw=20, emax=4, e2max=8` 相互作用，
      对 `He4 Nrefmax=2` 的嵌入参考实际运行一次 C++ J-scheme MR RHS 和
      固定短流；记录 wall time、峰值 RSS、各 profiler contraction、输出
      对称性和 1/64 线程数值一致性。生产运行不得构造完整 m-scheme 张量。
      实测 30 个 J-orbit、34320 条 TBME，单线程 RK4 短流 `1.37 s`、峰值
      RSS `54,360 KiB`；1/64 线程两体最坏差 `4.44e-16 MeV`。V/VI 收缩
      分别占 `0.628/0.492 s`，见 `docs/MR-IMSRG-Jscheme-large-space.json`。
- [ ] 根据 emax4 profiler 优化真实瓶颈；每次优化必须先过 emax2 的 SR、
      相关参考逐项、完整流和 NCSM 回归，再重跑 emax4 门禁。
      第一步已把 VI 的独立 `(J2,w)` 迹提出循环，emax4 VI 从 `0.492 s`
      降至 `0.022 s`，短流从 `1.37 s` 降至 `0.88 s`；输出相对优化前逐元素
      为零，emax2 Python oracle 仍为 `3.11e-12 MeV`，完整相关参考 CTest
      通过。第二步让 V 的 `X/Y/lambda2` 共用同一个 Pandya 求和、复用
      `imsrg++` 的 dense Six-J cache 和双算符 TBME accessor；V 从
      `0.616 s` 降至 `0.332 s`，短流从 `0.88 s` 降至 `0.602 s`，相对最初
      快 `2.28` 倍。单线程 J64 与优化前 bitwise 相同，1/64 线程最坏仍为
      `4.44e-16 MeV`，emax2 oracle 仍为 `3.11e-12 MeV`，完整相关参考 CTest
      用时 `351.93 s` 并通过。块内 OpenMP 在 emax4 无实测收益，已撤回。
      本项只等待 point7 完整流/NCSM 回归后勾选。
- [ ] 只有验证 emax6 interaction 的来源、header 和校验和后，才重复
      emax6 RHS/短流；不得使用文件名含 `candidate` 的核力产生验收结果。
- [ ] 在 emax4/6 实测内存与时间外推支持后，才开始大空间完整流和
      Magnus/重启策略；登录节点禁止重计算，统一由 point7 Slurm 运行。

“最终能量接近”不能替代 E--F；只有
`E/f/Gamma -> denominator -> eta -> named RHS contractions -> flow -> vacuum H -> NCSM`
整条链逐层通过，才能称为 C++ J-scheme 生产 MR 实现完成。

## 已完成基准

- [x] Python m-scheme MR 原型的 RDM、正规序、QCombo/显式 Fock-space
      commutator、White-NCSM、真空物化与 NCSM/no2bpack 闭环。
- [x] Python 生产 MR 路径在 `He4/O16, Nrefmax=0` 下对当前 C++
      `imsrg++` 逐正规序、分母、eta、命名 RHS 收缩、Euler/RK4 和
      `s=100` 完整流严格退化；证据见 `docs/MR-IMSRG-验收结果.md`。
