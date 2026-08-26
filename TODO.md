# MR-IMSRG 当前 TODO

## P1：当前唯一目标——生产 C++ J-scheme MR-Magnus(2)

已验收的 direct-flow MR-IMSRG(2) 保持不变并作为交叉验证。当前只实现
标量 `Omega` 及由它变换出的 Hamiltonian；一般张量观测算符、显式 3N、
`lambda3` 和 `Jref!=0` 不进入本阶段。

### M0. 文献与现有实现冻结

- [x] 复读 Morris--Parzuchowski--Bogner 2015 原始 Magnus 论文、Hergert
      2016 综述、Vobig 2020 Secs. 4.5--4.6/App. B.2 和 Mongelli 2022
      Secs. 5.4--5.6 的原始 PDF；确认 MR-Magnus(2) 的每个嵌套对易子就是
      已验收的 `lambda3=0` MR commutator，当前 `H(s)` 必须由当前
      `Omega` 经 MR-BCH 构造后再生成 `eta(s)`。
- [x] 审计现有 `BCH.cc`、`IMSRGSolver.cc` 和 driver：生产
      `method=magnus` 使用原有 `Solve_magnus_euler()`，以
      `BCH_Product(ds*eta,Omega)` 累积 `Omega`，以
      `BCH_Transform(Hstart,Omega)` 得到当前 Hamiltonian。缺口冻结为
      BCH/Product/NewOmega/GatherOmega/Transform 的 MR dispatcher 和
      driver scope，不另写 ODE 或复制 SR Magnus。
- [x] 在实现计划中冻结现有代码的符号/指数乘法顺序、分段 `Omega` 语义、
      Bernoulli/BCH 阈值、允许/拒绝的 correction 开关及验收容差。

### M1. Python m-scheme Magnus oracle

- [x] 基于 `prototype/mrimsrg/commutator.py` 实现最小 scalar MR-BCH：
      对固定 `Omega,O` 保存每阶 `ad_Omega^k(O)/k!`，检查 Hamiltonian
      Hermiticity、`Omega/eta` anti-Hermiticity 和逐阶收敛。
- [x] 实现与生产 `BCH_Product(dOmega,Omega)` 相同乘法顺序和截断的 oracle；
      固定随机小模型逐阶验收 Bernoulli 项，且在 `Omega=0` 时严格满足
      `dOmega/ds=eta`。
- [x] 冻结 `He4/O16 Nrefmax=0` 以及 `Be8/C12 Nrefmax=0`、
      `He4/O16 Nrefmax=2` 的首步、短流 `Omega`、逐阶 BCH 和最终 MR 正规序
      Hamiltonian fixture；不只保存零体能量。

### M2. 生产 C++ 接入

- [x] 为现有 scalar `BCH_Transform` 与 `BCH_Product` 增加显式 MR reference
      上下文；未提供 MR reference 时必须逐调用保持原 SR 路径，MR 时每个
      嵌套对易子必须调用 `MRCommutator::Commutator`。
- [x] 将 `Solve_magnus_euler`、`NewOmega`、`GatherOmega`、Hamiltonian
      重建与最终 `Transform` 接到同一 dispatcher；删除 driver 对 MR+Magnus
      的拒绝，同时继续拒绝 MR correction/IMSRG(3)/一般流动算符等未验收组合。
- [x] 保持现有 Magnus 的 `ds_0/dsmax/ode_tolerance/omega_norm_max`、分段
      与 checkpoint 语义，不创造 MR 专用 ODE 参数；原一阶 `method=magnus`
      保留为有限步 BCH-product oracle。
- [x] 将 `gen_job.py --mr-jscheme` 的唯一生产入口切换为
      `method=magnus_adaptive`，启用已有 `write_omega`/scratch 输出并在 manifest 中
      明确记录为本段 `Omega`；生成器不暴露或覆盖 MR 专用步长参数，driver
      回归测试同时验证非空 `Omega` 文件确实写出。
- [x] 修复已有但未完成的 `magnus_adaptive`：按 Vobig Eq. (4.6.4) 逐阶计算
      完整 Bernoulli 导数，以 MR dispatcher 计算每个嵌套对易子，并用现有
      `ode_tolerance/dsmax/omega_norm_max` 做 RK45 与分段；原一阶
      `method=magnus` 保留为独立 Euler oracle，不用极小 `domega` 掩盖误差。

出处注：完整 Bernoulli ODE/BCH 见 Morris--Parzuchowski--Bogner 2015
Eqs. (28),(30)；Magnus(2) 截断、级数停止和自动步长 RKF45 见
Vobig 2020 Eqs. (4.6.4),(4.6.7),(4.6.8) 及 App. B.2。完整书目信息、
本地 PDF 与“文献事实/仓库工程决策”边界见实现计划 Sec. 8.1。

### M3. 分层验收

- [x] 随机 scalar J-scheme 输入展开到 m-scheme，逐阶比较 MR-BCH 与
      BCH product 的 0B/1B/2B，max-abs `<=1e-10 MeV`。
- [x] `lambda2=0` 的 MR 生产入口对 `He4/O16 Nrefmax=0` 逐元素退化到现有
      SR Magnus：每段 `Omega`、每阶 BCH、`eta` 和最终 Hamiltonian 均比较，
      不允许核特例或测试专用分支。
- [x] 六个相关/单参考固定体系的首步 `Omega`、零阈值 MR-BCH Hamiltonian
      和默认阈值生产 driver 已逐项复现 Python/C++ oracle。
- [x] Magnus 短流与 direct flow 已在同一非零 `lambda2`、有隙小模型中按
      `ds=1e-2,5e-3,2.5e-3` 比较，Hamiltonian 差为
      `1.338e-4/3.368e-5/8.447e-6`，比值 `0.252/0.251`，验证共同连续流
      极限和预期二阶局部差。
      Magnus(2) 与
      direct-flow IMSRG(2) 长流不要求 bitwise 相同，但差异必须随数值阈值
      稳定并有矩阵元级记录。
- [x] `He4 Nrefmax=2, s=1` 将 `dsmax/ode_tolerance` 与 BCH 阈值十倍收紧后，
      Hamiltonian 最大矩阵元变化 `0.2708 keV`，NCSM `Nmax=8` 三条低能级
      最大变化 `0.4668 keV`；所有分段 `Omega` anti-Hermitian、Hamiltonian
      Hermitian，真空物化/J64/no2bpack 均读回，float32 packing 误差最多
      `1.09 eV`。
- [x] 完整 Release CTest `21/21` 以 `411.32 s` 通过；ASan+UBSan 的
      `MRMagnus/MRSRDriver/MRJobGenerator` 以 `10.03 s` 通过。论文
      Bernoulli 停止条件下，emax4/6 的 `s=1e-4` 单线程生产短流分别为
      `0.64/5.59 s`、峰值 `61.5/390.1 MiB`，均写出 J64；相对固定做到
      `k=8` 的输出逐位相同，scalar commutator 数从 `65` 降到 `23`。
- [x] 将独立 Python m-scheme 原型的标准库 `unittest discover` 注册为
      CTest `MRPrototype`；54 项 RDM/正规序/commutator/J-coupling/
      generator/flow/Magnus/I/O 测试全部通过，不再依赖手工另跑 README 命令。

### M4. 生产收尾：级数失败不得静默接受

- [x] 将 Vobig Eq. (4.6.8) 从“仅禁止提前停止”升级为可诊断的
      Magnus 级数失败：生产 `MagnusDerivative` 必须区分已收敛、
      嵌套范数不再递减与达到最大阶仍未收敛；禁止固定
      `k=8` 后静默返回不可靠的导数。`relative_threshold=0` 仍保留
      固定 `k<=8` 的代数 oracle 语义。
- [x] `magnus_adaptive` 在 RK45 stage 内遇到级数失败时必须拒绝该
      trial step、缩小步长并从未修改的已接受状态重试；若起点
      `Omega` 本身已无法收敛，或同一段累计八次 stage 缩步仍不能
      安全前进，则先物化当前段并以 `Omega=0` 开新段，
      禁止无限重试。日志记录失败阶数、范数比、旧/新步长及分段。
- [x] 用可控小模型强制触发“stage 拒步后恢复”和“起点分段后
      恢复”，验证拒绝步不改写 `s/Omega/H`、恢复后与小步长 oracle
      一致；重跑随机 MR-BCH、He4/O16 SR 退化、sanitizer 和分段
      `Omega` 读写测试。独立 Debug
      `-fsanitize=address,undefined` 构建下
      `MRMagnus/MRSRDriver/MRJobGenerator` 均通过；完整 Release
      CTest `21/21` 以 `499.20 s` 通过。
- [ ] 按 `He4 -> Be8 -> C12 -> O16` 跑完整生产流；每个体系保存
      默认/十倍收紧的最终 0B/1B/2B 矩阵元差、解耦残差、
      Magnus 收敛/拒步计数、分段 `Omega` 及 NCSM 低能谱；不以零体项
      代替读回对角化。默认 profile 不写 ODE 参数；
      `--mr-tight-validation` 只用于独立验收复算，显式将既有
      `ode_tolerance=1e-6` 收紧为 `1e-7`，不改 `ds_0/dsmax`。
      tight 必须把 `target_s` 固定为 default 的实际停止点并用
      `eta_criterion<=1e-20` 禁用提前停止，禁止比较两个不同流时刻。
      单核 gate 直接从 metadata 强制上述 solver/profile/override/终点条件，
      不把生成器单元测试当成生产结果确实使用这些设置的替代证据。
      四核聚合只接受包含完整 13 项门的 `mrimsrg_magnus_gate_v2`，拒绝
      旧 schema 或仅手填顶层 `passed=true` 的报告。
      flow reader 按 `IMSRGSolver::WriteFlowStatus` 的固定列宽读取到 `Eta_3`，
      再读取完整十进制 `Ncomm`；不能用空白切分，也不能把 `setw(7)` 当成
      最大宽度，因为七位计数会紧接 `Eta_3`，八位计数还会继续扩展。
      每份 v2 JSON 保存 gate reader 自身 SHA-256，四核聚合要求一致。
- [x] `Be8 Nrefmax=0` 与 `C12 Nrefmax=0` 的完整 default/同终点
      tight Magnus 流和 NCSM 读回已通过。二者 default 分别在
      `s=330.08347/226.53612` 达到残差比
      `9.78935e-7/9.71606e-7`；十倍收紧的三态最大变化为
      `0.64841/0.76990 keV`，J64/no2bpack 最大 packing 差为
      `1.145/3.305 eV`。作业为 `100602/100604/100618/100616`，
      拒步计数均为零；point7 原生 gate JSON 已写入被忽略的结果目录
      `result/mr-jscheme-magnus/gates/`。加入落盘 Omega 数值检查后的
      Release 完整回归 `22/22` 以 `491.61 s` 通过；gate 同时要求 default/tight
      均实际写出至少一个非空 `Omega` 段，并由链接生产 `Operator` 的
      C++ checker 读回检查标量量子数、有限性及 0B/1B/2B
      anti-Hermiticity；测试同时确认人为破坏的对角 1B 元素会被拒绝，
      Be8/C12 真实生产段的独立 ASan+UBSan 读回也已通过。
- [ ] `He4/O16 Nrefmax=2` default 长流 jobs `100600/100606` 已启动并
      通过输入、参考态和 `Omega` 写出检查；完成后以各自实际
      停止 `s` 生成同终点 tight 复算。
- [ ] 四核均必须达到目标掩码内解耦残差相对初值 `<1e-6`，
      十倍收紧后 NCSM 验收能级变化 `<1 keV`，并保持 Hermiticity、
      anti-Hermiticity、J64/no2bpack 读回窗和已验收的 m-scheme/SR 极限；
      最后由 `aggregate_magnus_gates.py` 强制核种/Nrefmax、interaction、
      NCSM Nmax/三态、冻结的 production executable/library SHA-256、
      NCSM/Omega validator
      与正式阈值一致，同时逐核锁死 JREF 和 A-dependent bare J64 哈希并
      生成单一四核 JSON，不手工拼接完成结论。

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
- [x] 重生并核对 QCombo `MR_IMSRG2.ipynb`，固定指标方向、反对称、
      组合系数和每个命名 contraction。QCombo 0.2.0 / SymPy 1.14.0
      的 17 个 LaTeX display 与保存输出逐字一致；1B `lambda2` 八项、
      0B `1/4 C2 lambda2` 和无显式 `lambda2` 的 2B 方程均已逐项映射到
      Python oracle 与 C++ IV--VI，见公式映射文档 Sec. 3.1。
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
- [x] 根据 emax4 profiler 优化真实瓶颈；每次优化必须先过 emax2 的 SR、
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
      point7 `s=100, rtol=atol=1e-10` 完整流回归 jobs
      `100377/100379/100381/100383` 已通过；新旧真空 Hamiltonian 全局最坏
      `7.44e-13 MeV`，四核各三条 NCSM 低能谱全局最坏变化
      `2.13e-13 MeV`，本项验收完成。
- [x] 进一步把 V 的 Pandya `lambda2` 中间指标限制到 cumulant 的精确轨道
      支撑域。若参考态在全空间活跃，该实现自动回到原全矩阵乘法；若固定
      emax2 参考嵌入较大流空间，则只形成
      `X(:,A) Lambda(A,A) Y(A,:)`，不引入阈值或近似。emax4/emax6 单线程
      短流均与优化前 J64 逐字节相同；emax4 总时由 `0.602 s` 降至
      `0.408 s`，emax6 由 `11.70 s` 降至 `6.81 s`，其中 emax6 V 从
      `7.751 s` 降至 `2.851 s`、BLAS 从 `1.873 s` 降至 `0.072 s`。
      emax6 的 1/64 线程两体最坏差仍为 `8.88e-16 MeV`，emax2 Python
      oracle 仍为 `3.106e-12 MeV`；随机全活跃参考及四核 checkpoint 的
      `MRCorrelatedDriver` 用时 `350.49 s` 并通过。实现 commit
      `cef4fce6`。
- [x] 在活跃支撑域上继续只构造 V 真正进入乘法的
      `X/Y(:,A)` 与 `X/Y(A,:)`，不再先生成完整 `d*d` Pandya X/Y 后切片；
      全空间活跃参考保留原完整分支。emax4/emax6 短流相对 `cef4fce6`
      仍逐字节相同，总时分别 `0.408 -> 0.312 s`、`6.81 -> 4.15 s`；
      emax6 V/build 为 `0.170/0.111 s`，相对上一版快 `16.76/25.00` 倍。
      emax2 Python oracle 为 `3.106e-12 MeV`，随机全活跃参考与四核
      checkpoint 的 `MRCorrelatedDriver` 用时 `352.52 s` 并通过。实现
      commit `dea9125b`。当前 emax6 主瓶颈已转为既有 White-NCSM generator
      (`2.191 s`)，不是 MR `lambda2` addon (`0.717 s`)。
- [x] 只有验证 emax6 interaction 的来源、header 和校验和后，才重复
      emax6 RHS/短流；不得使用文件名含 `candidate` 的核力产生验收结果。
      输入前置现已通过：从冻结 emax14 母文件
      `118c4c27...b9d7d` 流式抽取 467032 条原始 records，得到 canonical
      emax6 文件 SHA-256 `199d5a7a...b913b`；emax4→emax2 校准逐字节相同，
      emax14→emax6 在 A=4/16 下逐 J-channel 严格零差。旧 candidate 最坏
      差 `3.8147e-6 MeV`，继续禁用。He4 Nrefmax=2 嵌入参考的真实 emax6
      固定 RK4 一步已完成：56 个 J-orbit、467032 条 TBME；活跃支撑域优化后
      单线程 `6.81 s`、峰值 RSS `280880 KiB`，1/64 线程两体最坏差
      `8.88e-16 MeV`，且未构造 dense m-scheme Hamiltonian。
- [x] 对 `cef4fce6` 的活跃支撑域优化重新执行 point7 四核
      `s=100, rtol=atol=1e-10` 完整流和 NCSM 谱回归。generator jobs
      `100385/100387/100389/100391` 全部通过；四核 C++/Python 全流最坏差
      为 `2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`。新旧 J64 的全局
      max-abs 为 `8.384e-13 MeV`，`Nmax=8/0/0/2` 的 12 条 simpleFCI
      低能谱全局最坏变化 `1.741e-13 MeV`。
- [x] 对 `dea9125b` 的 skinny Pandya 构造重新执行 point7 四核
      `s=100, rtol=atol=1e-10` 和 `Nmax=8/0/0/2` simpleFCI 谱回归。
      generator jobs `100393/100395/100397/100399` 全部通过；四核
      C++/Python 全流最坏差为
      `2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`。相对
      `cef4fce6` 的 J64 全局 max-abs 为 `8.527e-13 MeV`，12 条低能谱
      全局最坏变化 `1.990e-13 MeV`，因此 skinny 构造已通过完整流、
      真空物化和后 NCSM 三层门禁。
- [x] 复用 White-NCSM 分母中每个 RHS 内不变的有序轨道对 monopole，
      避免每个 2B bra/ket、正反方向重复跨 J 求和；保持 Vobig
      Eqs. (6.5.31--34) 的算术次序、`Delta e` 掩码、权重和 cutoff 不变。
      commit `6e02676c` 的 emax4/emax6 短流均与 `dea9125b` 逐字节相同，
      generator 分别由 `0.1112 -> 0.00324 s`、`2.191 -> 0.0426 s`，总时
      由 `0.312 -> 0.222 s`、`4.15 -> 1.98 s`。分母逐项、SR/MR driver、
      emax2 Python oracle (`3.106e-12 MeV`) 和随机全活跃参考 CTest
      (`354.53 s`) 均通过。point7 generator jobs
      `100401/100403/100405/100407` 完成四核 full-flow；新旧 J64 全局
      max-abs `9.308e-13 MeV`，12 条低能谱最大变化 `1.137e-13 MeV`。
- [x] 在当前 Python/QCombo 环境重新执行
      `refs/qcombo/examples/MR_IMSRG2.ipynb`，把重生的 0B/1B/2B 项数、
      系数和指标序与 `MRCommutator.cc` 的命名项逐项记录。规范化 LaTeX
      输出 SHA-256 为 `431b893e...1a629a57`，完整环境、项数和映射见
      `docs/MR-IMSRG-J-scheme公式与代码映射.md` Sec. 3.1。
- [x] 在根目录生产 `gen_job.py` 增加显式、一次只生成一个参数点的
      `--mr-jscheme` 模式；冻结 C++ direct-flow 参数、输入/reference
      SHA-256、累计 `start_s/target_s`、该段 `smax`、J64/no2bpack 输出和
      point7 Slurm 环境。既有 SR/VS 参数路径不得改变。生成器另冻结
      git commit 与 executable SHA-256、拒绝覆盖，并由 `MRJobGenerator`
      CTest 和真实 emax4 input/reference 的 `--generate-only` 检查通过。
- [x] 用同一 emax2/4 Hamiltonian 验证 direct flow 的连续运行与 J64
      分段重启：比较终点完整 vacuum 0B/1B/2B、流诊断和输出校验和；只有
      数值等价通过后，J64 才能称为 Hamiltonian checkpoint/restart。
      `0→0.02` 对 `0→0.01→0.02` 在 emax2/4 的完整 0/1/2B 最坏差为
      `6.15e-14/5.92e-14 MeV`，终点 `E/||H||/||eta1,2||` 在打印精度相同；
      可复现比较器和机器记录见 large-space JSON。
- [x] 用 `gen_job.py --mr-jscheme` 在 point7 提交 `He4 Nrefmax=2, emax=4`
      direct flow；两段均先检查生成脚本并通过 `sbatch --test-only`。job
      `100413` 完成累计 `s=0→100`（`159.38 s`, `166.99 MB`），job
      `100423` 以 lossless J64 继续到累计 `s=1000`（`187.56 s`,
      `167.19 MB`）。两段均正常结束、MR 反正规序零体核对在打印精度为零，
      且同时产生 J64/no2bpack；但没有达到 `eta_criterion=1e-6`，所以只把
      它们称为已完成的实测流段，不能称为完整收敛流。
- [x] 对 `s=0/100/1000` 的每个 `Delta e != 0` J-scheme 通道分解
      Hamiltonian、lambda-free White-NCSM 分子、正反方向占据权重、EN
      分母和实际 `eta`。新增可复现工具
      `prototype/mrimsrg/diagnose_white_ncsm_tail.py`。`Rgen/Rgen(0)` 从
      `1` 降到 `1.0265e-2/1.3102e-3`，`Rnum/Rnum(0)` 降到
      `3.4784e-3/4.4032e-4`；主导分母约 `24.1 MeV`，没有小分母或 cutoff
      通道。慢尾来自 `2e-4--2e-2` 的小自然占据权重使不同通道依次成为
      主导，而不是 ODE 容差或 generator 计算性能。
- [x] 扫描 emax4 的下游谱稳定流窗，而不是只看 MR 零体项或 `eta`。
      `s=0,0.02,0.1` 已完成 `Nmax=2,4,6,8,16` NCSM；`s=1,10,100,1000`
      已完成 `Nmax=8`。裸 Hamiltonian 的全空间 (`Nmax=16`) 基态为
      `-25.2912232776 MeV`；`s=0.02/0.1` 的全空间漂移为
      `-45.476/-187.528 keV`。`s=100/1000` 的 Nmax8 基态已漂到
      `-30.588/-42.621 MeV`，因此明确拒绝作为下游生产 Hamiltonian。
- [x] 对 emax4 `s=0.02` 的 J64 与 no2bpack 做独立 NCSM 读回。前三态
      float32 packing 差值最大 `0.000726 keV`；ODE 容差由 `1e-9` 收紧到
      `1e-10` 后，J64 0/1/2B 最坏变化 `3.675e-14 MeV`、前三态最坏变化
      约 `3e-11 keV`，通过 `<1 keV` 门槛。
- [x] 冻结当前 emax4 判定：`s=0.02` 只称为“有限短流 NCSM 收敛加速器
      候选”。它把 `Nmax=2/4/6/8` 相对同一流后全空间的截断误差分别改善
      `412.148/147.117/39.306/9.869 keV`，但
      `Rgen/Rgen0=0.9659`，远未通过严格 `1e-6` 脱耦门槛。严格脱耦结果
      仍记为未获得，不得用短流候选替代。
- [x] 在进入 emax6 生产流前冻结 `prototype_downstream_stability_v1`：选择
      通过全部门槛的最大已测 `s`，要求全空间基态漂移 `<=100 keV`、每个
      已测 Nmax 的同终点截断误差均不得变差、packing/ODE 误差各 `<1 keV`
      且生成元 norm 必须下降。100 keV 是快速原型的显式工程预算，不冒充
      文献阈值；Vobig Sec. 6.5.5 只提供用后 NCSM 平台识别诱导高体项失控
      的方法依据。可复现选择器和 CTest 为
      `select_downstream_flow_window.py`/`MRDownstreamWindow`，当前唯一通过点
      是 `s=0.02`。严格 `Rgen/Rgen0<=1e-6` 继续作为独立失败门禁。
- [x] 只在上述 finite-s rule 冻结后运行 emax6 短流窗。canonical
      NNLOopt `emax=6,e2max=12` 输入 SHA-256 为
      `199d5a7a...b913b`；嵌入的 He4 `Nrefmax=2` 参考保持原 emax2
      占据/`lambda2`，新增轨道严格为空。point7 job `100432` 已完成
      `s=0->0.02` direct flow（`10.48 s`, `886.6 MB`），实际
      `||eta||/||eta(0)||=0.96099`；它是 finite-s pilot，不是严格脱耦。
- [x] emax6 的 lossless J64、独立 no2bpack exporter 和生产 no2bpack
      已通过 `Nmax=8` 三态谱闭环；两条 no2b 写路径的三态能量在
      double 数值噪声内相同。job `100434` 把 ODE 容差从 `1e-9`
      收紧到 `1e-10`，J64 0/1/2B 最坏变化
      `1.18e-13/1.13e-13/4.16e-13 MeV`，三态最坏变化约
      `1.2e-10 keV`，通过 `<1 keV` 门禁。
- [x] 完成 emax6 最大可行 NCSM 代理门禁。裸 Hamiltonian 与
      `s=0.02` 的 `Nmax=2,4,6,8` 已完成；流后相邻 Nmax gap 均小于
      裸值。point7 Nmax10 jobs `100441/100442` 均以 exit `0:0` 完成，
      维数 `183866`，裸/流后基态为
      `-26.790750086/-26.860681522 MeV`；`Nmax8->10` gap 从
      `618.748` 降为 `585.962 keV`，Nmax10 流/裸差为 `-69.931 keV`，
      因此预先冻结的最大可算空间代理通过。两作业各约 `174--175 s`、
      峰值 RSS `438--439 MiB`。这不是 emax6 全空间 (`Nmax=24`) 漂移
      门禁，不能冒充 `prototype_downstream_stability_v1` 的严格通过。
- [x] 根据上述实测冻结 emax6 `s=0.02` finite-s pilot、完整三态谱、
      ODE/格式误差、输入/可执行文件/作业哈希和 claim limit；机器记录见
      `MR-IMSRG-Jscheme-large-space.json` schema v10。
- [x] 建立并独立核对 emax8 canonical interaction：从冻结 emax14 母
      minipack 流式保留 3526624 条原始 records，child SHA-256 为
      `3fd1a003...d4bda`；A=4/16 的逐 channel 比较均严格零差。对应 bare
      J64/ref/no2bpack 已冻结校验和，嵌入参考的收缩误差 `1.85e-15`。
- [x] 完成 emax8 真实资源门禁。point7 job `100444` 对 90 个 J-orbit、
      3526624 条 TBME 执行固定 `ds=1e-4` RK4 一步（4 次 RHS），墙钟
      `17.282 s`、峰值 RSS `2070.297 MB`，真空物化零体核对为零；全程只走
      J-scheme。稠密 m-scheme 两体张量解析需求约 `1.518 TB`，lossless
      J64 validator 因该旧路径 `bad_alloc`，明确不作为 emax8 验收路径。
- [x] 用隔离的 `simple-ncsm` native no2bpack reader 完成 emax8 `Nmax=8`
      资源/格式 smoke test；该简化求解器不作为正式 NCSM 谱验收。
      bare/固定一步流 jobs `100447/100448` 均 `COMPLETED 0:0`，维数
      `44838`，三态谱分别为
      `[-26.2170529924,-24.7096273755,-24.7096273017]` 与
      `[-26.2177863230,-24.7099878250,-24.7099679677] MeV`；墙钟
      `81/75 s`、峰值 RSS `2385856/2386220 KiB`。基态变化
      `-0.733331 keV` 只证明生产物化/读回与资源可行，不是有限流物理门禁。
- [x] 解决 `simple-ncsm` 的大轨道表占据分区未定义行为，并用同一
      Hamiltonian 与 BIGSTICK 做严格小空间闭环。原先把 90 个 J-orbit
      occupation 塞进 `uint64_t` 会产生 `shift exponent 65`；隔离分支现用
      `std::vector<int>` collision-free key。原 emax8、660 m-orbit、
      `Nmax=2` 故障输入已正常 `exit 0` 并逐位复现三态。随后冻结 He4
      N2LOopt A=4、hw20、emax2/e2max4、BetaCM=0 的同一 no2bpack
      (`SHA-256 5cadb860...e01de`) 给两个求解器；二者行列式和完整谱均为
      59 维/59 态，最大谱差 `4.165e-6 MeV`，BIGSTICK 波函数在
      simple-ncsm Hamiltonian 下最大残差 `4.137e-6 MeV`。检查器随
      `simple-ncsm` commit `fcc6070` 提交；小空间 NCSM 问题已关闭，但
      大空间生产谱仍由 BIGSTICK 正式验收。
- [x] 用根目录单点生成器实际运行 emax8 `s=0.02` direct flow，并以
      `rtol=atol=1e-9/1e-10` 独立复算。point7 jobs `100451/100452` 均
      `COMPLETED 0:0`；25/31 次 scalar commutator 用时 `67.67/85.04 s`、
      峰值 RSS `5194.719/5227.641 MB`。lossless J64 的 0B/1B/2B 最坏差为
      `4.80e-14/6.53e-14/1.90e-13 MeV`。native no2bpack Nmax8 lightweight
      jobs `100458/100456` 的三态最大差 `4.97e-11 keV`，只作为格式/ODE
      smoke test；正式谱门禁改由 BIGSTICK 重做。
      初次松容差读回 job `100454` 因未知非零子进程返回码失败，虽打印相同
      谱也不计入成功；正式结果来自 `max_iter=500` 的成功重跑 `100458`。
- [~] **按用户决定终止，不执行：** 不再对 emax8 bare/`s=0.02` 做
      BIGSTICK Nmax 序列，不再研究 finite-s IM-NCSM 收敛改善或长流谱漂移。
      已有结果只作为已知边界保留，不再据此选择流窗、生成元、ODE 或
      Magnus。当前任务重新收敛到 C++ J-scheme MR-IMSRG 本身的实现、验证
      与性能优化。
- [x] **MR-P0：干净构建与生产测试矩阵。** 从当前提交新建 Release 和
      sanitizer build，运行全部 MRReference、正规序往返、命名
      `lambda2` contraction、随机 J↔m oracle、White-NCSM denominator/
      generator、SR dispatcher 退化、六个真实参考态 RHS/checkpoint 与
      driver 测试。必须区分已有覆盖与实际缺口；不得用后 NCSM 能量替代
      `E/f/Gamma -> eta -> named RHS` 的逐层断言。
      干净 Release build 的 8 个 MR CTest 已全部通过，总用时 `403.99 s`；
      优化后的完整四核相关参考逐命名 contraction/checkpoint oracle 再次以
      `374.41 s` 通过。ASan+UBSan 首轮在 MR 必经的既有 SR
      `comm220ss` 空 channel 上发现 Armadillo 空指针引用，已由 commit
      `5dfa32a4` 加严格零贡献保护；修复后 sanitizer 下 MR driver、SR
      退化和 denominator 三项通过；本轮两个优化后的完整四核相关参考
      sanitizer oracle 也以 `544.38 s` 通过。剩余的 `WriteBinary` 重复写
      异常已最小复现为 host 环境问题：纯 C `/usr/bin/python3` 启动时尚未
      加载 libstdc++，预加载 ASan 因而没有 `__cxa_throw` 的真实入口。
      commit `04202213` 让 sanitizer CTest 保持 libasan 第一并同时预加载
      libstdc++，还强制 `PYTHONPATH/LD_LIBRARY_PATH` 指向当前 build，避免
      误载 Release 模块。只对 Python-hosted 测试关闭非 instrumented
      CPython 的进程退出 LeakSanitizer；AddressSanitizer/UBSan 仍为
      `halt_on_error=1`。正式 sanitizer `MRReference`（含重复写拒绝）用时
      `15.51 s` 通过，六项短门禁 `46.85 s` 全通过，完整四体系
      `MRCorrelatedDriver` 又以 `465.57 s` 通过；普通 build 七项回归
      `41.40 s` 全通过。因此 MR-P0 的数值、生产入口和 sanitizer 环境均
      已闭合。
- [x] **MR-P1：补齐生产入口验收缺口。** 审计 `imsrg++ -> IMSRGSolver ->
      MRCommutator/Generator` 的真实调用链；对任何只在 helper/pybind 层测试、
      未覆盖生产 dispatcher 的路径增加最小回归。`lambda2=0` 必须逐元素
      退化到现有 SR，相关参考态必须逐命名 contraction 退化到 Python/QCombo
      oracle；优化前先冻结基准输出与容差。
      审计确认真实可执行程序已强制 HO 基、`white-ncsm`、direct flow、
      单一步空间、NN/NO2B、无额外流动算符，并已有相关参考 driver、
      `He4/O16 lambda2=0` 逐元素 SR 退化、四个相关参考逐命名 contraction
      和连续 RK4 checkpoint 门禁；没有发现只在 helper 层成立的核心路径。
      sanitizer 暴露并修复的 `comm220ss` 空 channel UB 是本轮唯一实际生产
      调用链缺陷。
- [x] **MR-P2：只剖析和优化 MR-IMSRG RHS。** 在 emax2/4/6/8 固定输入上
      测量单次完整 `MRCommutator + White-NCSM generator` 的 profiler、峰值
      RSS、临时张量尺寸和 1/多线程复现误差，定位当前真实瓶颈后复用
      `imsrg++` 的 channel/block/cache 基础设施做等价优化。禁止以改变物理
      掩码、舍弃小 `lambda2`、降低空间或缩短流来换性能。
      emax6 单线程基线（commit `5dfa32a4`，4 次 RHS）为 profiler real
      `2.14910 s`、峰值 `281084 KiB`，其中 MR setup/VI/addon 为
      `0.26468/0.26566/0.77661 s`。commit `e3965232` 仅用 `lambda2` 的
      精确标准耦合 pair 支撑计算 `lambda*X/Y`，setup 降到 `0.03311 s`；
      commit `6a7227c4` 再跳过 VI 中严格为零的 trace 元素，VI/addon/real
      降到 `0.14228/0.43116/1.77368 s`。优化前后 emax4/6 流后 J64 均
      bitwise 相同；emax4 最终 real/MR addon/峰值为
      `0.19801 s/0.06697 s/54128 KiB`，emax6 的 1/16 线程输出也 bitwise
      相同。
      当前 emax8 单线程实测 90 个 J-orbit、3526624 TBME、profiler real
      `12.33056 s`、峰值 `1377824 KiB`；MR addon `1.98985 s`，而既有 SR
      commutator `5.67302 s`（`comm222_phss=4.13714 s`）已是主瓶颈。
      emax2 对应 real/MR addon/峰值为
      `0.01889 s/0.00997 s/12484 KiB`。
      commit `b0c832a8` 进一步把标准耦合 IV/VI 共用的
      `lambda2*Y/lambda2*X` 从完整 channel 矩阵改为 cumulant 精确 pair
      支撑上的活跃行存储；VI 用保持原交换相位和 normalized-pair
      `sqrt(2)` 因子的 accessor 读回，不使用数值阈值，也不改变全活跃
      参考态的原矩阵乘法分支。新增 profiler 统一记录显式 standard-product
      存储和避免的 dense `LY/LX` 字节：emax4/6/8 分别避免
      `907/14067/108943 KiB`，对应新 product 存储为
      `1198/14949/110935 KiB`。三个空间的流后 J64 相对冻结基线均
      bitwise 相同，emax6 的 1/16 线程 SHA-256 也完全相同；emax8 当前
      MR addon 为 `1.91721 s`，其中 IV/V/VI/setup 为
      `0.25080/0.69959/0.77368/0.14675 s`。总 RSS 仍由完整算符和既有
      SR Pandya 路径主导，不能把消除 106 MiB MR 临时量误写成同等幅度的
      进程峰值下降。下一步只继续审计 MR IV/VI 的 channel 遍历、缓存和
      输出临时量；共享 SR `comm222_phss` 虽是 emax8 总体主瓶颈，但不在
      没有独立 SR/VS 性能与逐元素回归时贸然改动。
      commit `21760bbf` 完成上述 MR 查询审计：IV 对 `XLY/YLX`、VI trace
      对 `LY/LX` 各只做一次 channel/ket/交换相位/normalized-pair 查询；
      VI 主和复用现有 `GetTBME_J_twoOps()` 同时读取 `X/Y`，并在 accessor
      前只跳过严格违反角动量三角条件或不属于同一 scalar one-body channel
      的零贡献。有效项次序与算术表达式不变。emax4 的 IV/VI/addon 从
      `0.00860/0.01415/0.06659 s` 降到
      `0.00404/0.00562/0.05402 s`；emax6 从
      `0.07247/0.15042/0.46645 s` 降到
      `0.02706/0.04122/0.25384 s`；emax8 从
      `0.25080/0.77368/1.91721 s` 降到
      `0.11566/0.19998/1.12244 s`。emax8 单次短流 wall 由本轮相邻实测
      `13.86 s` 降到 `11.36 s`，峰值 RSS 保持约 `1.35 GiB`。当前 MR
      最大单项已变为 V/build；下一轮先审计 ordered-cross-block 构造是否
      仍重复生成精确相同的 active-active Pandya 元素，再决定最小改动。
      审计中试验的 active-active 复制和块内 OpenMP 在 emax8 均无可测
      收益，且后者增加约 30 MiB 峰值，已全部撤回。真正的浪费是 partial-
      active V 先形成两个完整 `d*d` 输出、最终却只读取共同 spectator
      `t` 的少量元素。commit `287adea4` 保留全活跃分支原 BLAS，只保存
      `X*lambda/Y*lambda` 的 `d*a` skinny 左因子并按实际迹元素点积。
      emax4/6/8 每个最大块分别避免 `441/4390/26244 KiB` dense 输出；
      emax6 的 V/BLAS/addon 从 `0.16431/0.04214/0.25384 s` 降至
      `0.13445/0.00283/0.22717 s`，emax8 从
      `0.62761/0.25488/1.12244 s` 降至
      `0.34349/0.00590/0.86220 s`。emax8 单线程 wall
      `11.36 -> 10.95 s`，16 线程 wall `7.27 -> 6.95 s`；进程 RSS 仍由
      完整算符/SR scratch 主导。下一步对 IV 的标准耦合 `XLY/YLX` 做同类
      迹专用审计；它们当前仍是 standard products 中最大的 dense 输出。
      commit `460a009d` 已完成该审计：空 cumulant channel 不再分配零
      `d*d` 矩阵，partial-active channel 保存 `X/Y(:,A)` 与
      `lambda*Y/X(A,:)` 并只点积实际 IV 迹元素；全活跃 channel 保留原
      完整 BLAS。统一 pair 坐标 helper 固定 channel、交换相位及相同轨道
      `sqrt(2)` 归一化。emax4/6/8 分别避免
      `1050/14495/109900 KiB` IV dense 输出，standard products 实存从
      `1198/14949/110935 KiB` 降至 `291/883/1993 KiB`。emax8
      setup/IV/addon 从 `0.13488/0.12549/0.86220 s` 降至
      `0.04524/0.04068/0.67228 s`；单/16 线程 wall 为
      `10.95 -> 10.80 s`、`6.95 -> 6.85 s`。至此 emax8 MR addon 只占
      单线程 profiler real 的 `6.3%`。point7 没有可追溯的 canonical
      emax10 输入，因此没有虚构或临时裁剪一个；改用现存 canonical
      NNLOopt `emax=12,e2max=24` 做更强容量门。输入 job `100531` 对
      182 个 J-orbit、156 个 channel、73842036 条 TBME 完成 J64 零误差
      round-trip 和 jref 嵌入，参考收缩/Hermiticity 误差为
      `1.85e-15/0`。固定 `ds=smax=1e-4` 的 64-thread RK4 job `100535`
      在 node8 以 `COMPLETED 0:0` 结束：进程 wall `185 s`、峰值
      `26264640 KiB`，J64/no2bpack 均成功物化且反正规序零体差为 0。
      MR addon 为 `8.69811 s`，只占 profiler real 的 `4.73%`；V build
      为 `1.13381 s`（`0.6%`），standard products 实存峰值只有
      `6494 KiB`，同时避免约 `2.30 GiB` IV dense 输出。因此当前实测
      不支持新增 recoupling plan/cache，MR-P2 以“容量通过、无剩余 MR
      热点值得复杂化”完成。
- [x] **MR-P3：每个优化提交的强制回归。** 每个边界清楚的优化必须先过
      MR-P0/P1 的小空间代数与生产入口测试，再重跑至少一个 emax4/6 RHS
      性能点；要求数值误差保持既有门禁、Hermiticity/anti-Hermiticity 与
      线程复现不退化，并记录优化前后 wall/RSS。当前阶段不以 NCSM Nmax
      收敛或长流谱行为作为优化验收指标。
      本轮两个优化均通过 `MRReference` 随机 J↔m/参考实现、真实相关
      `MRDriver`、`MRSRDriver`，第一项另通过四核完整 oracle 和 sanitizer
      driver；第二项的 emax4/6 输出相对第一轮冻结基线均逐字节相同。
      活跃行存储提交 `b0c832a8` 又通过 Release `MRReference`、`MRDriver`、
      `MRSRDriver`、`MRDenominators` 和完整四体系 `MRCorrelatedDriver`
      （`354.88 s`）；ASan+UBSan 下后三项通过，并用 emax4 嵌入参考直接
      跑过稀疏活跃行短流。sanitizer 输出及 emax4/6/8 Release 输出均与
      冻结 J64 bitwise 相同。本轮只验证 RHS/短流实现，不把有限-s NCSM
      行为列为性能优化门禁。
      成对查询提交 `21760bbf` 又通过相同的 emax4/6/8 bitwise J64 回归、
      emax6 1/16 线程同 SHA、Release `MRReference`、`MRDriver`、
      `MRSRDriver`、`MRDenominators` 和完整四体系 oracle（`355.27 s`）。
      保留在 `/tmp` 的增量 ASan+UBSan 构建已通过 `MRDriver/MRSRDriver`
      与 emax4 稀疏嵌入短流，未发现 `GetTBME_J_twoOps` 索引或三角筛选
      问题。
      V 迹专用提交 `287adea4` 又通过 emax4/6/8 与冻结 J64 bitwise 回归、
      emax8 1/16 线程同 SHA、Release 四项 MR 门禁和完整四体系 oracle
      （`352.75 s`）；增量 ASan+UBSan 的生产 MR/SR driver 及 emax4
      partial-active 短流通过。
      IV 迹专用提交 `460a009d` 又通过相同三档 bitwise J64、emax8 1/16
      线程同 SHA、Release `MRReference/MRDriver/MRSRDriver/MRDenominators`
      及完整四体系 oracle（`361.73 s`）；增量 ASan+UBSan 的生产 MR/SR
      driver 和同时覆盖空/full/partial channel 的 emax4 短流通过。最后的
      emax12 容量门只放大已通过这些回归的相同 C++ 路径，没有再引入算法
      改动；输入/可执行文件/输出哈希、wall/RSS 与 profiler 已写入
      `docs/MR-IMSRG-Jscheme-large-space.json` schema v14，本轮 P3 完成。
- [x] **MR-P4：生产 C++ J-scheme 严格脱耦门禁。** 固定同一个
      `He4, Nrefmax=2, hw=20, emax=2/e2max=4` 相关参考态、NNLOopt 输入、
      White-NCSM/EN 分母和 `Delta e != 0` 掩码，从 `s=0` 先流到 `1000`，
      再由 lossless J64 checkpoint 继续，直到预先冻结的
      `Rgen/Rgen(0)<=1e-6`。point7 jobs `100548/100550` 用完全相同输入和
      executable，分别以 `rtol=atol=1e-10/1e-9` 单线程完成；独立逐通道
      重算得到终点比值 `9.94895896e-7/9.94895893e-7`，两条均严格通过，
      不是 flowfile 的打印舍入。终点 `Rgen` 为约 `7.0e-7`，`eta1` 仅
      `1.26e-11`，没有小分母通道；物理流在累计 `s≈4.0774e5` 过线，随后
      generator 被停止判据置零，ODE 只快速前进到名义输出点。
      两种容差的 lossless J64 最坏 0B/1B/2B 差为
      `5.80e-12/3.19e-12/4.28e-10 MeV`，Nmax2 三态最大变化
      `2.38e-8 keV`。生产 driver 从严格终点 J64 重启后执行 0 次
      commutator，0B/1B/2B round-trip 最坏差低于 `1.60e-14 MeV`；MR
      反正规序零体核对为零。J64/no2bpack 由同一个 NCSM reader 读回，
      59 维三态最大格式差 `9.35e-4 keV`。最后
      `MRReference/MRDriver/MRSRDriver/MRCorrelatedDriver/MRDenominators/
      MRJobGenerator/MRTailDiagnostic` 七项以 `417.82 s` 全通过，其中四核
      m-scheme oracle 用时 `373.60 s`。完整哈希、资源和诊断见
      `docs/MR-IMSRG-Jscheme-large-space.json` schema v15。该门禁只验收
      MR-IMSRG 实现、停止判据、重启和物化，不重新启动 finite-s IM-NCSM
      效果研究。
- [x] **MR-P5：完成审计与整库回归。** 逐条把当前 `CLAUDE.md`、本文件
      P0--P4 和实现计划映射到生产源码与机器证据；确认 production driver
      使用 `MRReference -> Generator::ConstructGenerator_WhiteNCSM ->
      IMSRGSolver::EvaluateCommutator -> MRCommutator::Commutator`，其中公共
      收缩复用原 SR 路径、`lambda2=0` 直接返回原 SR commutator，相关参考
      只增加经过 oracle 验证的 `lambda2` 项，生产 RHS 不展开 m-scheme。
      审计发现 `src/CMakeLists.txt` 的六个重复上游测试仍用相对 build 目录的
      错误脚本路径；改成 `${CMAKE_SOURCE_DIR}/test/...` 后，从当前 checkout
      重新配置并执行完整 CTest，上游 SR/通用测试和八项 MR 测试共 20/20
      通过，总用时 `658.08 s`，四核相关参考 oracle 用时 `505.11 s`；完整
      命令记录在实现计划 Sec. 6.4。
      TODO 中唯一 `[~]` 是用户已明确终止的 emax8 finite-s NCSM 研究，不是
      MR 实现遗留项；当前声明范围内没有未关闭的实现或验收任务。
- [x] **MR-P6：生产 flow 参数继承。** 根目录 `--mr-jscheme` 固定接入原有
      `method=flow`，删除 wrapper 的 method 选择以及
      `ds_0/dsmax/ode_tolerance` 字段和 CLI 参数，不再用 MR 专用默认值覆盖
      `imsrg++`。生成脚本 metadata 显式记录 ODE 参数来自 executable 的
      runtime defaults 且 override 列表为空，schema 因接口变化升为
      `mrimsrg_cpp_jscheme_slurm_v2`；`target_s` 仍作为本段作业边界，
      `eta_criterion` 仍作为物理解耦停止条件保留。

“最终能量接近”不能替代 E--F；只有
`E/f/Gamma -> denominator -> eta -> named RHS contractions -> flow -> vacuum H -> NCSM`
整条链逐层通过，才能称为 C++ J-scheme 生产 MR 实现完成。

## 已完成基准

- [x] Python m-scheme MR 原型的 RDM、正规序、QCombo/显式 Fock-space
      commutator、White-NCSM、真空物化与 NCSM/no2bpack 闭环。
- [x] Python 生产 MR 路径在 `He4/O16, Nrefmax=0` 下对当前 C++
      `imsrg++` 逐正规序、分母、eta、命名 RHS 收缩、Euler/RK4 和
      `s=100` 完整流严格退化；证据见 `docs/MR-IMSRG-验收结果.md`。
