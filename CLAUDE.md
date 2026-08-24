# MR-IMSRG 快速原型仓库说明

## 最高优先级强制规则

1. **及时 commit，不能丢代码。** 每完成一个边界清楚、已经验证的阶段就立即提交；提交前检查 `git status` 和 `git diff`，只提交本任务相关文件，不夹带工作树中原有的其他改动。不得让已完成的实现长期只留在未提交工作树中。
2. **遇到问题先自行调研并作出决定。** 在用户已经授权的任务范围内，主动上网查询官方文档、开源代码、论文及其补充材料；论文必须尽量下载并阅读原文 PDF，不能只看网页摘要。结合代码和证据选择方案并继续推进，不因一般技术选择或实现偏好停下来询问用户。只有缺少必要权限、任务范围发生实质变化或操作具有不可恢复风险时才请求用户决定。
3. **尽量复用已有代码，不重造轮子。** 开始实现前先查找并理解仓库中已有的程序、接口、脚本和测试路径；已有能力可以直接完成任务时必须复用。只有现有实现确实不能满足任务需求时才新增代码，并把新增范围控制在必要的最小边界内，避免复制默认程序已经负责的初始化、运行或验证流程。
4. **仔细核对 `refs/` 中的代码与公式，正确性优先。** 实现 MR 正规序、RDM/cumulant、MR-IMSRG(2) 对易子、生成元或真空正规序转换之前，必须阅读 `refs/` 中对应的论文 PDF、公式笔记和 QCombo 代码/notebook，核清指标顺序、反对称化、组合系数、相位及采用的密度截断。不能凭记忆或只依赖 PDF 文本提取/OCR 抄写长公式；关键公式必须使用 QCombo、显式小 Fock-space 矩阵或另一份独立实现做符号或数值交叉验证，未通过单参考极限、Hermiticity/anti-Hermiticity、RDM 恒等式和随机小模型测试的实现不得用于四核结果。

## 当前唯一目标

本分支不以实现生产级、高效或通用的 MR-IMSRG 为目标。当前目标是尽快得到一个物理定义清楚、可验证、可物化为零/一/二体 Hamiltonian 的 MR-IMSRG(2) 结果，然后交给 NCSM 做确定性验证，也可把同一输出交给 FCIQMC 做后续测试。

执行计划以 `docs/MR-IMSRG-快速结果计划.md` 为准。不要让性能优化、J-scheme 化、显式 3N、一般张量算符或大模型空间重构阻塞第一批结果。

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

## 第一版算法边界

- 独立、低效、便于检查的 m-scheme MR-IMSRG(2) 原型；不要先改写现有高效 J-coupled commutator。
- 参考态由 `Nrefmax` 截断 NCSM 波函数产生。
- 至少输入 `gamma1`、`gamma2`，构造并真正保留 `lambda2`。
- 第一版设 `lambda3=0`；这一近似必须写入输出元数据。
- Hamiltonian 和生成元始终截断到参考态正规序 0B/1B/2B。
- 直接积分 `dH/ds = [eta,H]_(0,1,2B)`；第一版不要求 Magnus。
- 使用 IM-NCSM 的放松解耦：只处理 `Delta e != 0` 的 1p1h 与 2p2h 通道，保留同 HO 量子数参考空间内部耦合。
- 首选一个有文献公式和 QCombo 输出可核对的生成元；第一批只实现一种。不要同时开发多种生成元。
- 第一版不做自然轨道、显式 3N、`lambda3`、奇核、非标量密度、观测算符演化、MPI 或大规模优化。

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

1. 冻结核力、轨道顺序、相位和普通 Hamiltonian 的现有 FCI/NCSM 基准。
2. 打通 NCSM 参考态到 `gamma1/gamma2/lambda2` 的接口和恒等式测试。
3. 实现 MR 正规序及其真空表示往返。
4. 实现并逐项验证 MR-IMSRG(2) commutator。
5. 实现单一生成元、`Delta e != 0` 掩码和直接流积分。
6. 输出普通 0B/1B/2B Hamiltonian，先由 NCSM 读回验证。
7. 按 `He4 -> Be8 -> C12 -> O16` 产生第一批结果。
8. 管线通过后做 `He4/O16, Nrefmax=2`，最后才考虑 FCIQMC 性能测试或代码优化。

## 最低验收门槛

- `Tr(gamma1)=A`，RDM 的 Hermiticity、反对称性和收缩关系通过。
- 单 Slater 参考态给出 `lambda2=0`，并复现现有 SR-IMSRG(2) 极限。
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
- 不要为这个快速原型修改 `gen_job.py`、生产 Slurm 约定或现有 SR/VS-IMSRG 行为，除非任务明确要求。
- 不要提交生成的 Hamiltonian、波函数、RDM、大日志或结果目录；提交小型测试 fixture 时必须说明来源与校验和。
- 工作区可能已有用户修改。只提交当前任务明确生成或修改的文件，不顺手清理、移动或覆盖其他改动。
