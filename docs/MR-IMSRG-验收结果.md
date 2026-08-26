# NNLOopt emax2 MR-IMSRG 验收结果

> 当前状态（2026-08-25）：生产 MR 路径对当前 `imsrg++` 的
> `He4/O16, Nrefmax=0` 单参考退化门禁已经通过正规序、逐分母、逐生成元
> 元素、逐对易子收缩、共同短流和固定 `s=100` 完整流。验收不依赖零体
> 能量偶合；完整 `E/f/Gamma`、`eta`、RHS 和真空 Hamiltonian 均逐元素
> 比较。

本文记录快速 m-scheme MR-IMSRG(2) 原型的可重复验收证据。大型
Hamiltonian、RDM 和日志均保留在被 Git 忽略的本机/集群结果目录，
不提交到仓库。表中能量单位均为 MeV。

## 1. 固定输入与数值定义

- 相互作用：NNLOopt（文件名 `N2LO_opt`），`hw=20`, `emax=2`,
  `e2max=4`，仅 NN，无 Coulomb/显式 3N。
- 文件 SHA-256：
  `76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84`。
- `mrimsrg_prepare`、`mrimsrg_validate` 和 Python runner 均会独立计算
  该摘要并拒绝其他文件；不只信任路径或 metadata 中的文字。
- MR-IMSRG(2) 对易子保留线性 `lambda2`，设 `lambda3=0`；直接积分
  `dH/ds=[eta,H]_(0,1,2B)`。
- 生成元是 Vobig Sec. 6.5.4 定义的 `White-NCSM`：分子舍去所有
  不可约密度项，Epstein--Nesbet 分母舍去文献标记的
  `O(lambda2)` 项；`|Delta|<1e-6 MeV` 时统一替换为正的
  `+1e-6 MeV`，与当前 `src/Generator.cc` 的实际实现一致。
- `Rgen` 是 Vobig 式 (6.5.28)--(6.5.29) 中带 EN 分母、反对称化后
  的实际掩码生成元在完整 m-scheme 系数张量上的范数，是正式固定点判据。
  现有 `IMSRGSolver` 监测 J-coupled `Eta.Norm()`，两者对相同矩阵元采用
  不同的二体简并度权重，数值并不相等；本原型按自己的初值采用相对
  `1e-6` 门槛。完整流必须同时报告两种范数，不能直接比较数值或据此声称
  停止点相同。`Rnum` 是不带分母的
  lambda-free White-NCSM 分子；`Rstrict` 是含线性 `lambda2`
  的 `D-D^dagger` 诊断。对分数占据，分母加权与反对称化不对易，
  因此三者始终分别报告。
- 放松 IM-NCSM 掩码只消去 `Delta e != 0` 的 1B/2B 通道；同 HO
  量子数的参考空间内部块保留，因此最终零体项不必等于后
  NCSM 基态能。

## 2. 实现正确性证据

### 2.1 生产路径对 `imsrg++` 的 SR 退化核对（已通过）

新增 `prototype/mrimsrg/sr_imsrgpp_check.py`，它不实现第二套 SR 公式。
同一个 float64 `jcoupled64` 普通 Hamiltonian 被装入当前构建的
`pyIMSRG.Operator`；C++ 侧直接调用 `Operator::DoNormalOrdering()`、
`Generator::Update()` 和 `Commutator::Commutator()`。所得 J-coupled
算符按已经通过 NCSM 谱验收的相位和 pair normalization 反耦合到完整
m-scheme，再与生产 `prototype/mrimsrg/` 路径逐元素比较。

代数级本机 oracle 是本仓库 `build-src/src/pyIMSRG.so`。运行前工具核对实际编译
所用的 `Generator.cc`、`Commutator.cc`、`IMSRGSolver.cc` 与当前 checkout；
三个 SHA-256 分别为
`3fbc8d005792a171684fd4ab48f559e719ea0b198669b33a72b852a63ac50a49`、
`8bc9c92bdb1fae4fcb709039a3ae698d6fbd130d9a54b127cad4b059f55b0e47`、
`4c2ee3f198542870b800a0ea41158db93ec361e1216b295f214d5ec90d53b436`。

`He4/O16, Nrefmax=0` 均有 `max|lambda2|=0`。He4 的 288 个 1B、
41472 个 2B SR m-scheme 元素和 O16 的 768 个 1B、294912 个 2B
SR 元素中，被生产 `Delta e != 0` 掩码额外删除的元素均为零。

| 核 | 层级 | 0B max-abs | 1B max-abs | 2B max-abs |
|---|---|---:|---:|---:|
| He4 | 初始普通 H | 0 | 0 | `2.66e-15` |
| He4 | 正规序 H | `3.55e-15` | `7.11e-15` | `2.66e-15` |
| He4 | `eta(s=0)` | 0 | `4.16e-17` | `2.78e-17` |
| He4 | RHS(s=0) | `1.78e-15` | `1.78e-15` | `2.66e-15` |
| O16 | 初始普通 H | 0 | 0 | `3.55e-15` |
| O16 | 正规序 H | 0 | `1.42e-14` | `3.55e-15` |
| O16 | `eta(s=0)` | 0 | `2.78e-17` | `2.78e-17` |
| O16 | RHS(s=0) | `1.42e-14` | `2.66e-15` | `1.55e-15` |

EN 分母由新增的只读 pybind 诊断直接调用 C++
`Generator::Get1bDenominator/Get2bDenominator`。He4 比较 20/26 个
1B/2B 分母，最大差为 `1.42e-14/2.84e-14 MeV`；O16 比较 36/116
个，最大差为 `2.13e-14/4.26e-14 MeV`。两核均没有分母触发
`1e-6 MeV` cutoff。静态核查曾发现 Python 对小负分母保留负号、而当前
C++ 统一取正 cutoff；生产 Python 现已改为 C++ 的实际约定并增加
`[-1e-12,0,+1e-12] -> +1e-6 MeV` 回归测试。

生产 `commutator()` 现在直接求和一份可诊断的命名收缩表，不存在另写的
测试公式。八个 SR 项分别调用当前 C++ 的
`comm110/220/111/121/221/122/222_pp_hh/222_phss` 比较；两个额外 MR
`lambda2` 项在闭壳参考下单独验证为零。He4 各项全 rank 最大误差不超过
`3.55e-15 MeV`，O16 不超过 `1.42e-14 MeV`，具体最坏项都是 0B
`comm220ss`。因此总 RHS 的符合不是不同收缩项误差互相抵消。

同一点也确认不能把两边的范数名称当成相同定义：He4 的完整 m-scheme
生成元范数为 `0.4638696070`，当前 C++ `Eta.Norm()` 为 `0.7628643319`；
O16 分别为 `0.6991543919` 和 `1.2913659721`。逐矩阵元仍在双精度噪声内
一致；完整流将以共同目标 `s` 比较矩阵元，并分别记录各自范数。

用相同 RHS 做一个 `ds=1e-4` Euler 步并在步后重新生成 `eta/RHS`，He4
的 `H/eta/RHS` 全 rank 最大差分别为 `7.11e-15/4.16e-17/2.22e-15`
MeV，O16 分别为 `1.42e-14/2.78e-17/1.07e-14 MeV`。这说明退化不只在
`s=0` 单点成立。

进一步让两边运行完全相同的固定步 RK4，`ds=1e-3`，在
`s=0.001,0.002,0.003` 比较完整 Hamiltonian。He4 三点的全 rank 最大差
均为 `7.11e-15 MeV`，终点 `eta/RHS` 最大差为
`4.86e-17/3.55e-15 MeV`；O16 三点的 Hamiltonian 最大差均为
`1.42e-14 MeV`，终点 `eta/RHS` 最大差为
`4.16e-17/1.07e-14 MeV`。以上已通过 P0 的输入、正规序、分母、生成元、
RHS、Euler 单步和共同 RK4 checkpoints 门禁；完整直接流结果如下。

完整流在 point7 的独立干净 checkout 中完成。生产固定-s 作业从 commit
`6974c93b` 启动；正式验收器 commit 为 `4badb01e`，物化 C++ 真空文件的
补充验收器为 `6ca85738`。后两次提交只改诊断/输出，不改 C++ oracle 或
生产流方程。编译与运行环境为 GCC `10.2.1`、Boost `1.81.0`、GSL
`2.7.1`、`openblas/0.3.10-single`、Python `3.12.9`、NumPy `2.2.4`、
SciPy `1.15.2`、SymPy `1.13.3`；`ldd` 没有 `not found`。三份 oracle
核心源码仍逐文件通过上列 SHA-256 核对。

两边都使用直接 flow、White/EN、正的 `1e-6 MeV` cutoff，并关闭提前
`eta` 停止；固定比较点为 `s=100`。Python 使用 DOP853，当前 C++ 使用
Boost odeint Dopri5，因此相同 `rtol=atol` 下预期只在 ODE 截断误差内
相同，而不是逐 bit 相同。所有初点和终点选中分母均未触发 cutoff。

| 核 | `rtol=atol` | Python 0B | C++ 0B | H 最大差 | eta 最大差 | RHS 最大差 | 真空 H 最大差 |
|---|---:|---:|---:|---:|---:|---:|---:|
| He4 | `1e-8` | -20.244646589352 | -20.244646543773 | `2.22e-7` | `3.51e-9` | `3.19e-7` | `2.22e-7` |
| He4 | `1e-9` | -20.244646550834 | -20.244646544175 | `1.11e-8` | `1.77e-10` | `1.60e-8` | `1.11e-8` |
| O16 | `1e-8` | -63.721496035277 | -63.721495961970 | `7.33e-8` | `7.47e-11` | `6.02e-9` | `8.03e-9` |
| O16 | `1e-9` | -63.721495973051 | -63.721495962337 | `2.20e-8` | `3.23e-10` | `3.00e-8` | `2.20e-8` |

表中“最大差”取 0B/1B/2B 的最大绝对值。tight He4 的最坏 H 元素是
`Gamma[0,1,12,13]`，Python/C++ 分别为
`1.1821796358e-6/1.1933277485e-6 MeV`；tight O16 是
`Gamma[12,13,29,30]`，分别为
`2.2322531542e-7/2.0122705856e-7 MeV`。这些接近零的流后非对角元采用
绝对门槛判断。tight 的分母最大差为 He4 `8.39e-9 MeV`、O16
`5.77e-9 MeV`。

同一点分别记录不同范数定义：tight He4 的 Python m-scheme/C++
`Eta.Norm()` 为 `1.31587e-7/1.78308e-7`，tight O16 为
`5.04772e-8/1.02354e-7`。它们反映同一组已经逐元素对上的 eta，但二体
简并度权重不同。Python 从 `1e-8` 收紧到 `1e-9` 后，He4/O16 最终 0B
分别只变化 `3.85e-8/6.22e-8 MeV`，即
`3.85e-5/6.22e-5 keV`，远低于 `1 keV`。

正式生产 Slurm job 为 `100142--100145`，正式 retry 验收 job 为
`100154--100157`。最初的 `100148/100149` 在启动计算前因计算节点没有
`git` 命令而退出；生成器随后改为在登录节点解析并把 commit 字面量写入
脚本，未手改 Slurm、未改变物理输入。代表性命令为：

```bash
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 \
  --smax 100 --checkpoint-s 50 --rtol 1e-9 --atol 1e-9 \
  --max-step 10 --residual-ratio 1e-14 --label sr_s100 \
  --partition c128m512 --result-root /tns/mengziyan/mr-imsrg-sr-results

python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 \
  --sr-check-flow <fixed-s-flow> --jcoupled64 <s0-float64-file> \
  --pyimsrg-dir build-src/src --sr-check-ode-tolerance 1e-9 \
  --label s100_retry --partition c128m512
```

tight 流随后分别由 job `100160/100161` 调用实际 C++
`UndoNormalOrdering()`，物化成验收专用 float64 J-coupled 普通
`E0+t+V`。与生产 Python `to_vacuum()` 的逐元素最大差仍为上表的
`1.11e-8/2.20e-8 MeV`；真空零体常数在两边均只剩约 `1e-13 MeV`
的舍入噪声。两份 float64 文件分别交给相同的 `simpleFCI` NCSM reader，
不是只根据矩阵元误差推断谱：

| 核/空间 | 态 | Python float64 J | C++ float64 J | 差值 |
|---|---:|---:|---:|---:|
| He4, Nmax=8 | 0 | -20.244646550839 | -20.244646544180 | `6.66e-9` |
|  | 1 | -16.319249831823 | -16.319249830124 | `1.70e-9` |
|  | 2 | -16.070458954444 | -16.070458952310 | `2.13e-9` |
| O16, Nmax=2 | 0 | -63.721495973051 | -63.721495962338 | `1.07e-8` |
|  | 1 | -54.826082264030 | -54.826082257894 | `6.14e-9` |
|  | 2 | -54.454755882982 | -54.454755876976 | `6.01e-9` |

所有差值单位为 MeV，且各态 `2J` 顺序也一致：He4 为 `0,4,0`，O16
为 `0,0,4`。生产 tight Hamiltonian 同时导出既有 float32
`no2bpack` 做格式验收；基态相对 float64 J readback 的量化差为 He4
`4.78e-7 MeV = 0.000478 keV`、O16
`2.87e-6 MeV = 0.002872 keV`。因此 double-precision 角动量/真空转换
和 NCSM 谱闭环通过，float32 差异继续只作为既有下游格式精度记录。

- Python 回归测试：`44/44` 通过。覆盖 RDM/cumulant、普通与 MR 正规序
  往返、QCombo/显式 Fock-space 对易子、SR 极限、生成元、掩码、
  Hermiticity/anti-Hermiticity、自然基协变及输出；其中专门的回归测试
  区分了自然基和 HO 基掩码，并确认正式验收使用实际 `eta`
  而非不带分母的分子。
- QCombo 重生的 MR-IMSRG(2) 0B/1B/2B 公式与随机直接求和最大差异
  `1.8e-15`。
- 六个参考态的 `Tr(gamma1)=A` 误差不超过 `5.0e-14`，
  `gamma2 -> (A-1)gamma1` 收缩误差不超过 `7.2e-15`。
- `s=0` 真空 Hamiltonian 经 C++ bridge 交给现有 `shell-model-obs/simpleFCI`
  后，全 `He4, Nmax=8` 3060 维空间得到
  `-20.3388325043`, `2J=0`，精确复现冻结基准。
- 球形 EN 分母专门测试用两个 `j=1/2` 多重态验证：由 `J=0,1` TBME
  `2,6 MeV` 得到的 m-scheme 对角元 `6,4,4,6 MeV` 被还原为
  `(1*2+3*6)/4=5 MeV` monopole，与 `GetTBMEmonopole()` 一致。
- `no2bpack` 导出器会在写文件前做完整 m-to-J-to-m 重构。`s=0` He4
  的最大二体误差为 `2.66e-15 MeV`；采用修正生成元的新 He4 短流
  `s=0.001` 的一体、二体误差均为 `3.55e-15 MeV`。
- 新短流的 He4 `Nmax=8` dense NCSM 基态为 `-20.3396323958 MeV`，
  `no2bpack` 原生 reader 得到 `-20.3396333376 MeV`，相差
  `9.42e-7 MeV = 0.000942 keV = 0.942 eV`；这是既有格式用 float32
  保存 OBME/TBME 的量化舍入误差。先前把该差值写成 `0.942 keV` 是单位
  换算错误，现已更正。
- 同一组 J-coupled OBME/TBME 另以验收专用的 float64 `jcoupled64` 写出、
  独立读回并反耦合到 m-scheme 后，He4 `Nmax=8` 的前三态与 dense 路径在
  能量上的最大差为 `1.8e-14 MeV`，属于 double 数值噪声。这把角动量投影/
  相位/相同粒子对归一化/channel 顺序的误差与生产 `no2bpack` 的 float32
  量化误差分离开；验收后正式下游格式仍使用既有 float32 `no2bpack`，
  不引入新生产格式。
- 同一文件已直接交给 `/home/mengziyan/fciqmc/fciqmc-mpi`：生产
  FCIQMC 可执行文件用 `int_format=no2bpack` 完成 Hamiltonian 初始化和一步
  演化；该仓库的独立 NCSM 可执行文件在 3060 维 `Nmax=8` 空间得到前三态
  `-20.3396333376,-16.6810264996,-16.6581354522 MeV`，对应
  `2J=0,0,4`，与本仓库 no2bpack 验证器逐态一致。

## 3. 旧流撤销与更正

2026-08-24 用新增的严格 m-to-J-to-m 导出验收检查正式流产物时，发现
旧生成元把球形公式中的 `Gamma_ijij` 直接解释为逐 m 的 Slater 对角元。
逐元素除以这种依赖磁量子数的 EN 分母会使 `eta` 不再是旋转标量；
Hermiticity 与二体反对称性仍可通过，因此此前只检查这些对称性没有发现问题。

旧四核最终产物的最大 J-重构误差为：

| 核 | 最大一体误差 (MeV) | 最大二体误差 (MeV) |
|---|---:|---:|
| He4 | 7.602e-4 | 2.230e-3 |
| Be8 | 1.406e-4 | 8.604e-4 |
| C12 | 3.016e-4 | 4.175e-4 |
| O16 | 7.178e-4 | 7.946e-4 |

He4 从同一标量 `s=0` Hamiltonian 仅流到 `s=0.001` 就出现旧实现
`1.247e-4 MeV` 的二体误差，证明这不是长时间 ODE 累积误差。上述旧流、
能量表、容差比较和“收敛加速”结论全部撤销，不作为验收结果；旧目录只保留
作问题复现。正式四核结果必须由 `white_ncsm_spherical_monopole_v1`
重新产生，并逐点通过默认 `1e-9 MeV` J-重构门槛。

## 4. 修正后的四核 `Nrefmax=0` 最终文件

2026-08-25 在 point7 从 `s=0` 重新运行四核。He4/O16 的单参考极限在
`rtol=1e-6, atol=1e-8` 下保持旋转标量到约 `1e-14 MeV`。Be8/C12 的
baseline 虽达到 generator 停止条件，但 m-scheme 独立积分产生与 ODE 容差
同阶的非标量数值漂移，默认 `1e-9 MeV` J-重构门禁分别以
`4.90e-7 MeV` 和 `6.25e-8 MeV` 拒绝导出。两核随后用
`rtol=1e-9, atol=1e-11` 从头重跑；正式文件只取通过门禁的 tight 结果，
没有放宽 exporter 容差或使用被拒绝的 baseline 文件。

| 核 | job | ODE `rtol/atol` | `s_final` | `Rgen/Rgen0` | MR 0B | 最大 J 重构误差 1B/2B |
|---|---:|---:|---:|---:|---:|---:|
| He4 | 100096 | `1e-6/1e-8` | 93.395214 | 6.822e-7 | -20.2446468544 | `7.11e-15 / 5.00e-15` |
| Be8 | 100104 | `1e-9/1e-11` | 331.558295 | 8.778e-7 | -27.6482167218 | `3.15e-10 / 4.96e-10` |
| C12 | 100106 | `1e-9/1e-11` | 233.828911 | 9.246e-7 | -43.3077100404 | `1.13e-11 / 1.03e-11` |
| O16 | 100102 | `1e-6/1e-8` | 81.532065 | 9.148e-7 | -63.7214966519 | `1.07e-14 / 1.86e-14` |

最终文件是现有 NCSM/FCIQMC reader 直接读取的普通真空正规序
`no2bpack`，zero body 用 float64，OBME/TBME 按既有格式用 float32。
Be8/C12 的普通真空 zero body 分别为 `15.1396526724` 和
`17.6301797320 MeV`，reader 日志确认均计入总 Hamiltonian；它们不能被
遗漏，也不能用上表的 MR 0B 代替后 NCSM 对角化结果。

| 核 | NCSM 验证空间 | 维数 | dense 基态 | float32 `no2bpack` 基态 | packing 差值 |
|---|---:|---:|---:|---:|---:|
| He4 | Nmax=8 | 3060 | -20.244646854393 | -20.244646072414 | 0.000782 keV |
| Be8 | Nmax=0 | 51 | -27.942508359012 | -27.942505918504 | 0.002441 keV |
| C12 | Nmax=0 | 51 | -43.530828159008 | -43.530831415028 | 0.003256 keV |
| O16 | Nmax=2 | 1201 | -63.721496651949 | -63.721497668049 | 0.001016 keV |

Be8/C12 另用验收专用 float64 `jcoupled64` 独立写出、读回和反耦合；
Nmax=0 基态与 dense 路径分别相差 `2.8e-14 MeV` 和打印精度内的零，
确认上表剩余差异来自正式格式的 float32 量化。四个 packed reader 的基态
均为 `2J=0`。

集中交付目录在 point7 和本机相同仓库相对位置：

```text
result/mrimsrg-final-nref0/
```

| 文件 | SHA-256 |
|---|---|
| `He4_N2LOopt_hw20_emax2_e2max4_Nrefmax0_MRIMSRG.no2bpack` | `ac69e114b9208036aff8720fd55918d9c6d75ebbf72adda6f35d5c72a6100111` |
| `Be8_N2LOopt_hw20_emax2_e2max4_Nrefmax0_MRIMSRG.no2bpack` | `31a33f6598f58325139ddf7cee892ba176e6cc5290112f6086bd256b9f7a9be9` |
| `C12_N2LOopt_hw20_emax2_e2max4_Nrefmax0_MRIMSRG.no2bpack` | `f0c322e1f4faa348920afcad0319de64d2b079d479f712b22638da20db01da51` |
| `O16_N2LOopt_hw20_emax2_e2max4_Nrefmax0_MRIMSRG.no2bpack` | `031f815bd00cf47086f5f710c3ca57ecd802e02dca473469fe53c17a00da0d13` |

四个文件均为 25392 bytes；它们是按目标质量数读取 bare minipack 后得到的
四份 A-dependent Hamiltonian，不得跨核复用。

## 5. `Nrefmax=2` 相关闭壳参考态

`He4/O16, Nrefmax=2` 均有非零 `lambda2`，且 HO 基中 `gamma1`
非对角。按照 Vobig Sec. 6.5.3，公式和 `Delta e` 掩码均在临时连通块
自然轨道基中计算；自然轨道继承输出槽位的 `e` 标签，最终 Hamiltonian
再变回原 HO 轨道顺序。

| 核 | `||lambda2||F` | HO 基 `max|gamma1_offdiag|` | `s_final` | `Rgen/Rgen0` | `Rnum/Rnum0` | `Rstrict/Rstrict0` | MR 0B | 后 NCSM 基态 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| He4 | 0.598982 | 0.0722519 | 待填 | 待填 | 待填 | 待填 | 待填 | 待填 |
| O16 | 1.082447 | 0.0429950 | 待填 | 待填 | 待填 | 待填 | 待填 | 待填 |

采用旧逐 m 分母的 point7 作业 `100081/100083/100085` 已在确认根因后
取消，不读取其结果。修正实现提交后从 `s=0` 重新提交；不能从旧流继续，
runner 通过 `generator_implementation` metadata 明确拒绝这种 continuation。

## 6. 文件与来源标识

修正后的第一批四核正式流位于 point7；He4/O16 使用 baseline 容差，
Be8/C12 使用 tight 容差：

```text
/tns/mengziyan/mr-imsrg/result/mrimsrg-flow/He4_Nrefmax0_rtol1em06_atol1em08_sphmon_v1/flow
/tns/mengziyan/mr-imsrg/result/mrimsrg-flow/O16_Nrefmax0_rtol1em06_atol1em08_sphmon_v1/flow
/tns/mengziyan/mr-imsrg/result/mrimsrg-flow/Be8_Nrefmax0_rtol1em09_atol1em11_sphmon_v1_final/flow
/tns/mengziyan/mr-imsrg/result/mrimsrg-flow/C12_Nrefmax0_rtol1em09_atol1em11_sphmon_v1_final/flow
```

相关参考态流位于：

```text
/tns/mengziyan/mr-imsrg/result/mrimsrg-flow/<Nucleus>_Nrefmax2_rtol1em06_atol1em08_white_ncsm_published/flow
```

参考波函数 SHA-256 定义为按顺序连接 `orbits.npy`,
`determinants.npy`, `coefficients.npy` 的原始 NPY 字节；RDM SHA-256 定义为
little-endian float64 C-order `gamma1`, `gamma2`, `lambda2` 原始数据字节的连接。

| 参考 | 波函数 SHA-256 | RDM SHA-256 |
|---|---|---|
| He4 Nrefmax0 | `c3cf86c0d2d2c4aa0fc5725c64cbeddaa2eee0fc0da26563c39e2b272aff7c42` | `8eb3c61c2be6979fbad04096153f2a0cd7781e81a4ab5c9f4c1a3b23af403a25` |
| Be8 Nrefmax0 | `20f36fd7da0461e2e7b7f92b98662d0241be69f4482df6a1655aaeae02be469b` | `9e8da10770cf143dcc87382c96aabf64961263fac86a02ac13ec6a468121c538` |
| C12 Nrefmax0 | `3030eb1abd91c6b3766786d9fbe518c2407e5b072fe5d1075793a46aaa2ac70a` | `cceb2d7ecabf6defba55ddef4b5a81e9c87a1012b7708c491edc6209921c297b` |
| O16 Nrefmax0 | `30e38714240e3b351c41392ed6e1cb68d87fe8c4afb6c014203fca17dad22143` | `d555a50d74fd3eb3a3eda59edc63e90f5037c88f2bc6cab1351341b7b529f585` |
| He4 Nrefmax2 | `1ea6d95bb9b8f62cbe9309a45a544a380610475c3a1fe48992e8c89e71832372` | `b160ba51f981af596d960a03ace8a53ba777093c77afe325ac740e6337b90df2` |
| O16 Nrefmax2 | `d00d4a1d87bc4affe6a96bc0cf79f5ea720551890f22cf8c5bb366d73122a003` | `7da1f322109ac2eb629b1e000b8bd829cb0f7f92032cffc7a524a96eadd9bd0f` |

旧 `Nrefmax=0` 与已取消的 `Nrefmax=2` 作业均早于
`white_ncsm_spherical_monopole_v1`，不得作为正式结果或 continuation
起点。正式 `Nrefmax=0` job、能量、门禁和文件摘要已记录在第 4 节。

## 7. C++ J-scheme 固定步退化检查点

生产 C++ driver 与 Python m-scheme oracle 使用相同 float64 Hamiltonian、
RDM、White-NCSM 生成元和 `ds=1e-3` RK4。对 `He4/O16 Nrefmax=2`、
`Be8/C12 Nrefmax=0` 分别积分到 `s=0.001,0.002,0.003`。每一点均从
driver 写出的 HO 真空 `jcoupled64` 重新读入，再做 HO→NAT 与 MR 正规序，
随后逐元素比较 `H`、重新计算的 `eta` 和 `RHS=[eta,H]`。36 个对象比较
全部通过 `2e-10 MeV` 门槛；同一 J 可表示输入展开后，全局最坏 H 为
O16 `s=0.001` 的 `1.279e-13 MeV`，最坏 RHS 为 O16 `s=0.002` 的
`2.665e-14 MeV`。原始 m-scheme 张量投影到严格旋转标量 J 表示的独立
输入门禁最坏为 `6.506e-11 MeV`，不与代数误差混算。

另以 `ds=1e-4` 显式构造 `H_1=H_0+ds RHS(H_0)`，四体系 J/m Euler
一步的全局最坏误差为 `1.776e-15 MeV`。`IMSRGSolver` 中历史名称
`flow_euler` 实际使用 Boost `runge_kutta4`，故本门禁没有借用该名称。

完整数值见 [`MR-IMSRG-Jscheme-checkpoints.json`](MR-IMSRG-Jscheme-checkpoints.json)；
自动门禁为 CTest `MRCorrelatedDriver`，本次完整运行耗时约 250 秒。

同一门禁还把四体系各 10 个命名 contraction 的 0/1/2B 有效输出写入
[`MR-IMSRG-Jscheme-contractions.json`](MR-IMSRG-Jscheme-contractions.json)，
包含 max-abs、绝对/相对 Frobenius、最坏 m-scheme 索引及该处两边数值。
全局 max-abs 为 `1.954e-14 MeV`，全局相对 Frobenius 为
`1.380e-15`。

White-NCSM 的 EN 分母另由 CTest `MRDenominators` 逐条核对生产
`Generator::Get1b/2bDenominatorWhiteNCSM` 与 Python 论文公式。每核覆盖
4 个 1B 和 704 个 2B 正/反向 `Delta e != 0` 通道，四核共 2832 条，
全局 max-abs `1.421e-14 MeV`。紧凑机器表及列定义见
[`MR-IMSRG-Jscheme-denominators.json`](MR-IMSRG-Jscheme-denominators.json)。

## 8. C++ J-scheme 完整流、SR 退化与下游谱

point7 上使用统一 `rtol=atol=1e-9` 和固定终点 `s=100`，对
`Be8/C12 Nrefmax=0`、`He4/O16 Nrefmax=2` 比较 Python DOP853 与 C++
Boost Dopri5 的终点 `H/eta/RHS`，随后反正规序并变回 HO 基比较普通真空
Hamiltonian。四核全局 max-abs 依次为
`4.40e-9/6.88e-8/1.87e-8/1.04e-8 MeV`，均通过 `1e-5 MeV` 的完整流
门限。float64 J64 由独立 `shell-model-obs` NCSM reader 读回后，再转成
下游 float32 `no2bpack` 并重复对角化；打包导致的基态差最大为 Be8 的
`1.225e-6 MeV = 0.001225 keV`。逐对象最坏秩、索引、相对 Frobenius、
三条最低能级、Slurm job 和环境见
[`MR-IMSRG-Jscheme-full-flow.json`](MR-IMSRG-Jscheme-full-flow.json)。

`He4/O16 Nrefmax=0` 的原生 SR 与同一 MR driver 零 cumulant 入口另作严格
退化。输入真空 Hamiltonian 最坏 `3.55e-15 MeV`，EN 分母最坏
`4.26e-14 MeV`，十个命名 RHS 项最坏 `1.15e-14 MeV`，共同 Euler/RK4
短流也保持在 `1.42e-14 MeV` 内。`Delta e` 掩码没有删去任何原生 SR
非对角元素；两个 MR `lambda2` addon 精确为零。固定 `s=100` 完整流最坏
为 O16 RHS 的 `3.00e-8 MeV`。完整命名项、相对 Frobenius 和最坏元素见
[`MR-IMSRG-SR-degeneration.json`](MR-IMSRG-SR-degeneration.json)。

真实 `He4 Nrefmax=2` 单线程 RHS 基准中，C++ J-scheme 三次中位数为
`6.88 ms`，Python dense m-scheme 为 `1.766 s`，加速 `256.8` 倍；两输入
算符加参考密度的主数值数组由 `61.48 MB` 降到 `41.9 kB`。解析 channel
尺寸而不分配算符的存储表延伸到 `emax=14`，证明生产主存储按 J-channel
而非完整 magnetic-substate 四指标张量增长。复现命令和机器数据见
[`MR-IMSRG-Jscheme-performance.json`](MR-IMSRG-Jscheme-performance.json)。
同一解析还发现 ordered particle-hole 收缩不应同时保留全部 J/parity/Tz
块；改成逐块构造和释放后，该项 `emax=14` scratch 估算由 `48.46 GB`
降到 `2.87 GB`，公式与数值结果不变。

十倍 ODE 容差门禁随后在 point7 完成：四核固定在同一 `s=100`，把
`rtol=atol=1e-9` 收紧到 `1e-10`，再用 C++ 物化的 float64 J64 交给独立
NCSM reader。每核比较三条最低能级，12 个差值中的最大值为 Be8 基态的
`1.683e-11 MeV = 1.683e-8 keV`，远低于 `1 keV`；完整 job、谱和误差见
[`MR-IMSRG-Jscheme-full-flow.json`](MR-IMSRG-Jscheme-full-flow.json)。

较大空间的第一个真实运行门禁也已接通：`He4 Nrefmax=2` 的固定参考态
嵌入 NNLOopt `hw=20, emax=4, e2max=8`，新轨道零占据/零 cumulant、NAT
单位块。生产 J-scheme driver 的一个 `ds=1e-4` RK4 步单线程用时
`1.37 s`、峰值 RSS `54,360 KiB`，1/64 线程输出最坏差
`4.44e-16 MeV`；运行未构造稠密 m-scheme Hamiltonian。V/VI 两项是实测
主要瓶颈，详见
[`MR-IMSRG-Jscheme-large-space.json`](MR-IMSRG-Jscheme-large-space.json)。

为了让这个判断可自动复现，`prototype_downstream_stability_v1` 在已测流点中
选择满足所有条件的最大 `s`：全空间基态漂移不超过 `100 keV`，所有已测
Nmax 的同终点截断误差都不变差，packing 与十倍 ODE 容差误差各小于
`1 keV`，并要求生成元 norm 确实下降。100 keV 是本快速原型的工程预算，
不是文献阈值；文献只支持使用后 NCSM 平台监测诱导高体项失控这一方法。
`prototype/mrimsrg/select_downstream_flow_window.py` 从机器记录重算该决定，
CTest `MRDownstreamWindow` 同时验证当前选中 `s=0.02`、`s=0.1` 因漂移
被拒绝，以及严格脱耦门禁仍为失败。

随后先对 VI 做严格等价的有限求和重排：预计算其 `(J2,w)` 迹后，emax4
VI 用时由 `0.492 s` 降至 `0.022 s`，整个短流由 `1.37 s` 降至 `0.88 s`；
优化前后 J64 逐元素无差，emax2 Python oracle 和完整相关参考 CTest 均
保持通过。V 的 ordered Pandya 块现为下一主瓶颈，本报告未把这一阶段写成
性能优化完成。

V 随后复用了三个张量共同的 Pandya recoupling、dense Six-J cache 和已有
双算符 TBME accessor，没有改变指标、相位、有限求和次序或 BLAS
contraction。emax4 V 由 `0.616 s` 降至 `0.332 s`，整个短流由 `0.88 s`
降至 `0.602 s`，相对最初版本快 `2.28` 倍。单线程输出与优化前 bitwise
相同，1/64 线程最坏差仍为 `4.44e-16 MeV`；emax2 Python oracle 的最坏差
仍是 `3.11e-12 MeV`，完整相关参考 CTest 在 `351.93 s` 后通过。曾试验的
单块 OpenMP 在 emax4 没有加速，已撤回；本阶段仍需 point7 完整流/NCSM
回归，不能仅凭短流性能宣称全部优化验收完成。

完整回归随后在 point7 由专用生成器提交 jobs
`100377/100379/100381/100383`。四核均固定原 `s=100, rtol=atol=1e-10`
Python 轨迹，以优化后 C++ 重做 full-flow check；C++/Python 的全局最坏差
分别为 `2.212e-8/3.898e-8/8.870e-10/1.071e-9 MeV`。新旧优化版真空
Hamiltonian 的 0B/1B/2B 全局最坏差为 `7.44e-13/3.11e-13/6.04e-14
MeV`。新 J64 下载后由 `simpleFCI` 实际重算 He4/Be8/C12/O16 的三条最低
谱，按 `Nmax=8/0/0/2` 得到的 12 条能级相对优化前最坏变化
`2.13e-13 MeV`。至此 V/VI 性能重排通过完整流和 NCSM 门禁。

较大空间门禁随后扩到真实 emax6。正式输入由冻结 emax14 母 minipack
（SHA `118c4c27...b9d7d`）流式抽取，canonical child SHA 为
`199d5a7a...b913b`；A=4/16 的 `Hamiltonian::truncate(6)` 对照均严格零差，
旧 candidate 因 `3.8147e-6 MeV` 差而未使用。He4 Nrefmax=2 嵌入参考在
56 个 J-orbit、467032 条 TBME 上完成一个 `ds=1e-4` RK4 步：单线程
`11.70 s`、峰值 RSS `281132 KiB`，流后零体 `-17.016801664428 MeV`，
导出零体一致性为零。1/64 线程 0B/1B 相同，2B 最坏 `8.88e-16 MeV`。
该结果验收的是 emax6 RHS/短流和存储路径，不是完整收敛流。

V 的第三轮优化只利用 `lambda2` 的精确轨道支撑域，把
`X Lambda Y` 恒等改写为 `X(:,A) Lambda(A,A) Y(A,:)`；不使用数值阈值，
全空间活跃参考自动恢复原全矩阵实现。emax4/emax6 短流输出相对上一版均
逐字节相同，墙钟分别降到 `0.408/6.81 s`，emax6 V 从 `7.751 s` 降到
`2.851 s`。SR、随机全活跃参考、emax2 Python oracle 和 1/64 线程门均
通过。point7 generator jobs `100385/100387/100389/100391` 又以
`s=100, rtol=atol=1e-10` 完成 He4/Be8/C12/O16 全流回归；新旧真空 J64
全局 max-abs `8.384e-13 MeV`。下载后由现有 `simpleFCI` 按
`Nmax=8/0/0/2` 重算 12 条低能谱，最大变化 `1.741e-13 MeV`。因此本轮
性能变换已通过完整流、真空物化和后 NCSM 谱三层验收。

随后 commit `dea9125b` 继续避免构造不会进入乘法的 X/Y Pandya 元素，
partial-support 分支只形成活跃列和活跃行，全活跃参考仍保留原完整块。
emax4/emax6 短流相对 `cef4fce6` 均逐字节相同，墙钟降为
`0.312/4.15 s`；emax6 V/build 为 `0.170/0.111 s`。emax2 Python oracle
仍为 `3.106e-12 MeV`，随机全活跃参考和四核 checkpoint CTest 用时
`352.52 s` 并通过。

该版本随后由 point7 专用生成器提交 jobs
`100393/100395/100397/100399`，固定复用四核原
`s=100, rtol=atol=1e-10` Python 轨迹。He4/Be8/C12/O16 的 C++/Python
全流最坏差分别为 `2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`。
新 J64 相对已验收 `cef4fce6` 版本的 0B/1B/2B 全局最坏差为
`8.527e-13/2.576e-13/8.527e-14 MeV`。下载后由现有 `simpleFCI`
按 `Nmax=8/0/0/2` 实际重算 12 条低能谱，最大变化
`1.990e-13 MeV`。所以 `dea9125b` 已通过完整流、真空物化和后 NCSM 谱
三层验收；后续性能工作应转向现有 White-NCSM generator，而不是继续改写
已经不再主导耗时的 MR V contraction。

White-NCSM generator 随后按 Vobig Eqs. (6.5.31--34) 做了纯缓存优化：
每个 RHS 只跨 J 计算一次完整有序轨道对的 `Gamma_ijij` monopole 表，原
分母表达式、算术次序、`Delta e` 掩码、占据权重和 cutoff 均不变。
commit `6e02676c` 的 emax4/emax6 短流相对 `dea9125b` 逐字节相同；
generator update 由 `0.1112/2.191 s` 降至 `0.00324/0.0426 s`，总时由
`0.312/4.15 s` 降至 `0.222/1.98 s`。分母逐项、SR/MR driver、emax2
Python oracle (`3.106e-12 MeV`) 和随机全活跃相关参考 CTest
(`354.53 s`) 全部通过。

point7 generator jobs `100401/100403/100405/100407` 随后重复四核
`s=100, rtol=atol=1e-10` 全流。He4/Be8/C12/O16 的 C++/Python 最坏差
为 `2.212e-8/3.898e-8/8.871e-10/1.071e-9 MeV`；相对 skinny 版 J64
的 0B/1B/2B 全局最坏差为 `9.308e-13/2.895e-13/6.484e-14 MeV`。
`simpleFCI` 按 `Nmax=8/0/0/2` 重算 12 条低能谱，最大变化
`1.137e-13 MeV`。因此 monopole 缓存通过完整流、真空物化和后 NCSM 谱
三层验收。

## 9. emax4 长流停止判据与后 NCSM 流窗

`He4 Nrefmax=2` 嵌入 emax4 后，point7 jobs `100413/100423` 分别完成
累计 `s=100/1000`，但 `Rgen/Rgen0` 只降到 `1.0265e-2/1.3102e-3`。
逐通道诊断排除了小分母、cutoff 和 ODE 容差问题；慢尾来自 White-NCSM
方向分子中的小自然占据权重。后 NCSM 又给出更强的物理否决：Nmax8 基态
从裸值 `-25.153026 MeV` 漂到 `-30.588114/-42.620511 MeV`。所以这两个
终点虽然程序正常完成、物化和校验和完整，却不是可接受的下游 Hamiltonian。

为定位有限流窗，使用新增的受校验 emax4 `simpleFCI` basis 重新对角化
`s=0,0.02,0.1` 的 `Nmax=2,4,6,8,16`。对四粒子 emax4，`Nmax=16`
已经是全空间；三点全空间基态依次为
`-25.291223278/-25.336699633/-25.478750954 MeV`。`s=0.02` 相对自身
全空间的 Nmax2/4/6/8 截断误差为
`7865.585/2782.326/806.046/128.328 keV`，相对裸 Hamiltonian 分别改善
`412.148/147.117/39.306/9.869 keV`，同时引入 `-45.476 keV` 全空间
IMSRG(2) 漂移。`s=0.1` 虽把 Nmax8 截断误差继续降到 `97.750 keV`，
全空间漂移已增至 `-187.528 keV`。

`s=0.02` 的前三条 Nmax8 谱由 lossless J64 和下游 float32 no2bpack
独立读取，最大 packing 差 `0.000726 keV`。同一点以 `rtol=atol=1e-9`
和 `1e-10` 重跑，J64 的 0B/1B/2B 最坏变化为
`1.243e-14/2.842e-14/3.675e-14 MeV`，三条 NCSM 能级最坏变化约
`3e-11 keV`，通过 `<1 keV` ODE 门槛。

验收结论必须分开表述：`s=0.02` 通过格式、积分和“有限短流改善 Nmax
收敛”门禁，可以作为下一阶段的加速器候选；但其
`Rgen/Rgen0=0.9659`，没有通过严格 `1e-6` 脱耦门槛，也不能称为收敛的
MR-IMSRG Hamiltonian。Magnus 可能改善积分和恢复机制，却不会自动消除
这里已经实测到的 IMSRG(2) 长流截断漂移。完整机器数据见
[`MR-IMSRG-Jscheme-large-space.json`](MR-IMSRG-Jscheme-large-space.json)。

## 10. emax6 finite-s 最大可算 Nmax pilot

emax6 使用从冻结 emax14 母文件严格抽取的 canonical NNLOopt
`emax=6,e2max=12` interaction（SHA-256 `199d5a7a...b913b`），参考态仍是
经过验收的 He4 `Nrefmax=2` emax2 NCSM 参考的严格 embedding。新增轨道
占据与 `lambda2` 为零、NAT 变换为单位块；本节不把它称为 emax6 重新求解
的相关参考态。

point7 job `100432` 完成 `s=0->0.02` direct flow，墙钟 `10.48084 s`、
峰值 RSS `886.555 MB`。实际生成元 norm 从 `1.1971838402` 降到
`1.1504815787`，比值 `0.9609899`，所以严格 `1e-6` 脱耦仍失败。job
`100434` 把 ODE 容差由 `1e-9` 收紧至 `1e-10`；J64 的 0B/1B/2B 最坏差
为 `1.18e-13/1.13e-13/4.16e-13 MeV`，Nmax8 三态最坏变化约
`1.2e-10 keV`。生产 no2bpack 与由 J64 独立重写的 no2bpack 虽字节不同，
Nmax8 三态谱在打印的 double 精度内一致，ODE 与格式门都通过。

裸与流后 Hamiltonian 的相同 NCSM 序列得到：

| Nmax | 维数 | 裸基态 (MeV) | `s=0.02` 基态 (MeV) | 流减裸 (keV) |
|---:|---:|---:|---:|---:|
| 2 | 59 | -17.0135028813 | -17.6498957026 | -636.393 |
| 4 | 952 | -22.3617937335 | -22.7103577949 | -348.564 |
| 6 | 7916 | -24.7096267300 | -24.9028182763 | -193.192 |
| 8 | 44350 | -26.1720022056 | -26.2747194774 | -102.717 |
| 10 | 183866 | -26.7907500860 | -26.8606815224 | -69.931 |

相邻 Nmax 的裸 gap 为
`5.348291,2.347833,1.462375,0.618748 MeV`，流后为
`5.060462,2.192460,1.371901,0.585962 MeV`，四段均改善。Nmax10 裸/流
jobs `100441/100442` 均为 Slurm `COMPLETED`, exit `0:0`，墙钟
`174/175 s`、峰值 RSS `448636/449508 KiB`。因此预先写明的两个代理条件
——最大 Nmax 处流/裸差不超过 100 keV、相邻 gap 不变差——均通过。

这个结论严格限于“当前最大可算 Nmax pilot”。四粒子 emax6 的全空间是
Nmax24，本轮没有计算，因此没有 emax6 全空间漂移，也没有宣称完整
`prototype_downstream_stability_v1` 通过。下一阶段可以评估 emax8 的
实际资源和最大可行 NCSM 空间，但不能把本 pilot 外推成 emax8 物理验收。

## 10.1 simple-ncsm 与 BIGSTICK 同 Hamiltonian 严格核对

此前 emax8/Nmax2 sanitizer 运行在 90 个 J-orbit 上触发
`shift exponent 65 is too large for 64-bit type`。原因不是 NCSM 物理空间，
而是 `CIspace` 用一个 `uint64_t` 打包任意长度 occupation partition。隔离的
`/home/mengziyan/simple-ncsm` 已改用 `std::vector<int>` 作为 collision-free
key。用原 SHA-256 为 `e1653533...6946` 的 emax8 bare no2bpack、660 个
m-orbit、He4 `Nmax=2` 复跑，维数仍为 59，三态为
`[-17.013502881259512,-13.761246467855941,-13.761246401684351] MeV`，
进程正常 `exit 0`，未再出现非法移位。

为排除“两个程序实际用了不同相互作用”，严格核对另冻结一份 He4
Hamiltonian。源文件是固定 N2LOopt `A=4,hw=20,emax=2,e2max=4,
BetaCM=0` minipack，SHA-256 为
`76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84`。
identity 导出的普通真空 0B+1B+2B no2bpack SHA-256 为
`5cadb860e75c58662832bc7cceb78b3d8bb778a268cc262a814748bf6c8e01de`；
simple-ncsm 和 BIGSTICK 都直接读取这一份文件，没有分别重写相互作用。

BIGSTICK 在正宇称、`2M=0,Nmax=2` 的 59 维空间执行 exact full
diagonalization 并保留全部 59 态；`simple-ncsm` 对同一空间做 double
precision full diagonalization。新增检查器还把 BIGSTICK `.bas/.wfn`
的每个行列式和振幅映射回 simple-ncsm 顺序，结果为：

| 检查量 | 最坏值 |
|---|---:|
| determinant set / dimension | 完全相同，59 / 59 |
| 完整谱差 | `4.164634e-6 MeV` |
| BIGSTICK 振幅 Rayleigh 商与 simple exact 谱差 | `3.339867e-7 MeV` |
| BIGSTICK 振幅在 simple Hamiltonian 下残差 | `4.136186e-6 MeV` |
| 波函数 norm 误差 | `2.491691e-8` |

门禁为 `1e-5 MeV`，覆盖 BIGSTICK version-1 `.wfn` 固有的 float32
能量与振幅；两边共享的 no2bpack 本身没有二次量化差。检查器位于
`simple-ncsm` commit `fcc6070`。因此 simple-ncsm 的占据分区故障与
small-space NCSM/BIGSTICK 一致性问题均已关闭。该结论不改变大空间规则：
生产谱仍以 BIGSTICK 为正式求解器，simple-ncsm 作为独立交叉检查。

## 11. emax8 canonical input 与资源门禁

emax8 canonical NNLOopt `hw=20,e2max=16` 输入从冻结 emax14 母 minipack
直接抽取 3526624 条原始 records，文件 SHA-256 为
`3fd1a0038e8c6f00cd93ec26113f02fcb19dd51268a95238da67b17e12fd4bda`。
独立 `Hamiltonian::truncate(8)` 比较在 A=4 与 A=16 下均逐 channel 严格
零差。对应 bare J64 SHA-256 为 `cd14cd4b...eb9cce`，嵌入参考为
`9af7a6e0...a68b3`，后者保持原 emax2 NCSM 参考并让新增轨道占据/cumulant
严格为零。

point7 job `100444` 对 90 个 J-orbit、3526624 条 TBME 做固定
`ds=smax=1e-4` RK4 一步。4 次完整 J-scheme RHS 用时 `17.28220 s`，峰值
RSS `2070.297 MB`；流后 MR 零体项 `-17.017091905266 MeV`，真空物化零体
核对为零，并产生 J64/no2bpack。若构造 660 个 m-orbit 的稠密两体 double
张量则需约 `1.518 TB`；旧 lossless J64 validator 因而 `bad_alloc`，这条
dense-m 诊断路径已明确排除，生产流本身没有构造该张量。

native no2bpack NCSM Nmax8 jobs `100447/100448` 均为 Slurm
`COMPLETED`, exit `0:0`，维数 44838。bare 与一步流后三态为：

| 输入 | E0 (MeV) | E1 (MeV), 2J | E2 (MeV), 2J | wall | peak RSS |
|---|---:|---:|---:|---:|---:|
| bare | -26.217052992370 | -24.709627375508, 4 | -24.709627301702, 0 | 81 s | 2385856 KiB |
| s=0.0001 | -26.217786322986 | -24.709987825043, 4 | -24.709967967694, 0 | 75 s | 2386220 KiB |

一步流相对 bare 的基态变化为 `-0.733331 keV`。这一节只验收 canonical
输入、pure J-scheme 资源、真空物化和 native reader 闭环；`s=1e-4` 不构成
有限流的物理验收。资源门通过后才进入下面的 `s=0.02` 双容差流；即使该流
通过 ODE 门，也仍须用相同 bare/flow Nmax 序列做到最大可行空间。若没有
做到 emax8 四粒子全空间 Nmax32，只能报告最大可行 Nmax proxy，不能宣称
完整 downstream stability 门禁通过。

同一 emax8 输入随后完成 `s=0→0.02` 双容差 direct flow。point7 jobs
`100451/100452` 分别使用 `rtol=atol=1e-9/1e-10`，25/31 次 scalar
commutator 的墙钟为 `67.67063/85.03911 s`，峰值 RSS
`5194.719/5227.641 MB`。两者终点 MR 零体项均为
`-17.70238931537 MeV`，真空导出零体核对为零；生成元 norm 比值为
`0.959464`，所以该点不是严格脱耦结果。

lossless J64 逐元素比较的 0B/1B/2B 最坏差为
`4.80e-14/6.53e-14/1.90e-13 MeV`。native no2bpack Nmax8 正式读回为：

| ODE 容差 | job | 三态 (MeV) | Slurm |
|---:|---:|---|---|
| 1e-9 | 100458 | -26.356103569741, -24.779063666931, -24.775134512563 | COMPLETED 0:0 |
| 1e-10 | 100456 | -26.356103569741, -24.779063666931, -24.775134512563 | COMPLETED 0:0 |

三态最大差 `4.97e-14 MeV = 4.97e-11 keV`，ODE 门通过。初次松容差作业
`100454` 虽打印相同三态，但子进程最终返回 1、Slurm 标为 FAILED；它不计入
正式成功，重跑 `100458` 将 `max_iter` 提到 500 后以 exit `0:0` 完成。
下一门禁仍是 bare/flow 匹配 Nmax 序列和最大可行 Nmax proxy；目前不能据
Nmax8 单点宣称 emax8 downstream physics 通过。
