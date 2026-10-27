# MR-IMSRG（IM-NCSM 型）文档

本仓库在 imsrg++ 中实现显式多参考 IMSRG：参考态是 Nmax_ref 截断空间里的 NCSM（以后是 FCIQMC）本征态，
以它的密度做 MR 正规序，在 J-scheme 下用 Magnus 流动 Hamiltonian 与任意算符，把 Nmax_ref 盒子与外部解耦，
导出零、一、两体真空形式的 H(s)、O(s) 给下游 NCSM/FCIQMC 后对角化。

- 截断 MR-IMSRG(2)：对易子保留 λ2；参考态带 λ3 时进入零体项；生成元分子不含 λ。
- 三体力：MR-NO2B（相对相关参考态正规序，保留零、一、两体），full 与 no2b 文件都可用。
- 参考态：总 J = 0、正宇称，只用秩 0 密度；基可以是 HO、HF（imsrg++ 求）或参考态自身的自然轨道。
- 生成元：White-NCSM（Vobig 2020）；求解：`magnus_adaptive`（Bernoulli Magnus ODE + Dormand–Prince 5(4)）。
- 已核对：SR 极限逐位等于 imsrg++；相关参考态对 m-scheme 精确结果；同力文献（¹⁶O NNLOopt、¹²C/²⁰O 3NF400）。

## 流程

```
参考态：BIGSTICK（HO 或 HF 基）/ simpleFCI 原型 / 以后 FCIQMC
        │  wfn → 秩 0 的 γ1、γ2（、λ3）
        ▼
jref 文件：自然轨道占据、基变换 C_total、λ2，可选 λ3、E_ref         guide/参考态输入.md
        ▼
imsrg++：换到 MR 基、MR 正规序（3N：MR-NO2B）、White-NCSM 流、    guide/驱动参数与输出.md
         Ω 变换算符、撤销正规序、导出（HO 基或 MR 基）              guide/三体力与λ3.md、guide/算符演化.md
        ▼
H(s)、O(s)：no2bpack / jcoupled64 / opbin + .mrmeta.json
        ▼
后对角化：BIGSTICK（postdiag_bigstick.py，含标量算符期望值）/ FCIQMC
```

## 快速上手（He4，Nmax_ref = 2，emax2）

1. 参考态（输入在 `test/data/bigstick_He4_Nrefmax2/`，`spectrum.in` 读当前目录的 `emax2.MFDn`）：

   ```bash
   ln -s TwBME_N2LO_opt_bare_hw20_emax2_e2max4.MFDn emax2.MFDn
   bigstick-openmp.x < spectrum.in     # 菜单 n，截断 y / 2，写 he4_nref2.wfn
   bigstick-openmp.x < basis.in        # 菜单 ba，同一空间与截断，写 he4_nref2.bas
   python3 -m prototype.mrimsrg.bigstick_jref --bas he4_nref2.bas --wfn he4_nref2.wfn \
       --interaction emax2.MFDn --nrefmax 2 --hw 20 --output He4.jref
   ```

2. 同一相互作用转成 imsrg++ 读的 J64（含完整 A 依赖真空 H）：

   ```bash
   mrimsrg_minipack_to_j64 --interaction N2LO_opt_emax2.minipack --output He4_bare.jcoupled64 --A 4
   ```

3. 流动 H 并变换算符：

   ```bash
   mkdir omega_scratch
   imsrg++ 2bme=He4_bare.jcoupled64 fmt2=jcoupled64 3bme=none \
     reference=He4 valence_space=He4 A=4 hw=20 emax=2 basis=oscillator \
     method=magnus_adaptive nsteps=1 core_generator=white-ncsm \
     denominator_partitioning=Epstein_Nesbet denominator_delta=0 \
     BetaCM=0 nucleon_mass_correction=false \
     mr_reference_file=He4.jref smax=1e6 eta_criterion=1e-6 \
     scratch=omega_scratch Operators=HCM,Rp2,E2 flowfile=flow.dat intfile=He4 \
     write_H_no2bpack=H.no2bpack
   ```

   日志 `MR initial natural-basis Hamiltonian: E=` 应等于参考态 NCSM 能量（−17.01349 MeV）。这个体系有慢尾
   （到 η 判据约 8×10^5 步、单线程 44 h，见 method/生成元.md），试跑先用 `smax=0.5`。

4. 后对角化与标量算符期望值：

   ```bash
   python3 prototype/mrimsrg/postdiag_bigstick.py --no2bpack H.no2bpack --Z 2 --N 2 --emax 2 --hw 20 \
       --nmax 0 2 4 --operator HCM=H_HCM.no2bpack --operator Rp2=H_Rp2.no2bpack --workdir pd --output pd.json
   ```

HF 基参考态的完整链见 guide/参考态输入.md 2.1 节；集群作业用 `gen_job.py`（guide/驱动参数与输出.md 第 5 节），
文献对照用过的作业脚本见 results/文献对照.md。

## 文档地图

| 分类 | 文档 | 内容 |
|---|---|---|
| 使用 | [参考态输入](guide/参考态输入.md) | jref 格式与约定；BIGSTICK（含 HF 基）、simpleFCI、嵌入三种来源；读入校验 |
| | [驱动参数与输出](guide/驱动参数与输出.md) | MR 模式必设参数与拒绝项、method、Ω 存放、导出文件与元数据、日志、BIGSTICK 后对角化、集群作业 |
| | [算符演化](guide/算符演化.md) | 标量、张量算符经 Ω 的变换、导出与读回 |
| | [三体力与 λ3](guide/三体力与λ3.md) | MR-NO2B、λ3 的两处用法、no2b 文件的限制 |
| | [性能与并行](guide/性能与并行.md) | 并行方式、环境变量、可复现性、内存 |
| | [测试](guide/测试.md) | 每项 ctest 检查什么、怎么跑、需要哪些数据 |
| 方法 | [MR 对易子与正规序](method/MR对易子与正规序.md) | 标量、张量的 λ2 项，MR 正规序，BCH/Magnus 中的 MR 防护 |
| | [生成元](method/生成元.md) | White-NCSM 的定义与慢尾问题 |
| | [公式与代码映射](method/公式与代码映射.md) | MR-IMSRG(2) 收缩、J 耦合约定与代码位置 |
| 结果 | [文献对照](results/文献对照.md) | ¹⁶O NNLOopt（Gysbers）、¹²C/²⁰O 3NF400（Gebrerufael）同力对照 |
| | [验收结果](results/验收结果.md) | 2026-08 至 10 的验收证据（历史） |
| 数据 | `data/literature/`、`data/acceptance/` | 文献取点与验收的 JSON（脚本和测试会读） |
| 记录 | [工作记录](records/工作记录.md) | 2026-10-27 起的进展；待办在根目录 `TODO.md` |
| | [实现计划](records/实现计划.md) | 各阶段的决策与数据（历史，不再维护） |
| | [远程目录归档](records/远程目录归档.md) | 2026-09-08 集群目录整理的路径映射 |

m-scheme Python 原型与输入桥见 `prototype/mrimsrg/README.md`；上游 imsrg++ 的 Doxygen 文档是本目录下的 `*.dox`。

## 代码位置

| 内容 | 位置 |
|---|---|
| 参考态类、读入与校验、MR 正规序 | `src/MRReference.{hh,cc}` |
| 3N 的 MR-NO2B | `src/MRThreeBody.{hh,cc}` |
| MR 对易子（标量与张量 λ2 项、λ3 零体项） | `src/MRCommutator.{hh,cc}` |
| BCH/Magnus 的 MR 分派与防护 | `src/BCH.{hh,cc}`、`src/IMSRGSolver.cc`（`Solve_ode_magnus`） |
| White-NCSM 生成元 | `src/Generator.cc`（`ConstructGenerator_WhiteNCSM`） |
| 驱动（参数检查、正规序、导出、算符） | `src/imsrg++.cc` |
| jref 生成与 λ3 | `prototype/mrimsrg/bigstick_jref.py`、`export_jref.py`、`embed_jref.py`、`lambda3.py` |
| 后对角化、算符读回、no2bpack 组合 | `prototype/mrimsrg/postdiag_bigstick.py`、`no2bpack_combine.py` |
| 集群作业生成 | `gen_job.py`（`--mr-jscheme`、`--mr-prepare-jscheme`、`--mr-ncsm-readback`） |
| 测试 | `test/UnitTestMR*.py` |

## 主要资料

- H. Hergert, *In-Medium SRG for Closed and Open-Shell Nuclei*, arXiv:1607.06882：广义正规序、不可约密度、MR-IMSRG 方程。
- E. Gebrerufael et al., PRL 118, 152503 (2017)，arXiv:1610.05254；Gebrerufael 博士论文（2017）：IM-NCSM 的提出与早期生成元。
- R. Vobig 博士论文（2020）：Magnus、任意秩张量的 MR 对易子、White-NCSM、观测量演化。
- M. Mongelli 博士论文（2022）：观测量的领头阶三体修正、HF 与 NAT 基比较。
- QCombo（github.com/chenlh73/qcombo）：广义 Wick 定理的符号工具，用于核对 MR-IMSRG(2) 收缩。
- 论文 PDF 在本地 `refs/`（不进仓库）。
