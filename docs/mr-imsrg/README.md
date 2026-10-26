# MR-IMSRG（IM-NCSM 型）使用说明

本目录是本仓库显式多参考 IMSRG 的使用说明，按模块写“为什么这样实现、怎么用”。实现计划、阶段记录和验收证据另见：

- [实现计划](../MR-IMSRG-C++-J-scheme实现计划.md)：各阶段的决策与数据（历史记录）；
- [公式与代码映射](../MR-IMSRG-J-scheme公式与代码映射.md)：MR-IMSRG(2) 收缩、J 耦合约定与代码位置；
- [验收结果](../MR-IMSRG-验收结果.md)：SR 退化、完整流、下游谱的证据；
- [原型说明](../../prototype/mrimsrg/README.md)：m-scheme Python 原型与输入桥。

## 做的是什么

IM-NCSM 式的 MR-IMSRG(2)：参考态是 Nmax_ref 截断空间里的 NCSM 本征态（J = 0），以它的一体、两体密度做多参考正规序，
在 J-scheme 下流动 Hamiltonian（及任意标量、张量算符），把 Nmax_ref 盒子与外部脱耦，导出真空正规序、HO 基的 H(s) 给下游 NCSM/FCIQMC。

- 截断：MR-IMSRG(2)，λ3 = 0（λ2 在对易子里保留，生成元分子里不含）；
- 参考态：总 J = 0、正宇称，只用秩 0（标量）密度；
- 生成元：White-NCSM（Δe ≠ 0 掩码、Epstein–Nesbet 分母）；
- 求解：生产用 `magnus_adaptive`（完整 Bernoulli Magnus ODE + Dormand–Prince 5(4)），算符经同一个 Ω 变换。

## 流程

```
参考态 NCSM（simpleFCI 原型 / BIGSTICK；以后 FCIQMC）
        │  密度：γ1、γ2 的秩 0 部分
        ▼
jref 文件（自然轨道占据、HO→NAT 变换、λ2）          ← 参考态输入.md
        │  可选：嵌入更大 emax
        ▼
imsrg++ 驱动：HO→NAT、MR 正规序、White-NCSM 流、       ← 驱动参数与输出.md
             Magnus Ω 变换算符、撤销正规序、NAT→HO          生成元.md / 算符演化.md
        │                                                    MR对易子与正规序.md
        ▼
H(s)、O(s)：jcoupled64 / no2bpack / opbin + .mrmeta.json
        │
        ▼
下游：NCSM 读回（mrimsrg_validate）或 FCIQMC
```

## 快速上手（He4，Nmax_ref = 2，emax2）

1. 参考态（BIGSTICK，单进程即可；输入在 `test/data/bigstick_He4_Nrefmax2/`，`spectrum.in` 读当前目录的 `emax2.MFDn`）：

   ```bash
   ln -s TwBME_N2LO_opt_bare_hw20_emax2_e2max4.MFDn emax2.MFDn
   bigstick-openmp.x < spectrum.in     # 菜单 n，截断 y / 2，写 he4_nref2.wfn
   bigstick-openmp.x < basis.in        # 菜单 ba，同一空间与截断，写 he4_nref2.bas
   python3 -m prototype.mrimsrg.bigstick_jref --bas he4_nref2.bas --wfn he4_nref2.wfn \
       --interaction emax2.MFDn --nrefmax 2 --hw 20 --output He4.jref
   ```

2. 同一个相互作用转成 imsrg++ 读的 J64（含完整 A 依赖真空 H）：

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
     mr_reference_file=He4.jref mr_validation_tolerance=1e-10 \
     smax=1e6 eta_criterion=1e-6 scratch=omega_scratch write_omega=true \
     Operators=Rp2,E2 flowfile=flow.dat intfile=He4 \
     write_H_jcoupled64=H.jcoupled64 write_H_no2bpack=H.no2bpack
   ```

   得到 `H.jcoupled64`、`H.no2bpack`、`H_Rp2.jcoupled64`、`H_Rp2.no2bpack`、`H_E2.opbin` 及各自的 `.mrmeta.json`。
   日志里 `MR initial natural-basis Hamiltonian: E=` 应等于参考态 NCSM 能量（这里 −17.01349 MeV）。
   `eta_criterion` 与 imsrg++ 一样是 ‖η‖ 的绝对值。He4 Nrefmax=2 emax2 验收时用的是 7e-7（‖η(0)‖ = 0.7036 的相对 1e-6），要约 8×10^5 步、
   s ≈ 4×10^5、53 段 Ω、单线程 44 h（慢尾，见 生成元.md）；试跑可先用 `smax=0.5`。

4. 下游读回：`mrimsrg_validate --no2bpack H.no2bpack --Z 2 --N 2 --nmax 10 --states 3`（见 驱动参数与输出.md）。

集群上用 `gen_job.py --mr-jscheme` 生成同样参数的 Slurm 作业（驱动参数与输出.md 第 5 节）。

## 文件位置

| 内容 | 位置 |
|---|---|
| 参考态类、读入与校验、MR 正规序 | `src/MRReference.{hh,cc}` |
| MR 对易子（标量与张量 λ2 项） | `src/MRCommutator.{hh,cc}` |
| BCH/Magnus 的 MR 分派与防护 | `src/BCH.{hh,cc}`、`src/IMSRGSolver.cc`（`Solve_ode_magnus`） |
| White-NCSM 生成元 | `src/Generator.cc`（`ConstructGenerator_WhiteNCSM`） |
| 驱动（参数检查、正规序、导出、算符） | `src/imsrg++.cc` |
| SR 张量对易子（MR 复用并做了提速） | `src/Commutator/TensorCommutators.cc` |
| jref 生成：simpleFCI 原型 / BIGSTICK / 嵌入 | `prototype/mrimsrg/export_jref.py`、`bigstick_jref.py`、`embed_jref.py` |
| 集群作业生成 | `gen_job.py`（`--mr-jscheme`、`--mr-prepare-jscheme`、`--mr-ncsm-readback`） |
| 测试 | `test/UnitTestMR*.py`（见 测试.md） |

## 各模块说明

- [参考态输入](参考态输入.md)：jref 格式与约定、三种生成方式、读入校验；
- [驱动参数与输出](驱动参数与输出.md)：MR 模式的参数、限制、日志、导出文件、下游与集群作业；
- [MR 对易子与正规序](MR对易子与正规序.md)：标量、张量的 λ2 项，MR 正规序，BCH/Magnus 中的 MR 处理与防护；
- [生成元](生成元.md)：White-NCSM 的定义、已知的慢尾问题；
- [算符演化](算符演化.md)：标量、张量算符经 Magnus Ω 的变换与导出；
- [性能与并行](性能与并行.md)：并行方式、缓存与环境变量、可复现性、内存；
- [测试](测试.md)：每项测试覆盖什么、怎么跑。
