# NNLOopt emax2 MR-IMSRG 验收结果

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
  `O(lambda2)` 项，保持分母符号并用 `1e-6 MeV` cutoff。
- `Rgen` 是 Vobig 式 (6.5.28)--(6.5.29) 中带 EN 分母、反对称化后
  的实际掩码生成元范数，是正式固定点判据。现有 `IMSRGSolver`
  监测同一个 `Eta.Norm()`，但使用绝对阈值；本原型按照计划采用更严格的
  相对初值 `1e-6` 门槛。`Rnum` 是不带分母的
  lambda-free White-NCSM 分子；`Rstrict` 是含线性 `lambda2`
  的 `D-D^dagger` 诊断。对分数占据，分母加权与反对称化不对易，
  因此三者始终分别报告。
- 放松 IM-NCSM 掩码只消去 `Delta e != 0` 的 1B/2B 通道；同 HO
  量子数的参考空间内部块保留，因此最终零体项不必等于后
  NCSM 基态能。

## 2. 实现正确性证据

- Python 回归测试：`42/42` 通过。覆盖 RDM/cumulant、普通与 MR 正规序
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
  `no2bpack` 原生 reader 得到 `-20.3396333376 MeV`，相差 `0.942 keV`；
  这是既有格式用 float32 保存 OBME/TBME 的量化舍入误差。
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

## 4. `Nrefmax=2` 相关闭壳参考态

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

## 5. 文件与来源标识

第一批四核流位于 point7：

```text
/tns/mengziyan/mr-imsrg/result/mrimsrg-flow/<Nucleus>_Nrefmax0_rtol1em06_atol1em08_white_ncsm_strict/flow
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
起点。新生产作业的 commit、job id 与能量将在完成后补入本文。
