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
