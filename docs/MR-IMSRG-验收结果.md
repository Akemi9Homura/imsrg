# NNLOopt emax2 MR-IMSRG 验收结果

本文记录快速 m-scheme MR-IMSRG(2) 原型的可重复验收证据。大型
Hamiltonian、RDM 和日志均保留在被 Git 忽略的本机/集群结果目录，
不提交到仓库。表中能量单位均为 MeV。

## 1. 固定输入与数值定义

- 相互作用：NNLOopt（文件名 `N2LO_opt`），`hw=20`, `emax=2`,
  `e2max=4`，仅 NN，无 Coulomb/显式 3N。
- 文件 SHA-256：
  `76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84`。
- MR-IMSRG(2) 对易子保留线性 `lambda2`，设 `lambda3=0`；直接积分
  `dH/ds=[eta,H]_(0,1,2B)`。
- 生成元是 Vobig Sec. 6.5.4 定义的 `White-NCSM`：分子舍去所有
  不可约密度项，Epstein--Nesbet 分母舍去文献标记的
  `O(lambda2)` 项，保持分母符号并用 `1e-6 MeV` cutoff。
- `Rgen` 是 Vobig 式 (6.5.28)--(6.5.29) 中带 EN 分母、反对称化后
  的实际掩码生成元范数，是正式固定点判据，与现有
  `IMSRGSolver` 的 `Eta.Norm()` 停止语义一致。`Rnum` 是不带分母的
  lambda-free White-NCSM 分子；`Rstrict` 是含线性 `lambda2`
  的 `D-D^dagger` 诊断。对分数占据，分母加权与反对称化不对易，
  因此三者始终分别报告。
- 放松 IM-NCSM 掩码只消去 `Delta e != 0` 的 1B/2B 通道；同 HO
  量子数的参考空间内部块保留，因此最终零体项不必等于后
  NCSM 基态能。

## 2. 实现正确性证据

- Python 回归测试：`36/36` 通过。覆盖 RDM/cumulant、普通与 MR 正规序
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
- 四个 `Nrefmax=0` 流的一体 Hermiticity、二体 Hermiticity 和二体
  反对称误差在所有保存步均为数值零。

## 3. 第一批：`Nrefmax=0` 四核结果

这四个作业使用 `rtol=1e-6`, `atol=1e-8`。`He4/O16` 是
`lambda2=0` 的 SR-control；`Be8/C12` 是非零 `lambda2` 的真正 MR 流。
`Be8/C12` 作业还使用了比最终程序更严格的 `Rstrict<=1e-6`
停止门槛，因此同时满足已发表的 `Rgen` 标准。

| 核 | `||lambda2||F` | `s_final` | `Rgen/Rgen0` | `Rnum/Rnum0` | `Rstrict/Rstrict0` | MR 0B | 后 NCSM 基态 | 报告 `Nmax` | `s=0` 同空间 | 本征值变化 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| He4 | 0 | 87.5627 | 9.089e-7 | 8.015e-7 | 8.015e-7 | -20.2399184011 | -20.2399185504 | 8 | -20.3388325043 | +0.0989139539 |
| Be8 | 2.439406 | 341.4413 | 6.183e-7 | 6.268e-7 | 8.987e-7 | -27.6449723642 | -27.9391351682 | 4 | -25.1799853478 | -2.7591498204 |
| C12 | 2.185990 | 273.7649 | 2.116e-7 | 2.342e-7 | 9.693e-7 | -43.3080111039 | -43.5312362939 | 4 | -40.2376145487 | -3.2936217452 |
| O16 | 0 | 82.0767 | 7.639e-7 | 7.703e-7 | 7.703e-7 | -63.7200911877 | -63.7200911998 | 4 | -62.2980557824 | -1.4220354174 |

“本征值变化”是同一表列 `Nmax` 下演化后减 `s=0`，它包含
MR-IMSRG(2) 截断与有限 NCSM 空间的联合影响，不得解读为精确幺正
变换应保持的全空间谱。

### NCSM 收敛加速

| 核 | `Nmax` | 维数 | `s=0` 基态 | 最终基态 |
|---|---:|---:|---:|---:|
| He4 | 0 | 1 | -13.7612333298 | -20.2399184011 |
| He4 | 2 | 59 | -17.0134907986 | -20.2399185157 |
| He4 | 4 | 720 | -20.0960773346 | -20.2399185481 |
| He4 | 8 | 3060 | -20.3388325043 | -20.2399185504 |
| Be8 | 0 | 51 | -7.4286548545 | -27.9391351673 |
| Be8 | 2 | 4523 | -20.1837628509 | -27.9391351682 |
| Be8 | 4 | 73886 | -25.1799853478 | -27.9391351682 |
| C12 | 0 | 51 | -22.8877671396 | -43.5312362938 |
| C12 | 2 | 15897 | -35.9965187284 | -43.5312362939 |
| C12 | 4 | 552290 | -40.2376145487 | -43.5312362939 |
| O16 | 0 | 1 | -46.4602252394 | -63.7200911877 |
| O16 | 2 | 1201 | -58.3574679785 | -63.7200911995 |
| O16 | 4 | 231651 | -62.2980557824 | -63.7200911998 |

`Be8/C12/O16` 演化后的 `Nmax=0/2/4` 基态已在约 `1e-8 MeV`
或更好的量级一致；`He4` 的 `Nmax=0` 到全 `Nmax=8` 变化为
`0.149 keV`。这是放松 `Delta e` 解耦已实际改善后 NCSM 收敛的直接
验收，而不是只看 MR 零体项。

### 低激发态读回示例

`mrimsrg_validate --states 3` 直接复用 `simpleFCI` 的多态 Lanczos 与
`J^2` 计算。例如：

- `He4, Nmax=8`：`E0=-20.2399185504`; 随后两态
  `(Ex,2J)=(3.9191826814,4)` 和 `(4.1621749620,0)`。
- `Be8, Nmax=4`：`E0=-27.9391351682`; 随后两态
  `(Ex,2J)=(4.6314479391,4)` 和 `(8.8723744548,0)`。
- `C12, Nmax=2`：`E0=-43.5312362939`; 随后两态
  `(Ex,2J)=(4.8766987676,4)` 和 `(9.6377808031,0)`。
- `O16, Nmax=2`：`E0=-63.7200911995`; 随后两态
  `(Ex,2J)=(8.8947907231,0)` 和 `(9.2666575955,4)`。

### ODE 容差验收

把 `rtol/atol` 从 `1e-6/1e-8` 收紧到 `1e-7/1e-9`：

- `He4, Nmax=8` 后 NCSM 基态变化 `0.203 keV`；
- `C12, Nmax=2` 后 NCSM 基态变化 `0.217 keV`。

两者都通过 `<1 keV` 门槛，同时覆盖 SR-control 和非平凡 MR 案例。

## 4. `Nrefmax=2` 相关闭壳参考态

`He4/O16, Nrefmax=2` 均有非零 `lambda2`，且 HO 基中 `gamma1`
非对角。按照 Vobig Sec. 6.5.3，公式和 `Delta e` 掩码均在临时连通块
自然轨道基中计算；自然轨道继承输出槽位的 `e` 标签，最终 Hamiltonian
再变回原 HO 轨道顺序。

| 核 | `||lambda2||F` | HO 基 `max|gamma1_offdiag|` | `s_final` | `Rgen/Rgen0` | `Rnum/Rnum0` | `Rstrict/Rstrict0` | MR 0B | 后 NCSM 基态 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| He4 | 0.598982 | 0.0722519 | 待填 | 待填 | 待填 | 待填 | 待填 | 待填 |
| O16 | 1.082447 | 0.0429950 | 待填 | 待填 | 待填 | 待填 | 待填 | 待填 |

这两个作业的正式验收使用 `Rgen`。`Rnum/Rstrict` 会保留为生成元截断
诊断；若它不趋零，结果只声称收敛到已发表 White-NCSM 近似的固定点。

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

`Nrefmax=0` 生产作业在数值代码修订 `4278fc8a` 上运行；最终
仓库只改变了停止语义、metadata 和多态 NCSM 报告，生成元/对易子
数值公式未变。`Nrefmax=2` 作业使用当前论文一致语义的代码。
