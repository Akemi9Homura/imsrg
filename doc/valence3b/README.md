# 价空间三体（valence_3b_file）

分支 `vs-valence3b`，基于点7 部署版 0f2a34ed（含 Tokyo 算符写出的修复）。目的：VS-IMSRG(2) 流完后，用最低阶的办法生成价空间三体 W，
连同它在芯正规序下落下的低阶项一起写出，交给能带三体的价空间求解器（计划为 FCIQMC）。默认行为不变，只由参数 `valence_3b_file` 打开。

## 1. 参数

| 参数 | 取值 | 说明 |
|---|---|---|
| `valence_3b_file` | `none`（默认）或输出路径前缀 | 写 `<前缀>.snt`（H(s) 加上 W 的低阶项，相对芯正规序）和 `<前缀>_3b.int`（价空间三体）；照常的 `<intfile>.snt` 不变，作对照 |
| `valence_3b_source` | `perturbative`（默认）/ `flow` | `perturbative`：IMSRG(2) Magnus 流后由 Ω 生成；`flow`：取 IMSRG(3) 流（`IMSRG3=true`）中 H(s) 自己的三体 |
| `valence_3b_omega` | `gathered`（默认）/ `segments` | 微扰来源下，先把各段 Ω 合成一个，还是逐段累加，见第 3 节 |

要求 `method=magnus`。Ω 留在内存或写到 scratch 都可以，后者按 `Transform_Partial` 的文件名读回（`IMSRGSolver::GetOmegaSegment`）。

## 2. 做法

**W 的定义**（Johnson，Flow & Tell 2026，报告第 22 页）：

W = [Ω, H̃] 的三体部分，H̃ = H₀ + ½[Ω, H₀]（只取一体、两体），只保留六个指标都在价空间的元素。

这是 BCH 展开中三体在最后一次对易子里产生的那几项；W 自身再被 Ω 转动（[Ω, W] 等）不算。实现：`Commutator::comm223ss` 加已有的
`imsrg3_only_vvv` 开关；三体存储用 `ModelSpace::SetThreeBodyOrbits(valence)` 限制在价空间轨道上，E3max 不足 3·e_max(价) 时调大。

**为什么只要 vvv|vvv**：W 相对参考态（系综）正规序，求解器用相对芯的正规序。两者的密度差只在价空间里，芯的收缩在两种正规序下相同、相互抵消，
所以 P 空间里只有六个指标都在价空间的元素起作用。

**换回芯正规序**：三体矩阵元不变，落下的低阶项是用价空间部分的参考态密度（价空间里占据 > 0 的轨道）收缩出来的：

- 两体 −Σ_c n_c W_{abc,dec}，一体 +½ Σ n n W，零体 −⅙ Σ n n n W；
- 用 `DoNormalOrdering3(-1, valence_holes)` 一次算出（它只返回零、一、两体部分）。不能分两步"先对参考态撤销、再对芯正规序"，因为中间三体已经丢掉，
  芯的收缩就抵消不了；只算 vvv 时两种做法一样，取 IMSRG(3) 的 H(s) 时只有一次算出才对。

## 3. 两种 Ω 取法

Be8 emax2 上（第 5 节）两者差几十倍。芯解耦那一步的 Ω 只有"两个芯指标"的分量（cc→vv、qv、qq），在那一步生成的三体都带芯指标；
它们要到价空间那一步才被转成 vvv。逐段累加只算 vvv，这部分丢了；合成一个 Ω 时，BCH 的交叉项 ½[Ω_价, Ω_芯] 在二阶上包含它。
合成 Ω 的结果与 Johnson 报告的 Be8 emax2 数接近，所以默认用 `gathered`。

## 4. 代码改动

| 文件 | 改动 |
|---|---|
| `src/Parameters.hh` | 三个参数 |
| `src/imsrg++.cc` | 流后、换回芯正规序前算 W；换回芯正规序后写 `<前缀>.snt`、`<前缀>_3b.int` |
| `src/ModelSpace.{hh,cc}` | `SetThreeBodyOrbits`：限制张成三体空间的轨道并重建三体态 |
| `src/IMSRGSolver.{hh,cc}` | `GetNumberOfOmegaSegments`、`GetOmegaSegment`：按序号取 Ω 段，scratch 与内存都支持 |
| `src/ThreeBodyChannel.cc` | `GetLocalIndex`：不在三体空间里的态返回 −1（调用方本来就跳过 −1），不再抛 `map::at` |

## 5. 三体文件格式

`ReadWrite::WriteValence3body`：注释头列出轨道编号（质子在前、中子在后，从 1 起，与 `.snt` 同一套规则，但读的时候按 (n, l, 2j, 2tz) 对应更稳），
每行 `a b c d e f Jab Jde 2J V`，a ≤ b ≤ c、d ≤ e ≤ f，同轨道对的奇 Jab 不写。V 是反对称 m-scheme 元的不归一化 J 耦合：先 (ab) 耦合到 Jab，
再与 c 耦合到 J，四个 CG 系数（`ThreeBodyME::GetME_pn` 的约定），与 M 无关，6 位小数。展开回 m-scheme 的参照实现见 `tools/valence3b/pspace.py`
（耦合态写在规范行列式上，用伪逆反解，重复轨道自动处理）。

## 6. 检验（`tools/valence3b/check_be8.py`，本机约 5 秒）

N2LO_opt bare ħω20，emax2（me2j 文件截断），只用 NN，BetaCM = 0，basis=HF，He4 芯 + p 壳。`pspace.py` 是只作检验用的稠密 m-scheme 求解器。

| 检验 | 结果 |
|---|---|
| Ω 写 scratch 与留在内存 | `v3b.snt`、`v3b_3b.int` 逐字相同 |
| 不含 W 的 `.snt` 由 `pspace.py` 对角化 | E −26.3458、2⁺ 4.4066、⟨Rp2⟩ 4.1926，与 KSHELL/BIGSTICK 一致 |
| imsrg++ 写出的低阶项（两个 `.snt` 之差）与 m-scheme 公式 | Be8、C12 一致到 1e-6（文件精度） |
| C12：系综参考态即 p3/2⁴ × p3/2⁴ 行列式，⟨Φ₀\|:W:\|Φ₀⟩ | 2e-6 MeV，为零 |

Be8 基态中各项（MeV；精确 E −26.325、2⁺ 3.340、⟨Rp2⟩ 4.209）：

| | 三体 | 一体 + 两体 | 零体 | 合计 | 对角化后 E、2⁺ 的变化 |
|---|---|---|---|---|---|
| gathered | −0.0246 | +0.0910 | −0.0061 | +0.0603 | +0.060、−0.018 |
| segments | +0.0004 | +0.0008 | −0.0003 | +0.0009 | +0.001、−0.000 |
| Johnson 报告，Be8 emax2 | −0.0064 | +0.0930 | −0.0049 | 约 +0.08 | |

emax2 下这一项很小，符合预期（Johnson：emax10 下各项几 MeV、相消后约 1 MeV）。这里只说明工作流正确，不说明物理效应大小。
`--with-imsrg3` 可以加跑 emax2 的 IMSRG(3) 流（`valence_3b_source=flow`），本机两线程要几个小时，而且 emax2 下两者都很小，默认不跑。

## 7. 限制与下一步

- 只有 H 的三体；Rp2 等算符的三体部分没有做（要用它，求解器得能算三体期望值）。
- 裸 3N 的剩余部分 :V3:（NO2B 丢掉的正规序三体）没有包括，需要价空间那一块的完整 3N 矩阵元（NuHamil 生成）。
- 下一步：Ca 四轨道（NAT8）重跑并写出 W；先做三体单极诊断（每个跨壳核子的代价、4p4h 相对 0p0h、能隙随中子数的变化），
  再做一阶能量（KSHELL 密度），最后 FCIQMC 带三体求解。
