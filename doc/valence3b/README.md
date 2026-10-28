# 价空间三体（valence_3b_file）

分支 `vs-valence3b`，基于点7 部署版 0f2a34ed（含 Tokyo 算符写出的修复）。目的：VS-IMSRG(2) 流完后，用最低阶的办法生成价空间三体 W，
连同它在芯正规序下落下的低阶项一起写出，交给能带三体的价空间求解器（计划为 FCIQMC）。默认行为不变，只由参数 `valence_3b_file` 打开。

## 1. 参数

| 参数 | 取值 | 说明 |
|---|---|---|
| `valence_3b_file` | `none`（默认）或输出路径前缀 | 写 `<前缀>.snt`（H(s) 加上 W 的低阶项，相对芯正规序）和 `<前缀>_3b.int`（价空间三体）；照常的 `<intfile>.snt` 不变，作对照 |
| `valence_3b_source` | `perturbative`（默认）/ `flow` / `at_end` | `perturbative`：IMSRG(2) Magnus 流后由 Ω 生成；`flow`：取 IMSRG(3) 流（`IMSRG3=true`）中 H(s) 自己的三体；`at_end`：只作诊断，见第 3 节 |
| `valence_3b_omega` | `segments`（默认）/ `gathered` | 微扰来源下，逐段累加，还是先把各段 Ω 合成一个，见第 3 节 |

要求 `method=magnus`。Ω 留在内存或写到 scratch 都可以，后者按 `Transform_Partial` 的文件名读回（`IMSRGSolver::GetOmegaSegment`）。

## 2. 做法

**W 的定义**：每一段 Ω_k 上 W += [Ω_k, H̃_k] 的三体部分，H̃_k = H_{k−1} + ½[Ω_k, H_{k−1}]（只取一体、两体，H_{k−1} 是前 k−1 段变换后的 H），
只保留六个指标都在价空间的元素。这是流方程 dW/ds = [η, H]₃ 的最低阶积分 W = ∫ds [η(s), H(s)]₃ 的离散形式，即 V 的二阶；
W 自身再被转动（[η, W]₃）属于三阶，不算。单个 Ω 时就是 Johnson（Flow & Tell 2026，报告第 22 页）的 [Ω, H₀ + ½[Ω, H₀]]₃。实现：`Commutator::comm223ss` 加已有的
`imsrg3_only_vvv` 开关；三体存储用 `ModelSpace::SetThreeBodyOrbits(valence)` 限制在价空间轨道上，E3max 不足 3·e_max(价) 时调大。

**为什么只要 vvv|vvv**：W 相对参考态（系综）正规序，求解器用相对芯的正规序。两者的密度差只在价空间里，芯的收缩在两种正规序下相同、相互抵消，
所以 P 空间里只有六个指标都在价空间的元素起作用。

**换回芯正规序**：三体矩阵元不变，落下的低阶项是用价空间部分的参考态密度（价空间里占据 > 0 的轨道）收缩出来的：

- 两体 −Σ_c n_c W_{abc,dec}，一体 +½ Σ n n W，零体 −⅙ Σ n n n W；
- 用 `DoNormalOrdering3(-1, valence_holes)` 一次算出（它只返回零、一、两体部分）。不能分两步"先对参考态撤销、再对芯正规序"，因为中间三体已经丢掉，
  芯的收缩就抵消不了；只算 vvv 时两种做法一样，取 IMSRG(3) 的 H(s) 时只有一次算出才对。

## 3. Ω 的取法与 at_end

- **segments（默认）**：逐段累加，是上面的二阶结果，与 Ω 怎样分段无关（每段都小时就是那个积分）。
- **gathered**：先用截断到两体的 BCH 乘积把各段合成一个 Ω（与 hunter-gatherer 相同），再套单个 Ω 的公式。多出的
  [[Ω_价, Ω_芯], H]₃ 等是三阶项，而且只带上了一部分（"芯那一步生成、带芯指标的三体被价空间那一步转进 vvv"这一项是三阶，合并只借 Jacobi 恒等式近似地带上），
  结果依赖于表示；Ω 范数大时 BCH 乘积本身也不稳。曾因与 Johnson 报告的 Be8 emax2 数接近而设为默认，这个理由不成立，已改回 segments。
- **at_end（诊断）**：照 `imsrg3_at_end`，用 IMSRG(2) 的各段 Ω 以 IMSRG(3) 对易子重做 e^Ω H e^{−Ω}，三体在段间保留。
  Ω 是按 IMSRG(2) 截断求的，这个变换不再解耦：Be8 emax2 p 壳上价空间生成元的 ‖H_od‖ 从 H(s) 的 0.000 变成 21.2，
  P 块不再是有效的价空间哈密顿量（基态能量少束缚 3 MeV）。只作诊断。
- **Be8 p 壳为什么 segments 几乎为零**：p 壳 emax2 里价空间之外只有 sd 壳（和 0s 芯），价→外的单粒子跃迁改变宇称，
  vv↔vq 两体耦合被禁戒，二阶的 vvv 三体基本为零，剩下的都是三阶。所以 Be8 p 壳不适合检验 W 的大小。Ca 四轨道的价空间同时有 sd 和 pf，
  vv↔vq 是允许的，二阶 W 不为零。
- 真正的对照是完整的 IMSRG(3) 流方程（`valence_3b_source=flow`），它的 Ω 是带着三体求的，三阶以上自洽。

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
| segments | +0.0004 | +0.0008 | −0.0003 | +0.0009 | +0.001、−0.000 |
| gathered | −0.0246 | +0.0910 | −0.0061 | +0.0603 | +0.060、−0.018 |
| at_end（不解耦，只作诊断；另有两体回馈 E +3.92、2⁺ −0.80） | +0.625 | −1.867 | +0.394 | −0.848 | +3.065、−0.599（含两体回馈） |
| Johnson 报告，Be8 emax2 | −0.0064 | +0.0930 | −0.0049 | 约 +0.08 | |

at_end 的低阶项检验和 C12 行列式检验同样通过（1e-6），说明输出流程正确，问题在变换本身不解耦。

emax2 下这一项很小，符合预期（Johnson：emax10 下各项几 MeV、相消后约 1 MeV）。这里只说明工作流正确，不说明物理效应大小。
`--with-imsrg3` 可以加跑 emax2 的 IMSRG(3) 流（`valence_3b_source=flow`），本机两线程要几个小时，而且 emax2 下两者都很小，默认不跑。

组态重心（`tools/valence3b/centroid.py`，第 1 级诊断用）：对给定组态，把 H 在该组态全部行列式上的对角元取平均；k 体项等于
"k 粒子子集数 × 该轨道组合上的平均对角元"，不用求解器。C12（参考态 p3/2⁸ 行列式）上 W 的修正：参考组态 0，质子 1p1h 0，
质子 2p2h +0.014 MeV（与对 6 个行列式的暴力平均一致），4p4h +0.066 MeV。相对行列式正规序的三体项，对角元至少要三个准粒子才不为零，
所以 1p1h 组态上严格为零。Be8 的系综参考态就是 p3/2 组态的平均，对应组态上也为零。

推论：⁴⁰Ca 参考态是闭壳行列式，W 在一阶上不改变一粒子一空穴层面的量（质子能隙、3⁻），只作用于多粒子多空穴组态；
开壳核的系综参考态带分数占据，组态涨落使 W 也能改变那里的跨壳代价。

## 7. 限制与下一步

- 只有 H 的三体；Rp2 等算符的三体部分没有做（要用它，求解器得能算三体期望值）。
- 裸 3N 的剩余部分 :V3:（NO2B 丢掉的正规序三体）没有包括，需要价空间那一块的完整 3N 矩阵元（NuHamil 生成）。
- 下一步：完整 IMSRG(3) 流方程的对照（Be8 p 壳，以及一个两种宇称都有的小价空间，二阶 W 在那里不被禁戒），本机太慢，放点7 的 v3b worktree；
  确定近似可用后，再上 Ca 四轨道（NAT8）：三体单极诊断 → 一阶能量（KSHELL 密度）→ FCIQMC 带三体求解。
