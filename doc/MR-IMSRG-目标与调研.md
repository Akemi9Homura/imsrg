# MR-IMSRG 工作目标与现状调研

## 1. 这项工作的目的

本仓库从 `im-fciqmc-hf-export` 的已验证代码基线分出，用于独立推进 MR-IMSRG 方向，不与原 IMSRG 导出链路的维护混在一起。

这项工作服务于 IM-FCIQMC 的 B 路线。B 路线要回答的问题是：经过 in-medium 变换后，Hamiltonian 是否能在不破坏目标物理的前提下显著改善 FCIQMC 对 walker population 的需求，使原本难以收敛的大体系在可接受的 $N_w$ 下收敛。

已经确认，完全解耦的单参考 IMSRG(2) Hamiltonian 不能用于检验这一问题。对于参考行列式 $D_0$，二体 Hamiltonian 只能直接连接 1p1h 和 2p2h 行列式，而这些矩阵元正是单参考流要消去的部分。因此，物化的 rank-2 $H(s)$ 会使 $D_0$ 结构性孤立，随后进行 FCIQMC 只能得到平凡的单行列式结果，不能说明一般意义上的软化是否改善 $N_w$ 收敛。

MR-IMSRG 的候选价值在于：用相关的多行列式态或参考空间作为变换参考，解耦对象不再是单个 $D_0$。如果截断后的 Hamiltonian 仍保留参考空间内部以及它与其余行列式之间的非平凡结构，就可能成为检验“软化与 FCIQMC population 收敛关系”的有效代理。

这里的 MR-IMSRG 物化测试仍然只是 B 路线的诊断工具。IM-FCIQMC 的最终目标仍是处理未物化的 $U^\dagger H(0)U$，从而恢复原始 Hamiltonian 的目标 FCI 能量；不能把截断 MR-IMSRG Hamiltonian 的本征值直接当作最终答案。

## 2. 真正的 MR-IMSRG 指什么

MR-IMSRG 相对于普通单参考 IMSRG 的关键变化不是“参考态含有多个行列式”这一句描述，而是算符代数必须相对于相关参考态进行广义正规序。

相关参考态的收缩不能只用单轨道占据数 $n_p$ 表示。除了常规的一体密度，还需要不可约密度矩阵（density cumulants），其中 MR-IMSRG(2) 流方程显式依赖

$$
\lambda^{(2)},\qquad \lambda^{(3)}.
$$

当参考态退化为单个 Slater determinant 时，$\lambda^{(2)}$ 和 $\lambda^{(3)}$ 消失，MR 方程才退化为单参考 IMSRG 方程。因此，只有分数占据数、natural occupations 或 ensemble occupations，并不自动构成文献中 IM-NCSM/IM-GCM 所使用的完整 MR-IMSRG。

算符中的三体矩阵元与参考态的三体不可约密度也是两类不同对象。启用 IMSRG(3)、`ThreeBodyME` 或更高的流动算符 rank，不会自动补出 $\lambda^{(3)}$，因而不能把单参考代码变成 MR-IMSRG。

## 3. IM-NCSM 与 IM-GCM 中的用法

### IM-NCSM

IM-NCSM 先在较小的 NCSM 参考空间中求得相关的多行列式参考态及其密度矩阵，再以该参考态进行 MR-IMSRG 演化，最后在更大的 NCSM 空间中对演化后的 Hamiltonian 对角化。

对闭壳体系也不能只把 `Nmax_ref=0` 的单个闭壳行列式换一个名称叫作 MR。如果希望得到非平凡的相关参考态，参考空间本身必须允许产生多行列式成分，例如使用非零的参考空间截断。

### IM-GCM

IM-GCM 使用投影 HFB/GCM 态或相应的参考系综构造相关参考，并在 MR-IMSRG 中保留相关参考的不可约二体、三体密度信息。GCM 参考态的具体构造与本仓库当前阶段可以分开研究，但它进一步证明了：文献中的多参考 in-medium 变换依赖相关密度，而不只是改变一体占据数。

## 4. 当前 `imsrg++` 代码的能力边界

当前代码的官方能力范围是 Hartree-Fock、single-reference IMSRG 和 valence-space IMSRG，未声明支持 MR-IMSRG。

代码中的参考态数据主要是每个轨道的 `occ` 和 `occ_nat`。正规序和 commutator 使用 $n_p$、$n_pn_q$ 一类占据因子；没有发现用于输入、存储和收缩相关参考态 $\lambda^{(2)}$、$\lambda^{(3)}$ 的数据结构，也没有接受 CI/NCSM/GCM 参考波函数或 2-RDM/3-RDM 的接口。

代码确实包含以下相邻能力，但它们不能等同于真正的 MR-IMSRG：

- fractional occupation 和 ensemble normal ordering；
- targeted normal ordering；
- valence-space IMSRG；
- natural-orbital occupation；
- 二体、三体流动算符及 IMSRG(3) 的部分实现。

其中 ensemble normal ordering 只给出一体占据层面的混合参考描述。在当前实现里，更高阶参考密度被占据数乘积所替代，无法复现包含 $\lambda^{(2)}$、$\lambda^{(3)}$ 的 MR-IMSRG(2) 方程。因此，当前代码不存在一个输入参数或开关可以直接完成 IM-NCSM/IM-GCM 式 MR 计算。

## 5. 开源实现调研

截至 2026-08-24，没有检索到公开、可直接读入核相互作用并完成数值 MR-IMSRG 演化的通用核结构程序。已找到的公开项目如下。

### IMSRG++

[官方 IMSRG++ 仓库](https://github.com/ragnarstroberg/imsrg) 是本仓库的上游来源。其公开说明只列出 HF、single-reference IMSRG 和 valence-space IMSRG，与本地源码检查一致。

### QCombo

[QCombo](https://github.com/chenlh73/qcombo) 是面向多体算符 commutator 的开源符号工具，实现 generalized Wick theorem，并提供 `MR_IMSRG2.ipynb` 示例。它可以生成 MR-IMSRG(2) 相关的 LaTeX 表达式和供 AMC 进行 J-scheme 化简的输入。

QCombo 不是数值 MR-IMSRG 求解器：它不负责读取核力、构造参考态密度、积分流方程或输出可直接对角化的演化 Hamiltonian。它的价值是作为公式生成和独立代数校验工具。

### ADG

[ADG](https://github.com/adgproject/adg) 可以自动生成和求值多体图，包括 BIMSRG 图。它同样主要服务于图和代数表达式生成，不是可直接运行核结构 MR-IMSRG 计算的数值实现；BIMSRG 也不是本文所讨论的粒子数守恒 MR-IMSRG。

### 其他相邻实现

量子化学中存在开源的 multireference DSRG 实现，可以提供广义正规序和多参考流方法的软件设计参考，但其积分表示、对称性、参考态接口和核结构 J-coupled 约定均与当前程序不同，不能直接替代核 MR-IMSRG。

## 6. 当前结论与工作边界

1. MR-IMSRG 是值得独立验证的 B 路候选，因为它可能避免单参考 rank-2 $H(s)$ 将 $D_0$ 结构性孤立的问题。
2. 现有 `imsrg++` 不能直接运行真正的 MR-IMSRG；fractional occupation、ENO、VS-IMSRG 或更高算符 rank 都不能替代相关参考的 $\lambda^{(2)}$、$\lambda^{(3)}$。
3. 当前没有找到可直接复用的公开数值核 MR-IMSRG 求解器；论文工作所用数值实现至少没有以容易发现和直接使用的形式公开。
4. QCombo 和 ADG 可用于公式、图和约定检查，但不能直接产出演化后的数值 Hamiltonian。
5. 本文档只固定研究目的、物理判据和已确认的软件现状，不预先决定具体实现路径。是否获取作者代码、扩展当前 J-scheme 程序，或先建立独立的小空间原型，应在后续工作中根据可复用资源和验证成本决定。

## 7. 主要资料

- H. Hergert, *In-Medium Similarity Renormalization Group for Closed and Open-Shell Nuclei*, [arXiv:1607.06882](https://arxiv.org/abs/1607.06882)：广义正规序、不可约密度和 MR-IMSRG 方程综述。
- E. Gebrerufael et al., *Ab Initio Description of Open-Shell Nuclei: Merging No-Core Shell Model and In-Medium Similarity Renormalization Group*, [arXiv:1610.05254](https://arxiv.org/abs/1610.05254)：小参考空间、多参考流与后续 NCSM 对角化的组合。
- J. M. Yao et al., *Ab Initio Treatment of Collective Correlations and the Neutrinoless Double Beta Decay of $^{48}$Ca*, [arXiv:1908.05424](https://arxiv.org/abs/1908.05424)：投影 GCM 参考态与 MR-IMSRG 的结合。
- C. Li et al., QCombo 项目及文档：[GitHub](https://github.com/chenlh73/qcombo)、[PyPI](https://pypi.org/project/qcombo/)。
- A. Tichai et al., ADG/BIMSRG 项目：[arXiv:2102.10889](https://arxiv.org/abs/2102.10889)、[GitHub](https://github.com/adgproject/adg)。
