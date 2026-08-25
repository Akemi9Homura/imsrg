# MR-IMSRG(2) J-scheme 公式与 `imsrg++` 代码映射

## 1. 适用范围与来源

本表只覆盖当前生产目标：实数、粒子数守恒、`Jref=0` 标量参考态，
Hamiltonian 与生成元保留 0B/1B/2B，`lambda2` 保留而 `lambda3=0`。

公式以以下原文交叉核对：

- Hergert 2016 review，Eqs. (49)--(51)，m-scheme MR-IMSRG(2)；
- Gebrerufael 2017 学位论文，Eqs. (4.43)--(4.50)、(4.89)、
  (4.126)--(4.159)，标量 J-coupling、Pandya 变换和 BLAS 形式；
- Vobig 2020，Eqs. (6.5.13)--(6.5.16)、(6.5.26)--(6.5.34)，
  White-NCSM 掩码、分子和 EN 分母；
- Mongelli 2022，Eqs. (5.210)--(5.221)，修改后 White 生成元；
- `refs/qcombo/examples/MR_IMSRG2.ipynb` 的保存符号输出。

实现约定同时对照当前 `src/TwoBodyME.*`、`Commutator.cc`、
`Operator.cc`、`Generator.cc`、`IMSRGSolver.cc`，以及
`/home/mengziyan/fciqmc/fciqmc-mpi/src/libfci/fciqmc_obs.cpp` 的
`build_rdm2_mscheme_from_TBTD()`。不能脱离这些代码约定单独使用下式。

## 2. 二体矩阵与密度的冻结约定

Gebrerufael Eq. (4.45) 使用没有 identical-pair normalization 的
`^J O_abcd`。`imsrg++::TwoBodyME` 的矩阵中则保存正规化反对称 pair
矩阵元

\[
 O^J_{(ab)(cd)} =
 \frac{1}{\sqrt{1+\delta_{ab}}\sqrt{1+\delta_{cd}}}
 \frac{1}{2J+1}\sum_{M,m_a\ldots m_d}
 C^{JM}_{a m_a b m_b}C^{JM}_{c m_c d m_d}O_{m_a m_bm_cm_d}.
\]

`GetTBME_J()` 返回论文 flow 方程使用的未正规化矩阵元，即给存储值
乘回两个 `sqrt(1+delta)`；`GetTBME_J_norm()` 返回存储值。受限 pair
矩阵乘法直接使用存储矩阵，不再额外加入 `1/2`。

本项目的 `lambda2` 使用相同的正规化 pair block 存储。这时

\[
 \frac14\sum_{pqrs} O_{pqrs}\lambda^{pq}_{rs}
 =\sum_{J,\pi,T_z}(2J+1)
   \sum_{a\le b,c\le d}O^J_{ab,cd}\lambda^J_{ab,cd}.
\]

FCIQMC/shell-model-obs 的 rank-0 TBTD 是约化二体密度；其 m-scheme
重建式另含 `1/sqrt(2J+1)`。因此 reader 写入 `TwoBodyME` 前必须做

\[
 \lambda^J_{\texttt{TwoBodyME}}
 = \lambda^J_{\text{TBTD}}/\sqrt{2J+1}.
\]

不能把 TBTD 数值不经转换直接当作 Hamiltonian 型 TBME。

独立 oracle `prototype/mrimsrg/jcoupling.py` 不调用 `imsrg++`，而是显式
用 CG 系数完成上述投影。当前结果为：

| 输入 | `max|J->m - m|` | coupled/m-scheme 标量收缩差 |
|---|---:|---:|
| 随机严格标量小张量 | 回投影 `<=2e-14` | `<=1e-13` |
| `Be8, Nrefmax=0` 的 `lambda2` | `5.38e-12` | `7.11e-15 MeV` |
| `He4, Nrefmax=2` 的 `lambda2` | `1.29e-12` | `8.88e-16 MeV` |

真实态第一列是原 m-scheme RDM 的极小非标量数值残差；coupled block
回投影自身仍按随机严格标量门禁验收。reader 必须报告该投影残差，不能
静默丢弃超过容差的分量。

## 3. commutator contraction 表

令 `X=eta`、`Y=H`，一般接口仍实现 `[X,Y]`。Gebrerufael 的
`chi=-1` 对应 anti-Hermitian `eta`，所以 `[eta,H]` Hermitian。

| 输出 | m-scheme/QCombo | Gebrerufael J-scheme | 生产映射 |
|---|---|---|---|
| 2B | Hergert Eq. (51) | (4.158a) | 全部复用现有 `comm121ss`、`comm122ss`、`comm222_pp_hhss`、`comm222_phss`；无显式 `lambda2` |
| 1B 公共项 | Hergert Eq. (50) 前三组 | (4.117)--(4.125) | 复用 `comm111ss`、`comm121ss`、`comm221ss`，占据取自然轨道 `n_p` |
| 1B-IV | QCombo `1/4 X2 Y2 lambda2` direct | (4.126)--(4.127) | 每个普通 pp channel 的 `X^J lambda^J Y^J`，再以 `(2J+1)/(4(2j_1+1))` 缩并 spectator `t` |
| 1B-V | QCombo crossed `X2 Y2 lambda2` | (4.128)--(4.129) | 对 `X/lambda/Y` 分别作当前 Eq. (4.49) 符号的 Pandya 变换，在 cross-coupled channel 做三矩阵乘积 |
| 1B-VI | QCombo 最后两组 `1/2 X2 Y2 lambda2` | (4.130)--(4.149) | 先在普通 channel 形成 `lambda*Y`、`lambda*X`，再只取正宇称 `J=0` Pandya block，与 `Xbar/Ybar` 相乘 |
| 0B 公共项 | Hergert Eq. (49) 前两组 | (4.150)--(4.155) | 复用 `comm110ss`、`comm220ss` |
| 0B MR | `1/4 C2 lambda2` | (4.156)--(4.158c) | 在完整 2B commutator 已累加后做 `(2J+1) trace(C2^J lambda2^J)`；不能只收缩 2B--2B 子项 |

1B-IV--VI 最终都必须加入论文的 `-chi[1<->2]`。生产代码应通过
`Operator` 的 Hermitian/anti-Hermitian 属性填写另一半矩阵，不能在三组
内部各自采用不一致的手写共轭规则。

Hergert Eq. (50) 的四个显式 m-scheme `lambda2` 指标式，在上述
J-coupling 后合并成 Gebrerufael 的 IV、V、VI 三种矩阵拓扑；不是少了一项。
第一版 C++ 应先保留三种独立 profiler/test 名称，最后才允许融合中间量。

按 `TwoBodyME` normalized-pair block 数值实现 Eq. (4.89b) 的未限制
spherical-orbit 慢循环时，IV、V 在最终 `1<->2` permutation 之前的系数
分别为 `1/8`、`1/2`，VI 保持 `-1/2`。这两个额外的 `1/2` 是论文
unnormalized/unrestricted pair 写法转为当前存储后的结果，不是物理截断。
独立 CG 展开已把三种拓扑分别对到 m-scheme；C++
`MRCommutator::comm221_lambda2_reference()` 又分别对到该 Python J-scheme
oracle，max-abs 依次为 `0`、`0`、`2.84e-14`。生产实现
`MRCommutator::comm221_lambda2()` 使用普通 pair block 和有符号
`Delta Tz` 的有序 particle-hole block 做矩阵乘法；它对慢速 C++ reference
的 IV、V、VI max-abs 分别为 `1.07e-14`、`3.55e-14`、`5.68e-14`。
生产 RHS 中不构造 m-scheme 张量。

需要区分论文前因子和归一化 block 乘积前因子。IV 的慢速、未限制
spherical-orbit 指标式在最终 permutation 前是 `1/8`；生产代码先以
normalized-pair block 形成 `X*lambda*Y`，其中对中间相同粒子 pair 的
限制求和已经吸收另一组归一化/计数因子，所以 spectator contraction 的
系数是 `1/2`。VI 中完整未限制 `(r,v)` 求和等于 normalized-pair
`lambda*Y` block 乘积的两倍，因此矩阵乘积实现使用 `-1`，而慢速指标式
保持 `-1/2`。这些不是改变物理公式；逐拓扑的 J-to-m 和慢/快双重比较
是当前系数的验收依据。

完整入口 `MRCommutator::Commutator()` 先调用原
`Commutator::Commutator()`，只加入上述 1B 项及完整 2B 输出与
`lambda2` 的 0B 收缩。`lambda2.Norm()==0` 有显式直接返回门禁，随机
算符测试中 MR 与原 SR 输出的 0B/1B/2B 差为逐位 `0`。

## 4. 现有代码可直接复用的边界

- `ModelSpace::SetupKets()` 已根据 orbit `occ` 构造分数占据的 pp/hh/ph
  权重；公共 MR 项不需要重写。
- `Commutator::DoPandyaTransformation*()` 使用与 Gebrerufael Eq. (4.49)
  相同的整体负号和 `(2J+1) 6j`。MR-V/VI 的新矩阵布局必须复用或逐元素
  对照这一实现。
- `Operator::DoNormalOrdering2()` 已完成 `gamma1=n_p delta_pq` 收缩。
  NN-only 的 MR 正规序只需额外给 ZeroBody 加
  `sum_J(2J+1) V^J lambda^J`；反正规序减去同一项。
- `Generator::Get1bDenominator()`、`Get2bDenominator()` 已实现当前
  scalar EN monopole 约定。现有 cutoff 对 `|Delta|<cutoff` 一律置为
  **正** `cutoff`；为满足严格 SR 退化，MR 生产入口首版必须复用这一行为，
  不得单独改成保号 cutoff。若以后修正，应作为 SR/MR 共同的独立变更并
  重做全部退化门禁。
- `IMSRGSolver` 已有直接 RK4 与 adaptive ODE。首版只增加显式 MR
  commutator context，不复制 solver。

## 5. 实现及验收顺序

1. 用 `jcoupling.py` 冻结真实 `gamma2/lambda2` 的 block、相位和收缩；
2. 新增显式 `MRReference`，reader 后直接做 coupled RDM 恒等式与标量
   投影残差检查；
3. 先实现慢速、只供测试的 J-scheme IV/V/VI 指标循环，对 m-scheme
   命名项逐元素验证；
4. 再用普通/Pandya channel 矩阵乘法实现生产 IV/V/VI，并逐项对慢速
   J reference；
5. 增加完整 `C2-lambda2` 0B 收缩；
6. `lambda2=0` 时确认 dispatcher 只执行原 SR contractions，然后进入
   White-NCSM generator 与 flow 集成。

截至 2026-08-25，第 1--6 步中的 commutator 部分均已实现并通过；下一
边界是 White-NCSM generator 与 solver/driver 显式上下文接线。

## 6. 生产参考态文件闭环

`prototype/mrimsrg/export_jref.py` 复用现有 simpleFCI/shell-model-obs
波函数桥、`compute_densities()` 和已验收的自然基/J-coupling oracle，写出
`mrimsrg_jref_v1`。文件使用 little-endian float64，包含来源 SHA-256、
`A/Z/Nrefmax/J2/parity/hw/emax/e2max`、完整 `(n,l,j2,tz2)` 轨道表、
自然占据、原 HO 到内部自然 J-orbit 的正交变换，以及 normalized-pair
`lambda2` blocks。C++ 按量子数映射轨道和 pair，不依赖两端原始编号。

reader 在进入流之前检查 `J=0+`、空间参数、自然变换的正交性和球形 block
结构、占据范围与粒子/质子 trace、`lambda2` Hermiticity 及 cumulant 收缩。
`Be8 Nrefmax=0` 的 17 KiB 实文件读回后，来源哈希与既有验收记录一致，
cumulant 收缩误差 `6.44e-15`，全部 block 对 Python CG 投影逐元素误差
`<1e-13`。`He4 Nrefmax=2` 的非单位径向自然轨道变换也通过读回与收缩检查。

## 7. White-NCSM 生成元的生产归约

`Generator(type="white-ncsm")` 实现 Vobig Sec. 6.5.4 中舍去所有
irreducible-density 项的方向分子，同时保留自然占据产生的 norm factors。
令 `bari=1-n_i`，则球形一体方向量为

```text
D1(i,j) = bari*n_j*f(i,j)
Delta1(i,j) = -bari^2*n_j^2*GammaMono(i,j)
              +bari^2*n_j*f(i,i) - bari*n_j^2*f(j,j)
              +E*(bari*n_j-1).
```

二体方向量为

```text
w = bari*barj*n_k*n_l
D2(ij,kl) = w*Gamma(ij,kl)
Delta2(ij,kl) = w*[bari*barj*GammaMono(i,j)
                    +n_k*n_l*GammaMono(k,l)
                    -bari*n_l*GammaMono(i,l)
                    -bari*n_k*GammaMono(i,k)
                    -barj*n_k*GammaMono(j,k)
                    -barj*n_l*GammaMono(j,l)
                    +bari*f(i,i)+barj*f(j,j)
                    -n_k*f(k,k)-n_l*f(l,l)] + E*(w-1).
```

生产代码先对每个方向应用与现有 `Generator.cc` 相同的正 `1e-6 MeV`
cutoff，再作 `eta=D/Delta-(D/Delta)^dagger`；不能先反厄米化分子再除一个
公共分母。只保留 `2n_i+l_i != 2n_j+l_j` 和
`e_i+e_j != e_k+e_l`。Slater He4 随机算符对原 `white` 入口的完整差
`<1e-12`；分数占据 Be8 随机算符对 Python m-scheme 的全部 1B/J-coupled
2B block 误差 `<2e-11`。

随机严格标量张量要求 J/m 转换 `<=1e-12`；真实 RDM 的输入标量投影
残差单独报告并暂以 `1e-10` 为拒绝阈值。所有能量/RHS contraction 的
coupled 对 m-scheme 误差仍要求 `<=1e-10 MeV`。
