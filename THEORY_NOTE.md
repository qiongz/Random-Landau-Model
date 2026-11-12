# Random-Landau-Model 理论推导笔记

> 对应代码仓库：[Random-Landau-Model](https://github.com/qiongz/Random-Landau-Model)  
> 主要参考文献：Q. Zhu *et al.*, Phys. Rev. B **99**, 024205 (2019) [[arXiv:1804.00398](https://arxiv.org/abs/1804.00398)]  
> 早期方法：Y. Huo and R. N. Bhatt, Phys. Rev. Lett. **71**, 2376 (1993)

---

## 1. 物理背景与目标

### 1.1 整数量子 Hall  plateau 过渡

强磁场下二维电子气的 Hall 电导量子化为 $\sigma_{xy} = \nu e^2/h$。调节磁场或 Fermi 能 $E_F$，系统可在不同 plateau 之间发生**整数 quantum Hall 相变**。

相变由**关联长度**（局域长度）$\xi(E)$ 控制：

$$
\xi(E) \sim |E - E_c|^{-\nu}, \tag{1}
$$

其中 $E_c$ 为临界点能量，$\nu$ 为**局域长度临界指数**（localization-length exponent）。非相互作用电子的理论共识约为 $\nu \approx 2.4$，但 CCN 网络模型等近期工作给出更大值；Zhu *et al.* (2019) 用本仓库实现的两个投影模型重新标度，得到 $\nu = 2.48 \pm 0.02$。

### 1.2 程序要算什么

1. 在**无序**下对角化投影哈密顿量，得到本征能量 $E_m$ 和本征态 $|\psi_m\rangle$。
2. 对每个本征态计算 **Chern 数** $C(m)$——边界条件平均的 Hall 电导（拓扑不变量）。
3. **$C(m) \neq 0$ 的态**判为**扩展（conducting）态**；$C(m) \approx 0$ 为**局域态**。
4. 对有限尺寸系统统计 conducting 态数目 $N_c$，做 finite-size scaling 提取 $\nu$。

Chern 数方法的关键优势（Arovas *et al.*）：扩展态的波函数节点对 twist 边界条件敏感，局域态则不敏感——无需 transfer matrix 的一维 crossover。

---

## 2. 公共框架：Twist 边界条件

对有限 $L_1 \times L_2$ 样品，引入**广义边界条件**（magnetic translation with twist）：

$$
\hat{t}(L_j \hat{e}_j)\,|\psi\rangle = e^{i\theta_j}|\psi\rangle, \qquad \theta_j \in [0, 2\pi). \tag{2}
$$

参数 $\boldsymbol{\theta} = (\theta_1, \theta_2)$ 张成一个二维环面 $T^2$。对每个固定的 $\boldsymbol{\theta}$，哈密顿量 $H(\boldsymbol{\theta})$ 有 $N_\phi$ 个本征态；随 $\boldsymbol{\theta}$ 平滑变化时，每个本征态携带一个 U(1) 联络。

**代码实现**：主程序在 $(\theta_1, \theta_2)$ 网格上扫描，网格大小为 `n_mesh`（即论文中的 $L_g$），格点 $\theta_j = 2\pi \cdot l_j / n_\text{mesh}$。

---

## 3. 模型 A：连续随机朗道模型（RLM）

对应目录 `rlm/`，可执行文件 `rlm_cpu_dp` / `rlm_gpu_magma` 等。

### 3.1 哈密顿量

对称规范 $\mathbf{A} = \frac{B}{2}(-y, x, 0)$，$\hbar = e = c = 1$，磁长度 $l_B = 1$：

$$
H = \frac{1}{2}(\boldsymbol{\Pi} - \mathbf{A})^2 + V(\mathbf{r}), \tag{3}
$$

其中 $\boldsymbol{\Pi} = \mathbf{p} + \mathbf{A}$ 为 kinetic momentum。定义**磁平移生成元**：

$$
\boldsymbol{\kappa} = \boldsymbol{\Pi} - \mathbf{A} = \boldsymbol{\Pi}(-B), \qquad [\kappa_1, \kappa_2] = -i. \tag{4}
$$

样品内总磁通量子数 $N_\phi = L_1 L_2 / (2\pi l_B^2)$（代码中 `-n` 即 $N_\phi$）。当 $N_\phi$ 为整数时，$\hat{t}(L_1\hat{e}_1)$ 与 $\hat{t}(L_2\hat{e}_2)$ 对易。

### 3.2 最低朗道能级（LLL）基

定义原胞磁平移算符 $t_j = e^{i\kappa_j L_j / N_\phi}$，构造 stripe 基 $| \phi_m \rangle$（Landau 规范下）：

$$
t_1 |\phi_m\rangle = |\phi_{m+1}\rangle, \qquad
t_2 |\phi_m\rangle = e^{-i 2\pi m / N_\phi} |\phi_m\rangle. \tag{5}
$$

$|\phi_m\rangle$ 构成 LLL 的一组完备正交基，$m = 0, 1, \ldots, N_\phi - 1$。

### 3.3 无序势：$\delta$ 点杂质

$$
V(\mathbf{r}) = \sum_{i=1}^{N_\text{imp}} \alpha_i\, \delta(\mathbf{r} - \mathbf{r}_i), \tag{6}
$$

- 杂质数 $N_\text{imp} = \texttt{impurity\_concentration} \times N_\phi / \texttt{quanta\_concentration}$（见 `main_rlm_cpu.cpp`）。
- 位置 $\mathbf{r}_i = (x_i, y_i)$ 在 $[0, L_1] \times [0, L_2]$ 均匀随机；强度 $\alpha_i = \pm 1$ 等概率（`disorder_potential.cpp`）。
- 论文取 $N_\text{imp} = 16 N_\phi$；代码默认 `impurity_concentration = 16`。

Fourier 分量：

$$
V_{\mathbf{q}} = \sum_i \alpha_i\, e^{-i \mathbf{q} \cdot \mathbf{r}_i}, \qquad
\mathbf{q} = \left(\frac{2\pi m}{L_1}, \frac{2\pi n}{L_2}\right). \tag{7}
$$

### 3.4 投影到 LLL

将无序势投影到 LLL 子空间：

$$
\tilde{V} = \sum_{j,k} |\phi_j\rangle\langle\phi_j| \sum_{\mathbf{q}} V_{\mathbf{q}}\, e^{i\mathbf{q}\cdot\mathbf{r}} |\phi_k\rangle\langle\phi_k|. \tag{8}
$$

利用 LLL 基下平面波矩阵元（Landau 形式因子），可得动量空间系数 $V_{mn}$。代码中 `get_potential_coeff` 计算：

$$
V_{mn} = \underbrace{\sum_i \alpha_i\, e^{-i(k_x x_i + k_y y_i)}}_{\text{impurity 求和}}
\times \underbrace{\frac{1}{N_\phi}\, e^{-\frac{\pi}{2N_\phi}\left(m^2\frac{L_2}{L_1} + n^2\frac{L_1}{L_2}\right)}\,
e^{-i\pi mn/N_\phi}}_{\texttt{coeff\_mn}[m,n]}, \tag{9}
$$

其中 $k_x = 2\pi(m - m_0)$，$k_y = 2\pi(n - n_0)$，$m_0 = n_0 = \texttt{off\_head} = \lfloor\sqrt{N_\text{imp}}/2\rfloor + 1$ 为 K 空间截断中心。动量截断 $|m - m_0|, |n - n_0| \le \texttt{off\_head}$，矩阵维度 $\dim_m = \dim_n = 2\cdot\texttt{off\_head} + 1$。

> **代码对应**：`matrix_coefficients.cpp` 中 `prepare_coeff` 预计算 `coeff_mn`、`coeff_m_theta`、`coeff_jm`；每个 disorder sample 调用 `get_potential_coeff` 得到 `v_mn`。

### 3.5 含 twist 的投影哈密顿量矩阵元

在 twist 边界条件 $(\theta_1, \theta_2)$ 下，投影哈密顿量在 LLL 基 $| \phi_j \rangle$ 上的矩阵元为（论文 Eq. (12)）：

$$
H_{jk}(\theta_1, \theta_2)
= \sum_{m,n} V_{mn}\,
e^{-\frac{\pi}{2N_\phi}\left(m^2\frac{L_2}{L_1} + n^2\frac{L_1}{L_2}\right)}\,
e^{-i\pi mn/N_\phi}\,
e^{i(m\theta_2 - n\theta_1)/N_\phi}\,
e^{-i 2\pi m j / N_\phi}\,
\langle \phi_j | \phi_{k-n} \rangle. \tag{10}
$$

磁平移代数给出 Kronecker $\delta$ 型选择定则：$\langle \phi_j | \phi_{k-n} \rangle$ 仅在 $k - n \equiv j \pmod{N_\phi}$ 时非零。等价地，对每个固定的 $m$，只有**至多两个** $n$ 值贡献：

$$
n = \texttt{off\_head} + (k - j), \quad \text{或} \quad
n = \texttt{off\_head} + (k - j) - N_\phi,
$$

当相应指标落在 $[0, \dim_n)$ 内时。

将 Eq. (10) 分解为代码中的三块预计算系数：

| 系数 | 公式 | 代码变量 |
|------|------|----------|
| 形式因子 | $\frac{1}{N_\phi} e^{-\pi(m^2 L_2/L_1 + n^2 L_1/L_2)/(2N_\phi)} e^{-i\pi mn/N_\phi}$ | `coeff_mn[m,n]` |
| $\theta_2$ 相位 | $e^{i(m-m_0)\theta_2/N_\phi}$ | `coeff_m_theta[θ₂,m]` |
| 基函数重叠 | $e^{-i 2\pi (m-m_0)(j+1)/N_\phi}$ | `coeff_jm[j,m]` |
| $\theta_1$ 相位 | $e^{i(n-n_0)\theta_1/N_\phi}$ | `coeff_m_theta[θ₁,n]`（取共轭） |

于是

$$
H_{jk}(\theta_1, \theta_2)
= \sum_m \underbrace{\overline{\texttt{coeff\_m\_theta}[\theta_2,m]}\,\texttt{coeff\_jm}[j,m]}_{\texttt{coeff\_theta\_jm}}
\sum_{n \in \mathcal{S}(j,k)} \texttt{v\_mn}[m,n]\,\overline{\texttt{coeff\_m\_theta}[\theta_1,n]}, \tag{11}
$$

其中 $\mathcal{S}(j,k)$ 为上述 Kronecker 选择定则允许的两个 $n$ 值集合。

**Hermiticity**：代码只算 $k \ge j$ 的上三角元，再令 $H_{kj} = H_{jk}^*$。

**对角化**：对每个 $(\theta_1, \theta_2)$，$H$ 为 $N_\phi \times N_\phi$ Hermitian 矩阵，调用 MKL `heevd` / MAGMA `cheevd_gpu` 求本征值和本征矢（`hamiltonian_rlm_cpu.cpp`，`mkl_diag.cpp`）。

### 3.6 粒子-空穴对称

短程无序 + 投影到 LLL 后，能谱关于 $E_c = 0$ 对称（论文 Sec. II B）。临界能量固定在零，便于标度分析。

---

## 4. 模型 B：Hofstadter 晶格模型（HLM）

对应目录 `hlm/`，可执行文件 `hlm_cpu_dp`。

### 4.1 晶格哈密顿量

$$
H = -\sum_{\langle i,j\rangle} t\, e^{i\theta_{ij}} c_i^\dagger c_j + \text{h.c.} + \sum_i \epsilon_i\, c_i^\dagger c_i, \tag{12}
$$

其中 $\theta_{ij} = \frac{2\pi}{q} \times (\text{沿 hops 累积的磁通相位})$，每 unit cell 磁通 $\phi = \phi_0/q$，$q$ 即代码参数 `-q`（默认 $q=3$）。

对 $q=3$，干净 Hamiltonian 有 3 个子带，Chern 数分别为 $+1, -2, +1$。程序投影到**最低子带** `-b 0`，并将 hopping 设 $t=0$、无序宽度 $W=1$（论文：$t/W \le 0.2$ 时临界行为与 $t$ 无关）。

### 4.2 子带投影

对每个 twist 点 $(\theta_1, \theta_2)$：

1. **求干净子带波函数** $\psi_{n,\mathbf{k}}(\mathbf{r})$：对角化 $q \times q$ 的 k-space Hofstadter 哈密顿量（`get_wfs_ktheta`）。
2. **乘无序势**：$\tilde{\psi}_{n,\mathbf{r}} = \epsilon(\mathbf{r})\, \psi_{n,\mathbf{k}}(\mathbf{r})$。
3. **投影矩阵元**：

$$
H_{nm} = \sum_{\mathbf{r}} \psi_{n,\mathbf{k}}^*(\mathbf{r})\, \epsilon(\mathbf{r})\, \psi_{m,\mathbf{k}}(\mathbf{r})
= \langle \psi_n | V | \psi_m \rangle. \tag{13}
$$

代码用 MKL `gemm3m` 实现：$H = \Psi^\dagger V \Psi$，其中 $\Psi$ 为 $N_\text{sites} \times N_\phi$ 的干净波函数矩阵（`hamiltonian_hlm_cpu.cpp`）。

4. 对角化 $N_\phi \times N_\phi$ 的 $H$，再变换回实空间存储完整波函数供 Chern 计算：

$$
|\Psi_m(\boldsymbol{\theta})\rangle = \sum_n \Psi_n \cdot U_{nm}. \tag{14}
$$

### 4.3 与 RLM 的关系

附录 A 表明：小磁通极限下 HLM 平滑过渡到连续 LLL 模型。两个模型给出一致的 $\nu$，互相独立验证。

---

## 5. Chern 数计算（Fukui–Hatsugai–Suzuki 方法）

### 5.1 定义

第 $m$ 个本征态在边界条件空间的**第一 Chern 数**等于平均 Hall 电导：

$$
\langle \sigma_{xy}(m) \rangle = \frac{e^2}{h}\, C(m), \tag{15}
$$

$$
C(m) = \frac{1}{2\pi i} \int_0^{2\pi}\!\!\int_0^{2\pi}
\left\langle \frac{\partial \psi_m}{\partial \theta_x} \Big| \frac{\partial \psi_m}{\partial \theta_y} \right\rangle
- \text{h.c.}\;\; d\theta_x\, d\theta_y. \tag{16}
$$

### 5.2 离散化：link variable 与 plaquette

将 $T^2$ 离散为 $L_g \times L_g$ 格点（代码 `n_mesh`）。定义 U(1) link：

$$
U_\mu^m(\boldsymbol{\theta}) = \frac{\langle \psi_m(\boldsymbol{\theta}) | \psi_m(\boldsymbol{\theta} + \hat{\mu}) \rangle}
{\big|\langle \psi_m(\boldsymbol{\theta}) | \psi_m(\boldsymbol{\theta} + \hat{\mu}) \rangle\big|}, \tag{17}
$$

其中 $\hat{\mu}$ 为 $\theta_x$ 或 $\theta_y$ 方向的格矢（长度 $2\pi/L_g$）。

每个 plaquette 的 Berry 相位：

$$
\gamma_m(\boldsymbol{\theta}) = \frac{1}{i} \ln
\Big[ U_x^m(\boldsymbol{\theta})\, U_y^m(\boldsymbol{\theta}+\hat{x})\,
\big(U_x^m(\boldsymbol{\theta}+\hat{y})\big)^{-1}\,
\big(U_y^m(\boldsymbol{\theta})\big)^{-1} \Big]. \tag{18}
$$

Chern 数：

$$
C(m) = \frac{1}{2\pi} \sum_{\boldsymbol{\theta}} \gamma_m(\boldsymbol{\theta}). \tag{19}
$$

### 5.3 代码实现（`chern.cpp`）

对每个本征态指标 $n = 0, \ldots, N_\phi - 1$，在 $\theta_2$ 方向的每个格点 $m$ 上：

1. 取 plaquette 四个顶点的波函数 $|\psi(\theta_1, \theta_2)\rangle$。
2. 沿四边计算 overlap $S_\ell = \langle \psi_l | \psi_r \rangle$（MKL `zdotc`/`cdotc`）。
3. 累乘 $\prod_{\ell=1}^{4} S_\ell$，取辐角：

$$
\gamma_m = -\frac{\arg\big(\prod_{\ell=1}^{4} S_\ell\big)}{2\pi}. \tag{20}
$$

4. 对所有 plaquette 求和，并加 $1/N_\phi$ 归一化（`main_rlm_cpu.cpp` 第 237–238 行）。

**拓扑求和规则**：每个 disorder sample 所有本征态 Chern 数之和为 $+1$（最低子带/LLL 的总拓扑 charge）。程序以此为 consistency check。

**扩展 vs 局域**：
- $C(m) \approx +1$：扩展（conducting）态
- $C(m) \approx 0$：局域态
- 偶尔 $C(m) \approx -1$：成对出现的 conducting pair（见论文 Fig. 6）

**RLM 与 HLM 符号**：HLM 在累加 Chern 贡献时多一个负号（`main_hlm_cpu.cpp` 第 186 行），源于波函数存储约定不同，物理结果一致。

---

## 6. 有限尺寸标度

### 6.1 Conducting 态密度

$$
N_c = \int_{-\infty}^{\infty} \rho_c(E)\, dE, \tag{21}
$$

其中 $\rho_c(E)$ 为 $C(m) \neq 0$ 的本征态密度。对每个 disorder sample，$N_c$ 即 conducting 态的个数。

### 6.2 标度律

有限尺寸 $L \sim \sqrt{N_\phi}$ 下，只有局域长度 $\xi(E) > L$ 的态才是 conducting。在 $E_c$ 附近：

$$
N_c \sim L^{2}\, \rho(E_c)\, |E_m - E_c| \sim N_\phi^{\,1 - 1/(2\nu)}. \tag{22}
$$

**纯幂律拟合**（论文 Eq. (22)）：

$$
N_c = a\, N_\phi^{\,1 - 1/(2\nu)}. \tag{23}
$$

**含 irrelevant 修正**（晶格模型，$y \approx 4.3$）：

$$
N_c = a\,\big(1 + b\, N_\phi^{-y/2}\big)\, N_\phi^{\,1 - 1/(2\nu)}. \tag{24}
$$

Zhu *et al.* 的主要结果：

| 模型 | 方法 | $\nu$ |
|------|------|-------|
| Hofstadter 晶格 | $N_c$ 标度 + irrelevant 修正 | $2.479 \pm 0.005$ |
| 连续 LLL | $N_c$ 标度（无修正） | $2.46 \pm 0.01$ |
| 合并 | — | **$2.48 \pm 0.02$** |

### 6.3 为何 Chern 数方法对标度扰动不敏感

晶格模型 irrelevant 指数 $y = 4.3 \gg 1$，说明拓扑数计算受 irrelevant 算符影响很弱——这解释了 Huo–Bhatt (1993) 在小系统上就能得到有效 $\nu$ 的原因，也解释了为何 CCN 模型需要谨慎处理 marginal irrelevant 算符。

---

## 7. 数值流程与代码映射

```
┌─────────────────────────────────────────────────────────────┐
│  初始化：n_phi, n_mesh, 线程数, 随机种子                      │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  prepare_coeff：预计算 coeff_mn, coeff_m_theta, coeff_jm     │
│  （与 disorder 无关，只做一次）                               │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
          ┌────────────────┴────────────────┐
          │  对每个 disorder sample n       │
          └────────────────┬────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  生成无序势 → 计算 v_mn（CPU 或 CUDA potential_coeff.cu）    │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  θ₁ = 0, 1, …, n_mesh：                                      │
│    对每个 θ₂ 并行（pthread）：                                │
│      组装 H(θ₁,θ₂) → 对角化 → 存波函数 + 能量               │
│    若 θ₁ ≥ 1：并行计算 Chern 贡献                           │
└──────────────────────────┬──────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  对 θ 平均能量；累加 Chern → 输出 energy_Chern_No.txt        │
└─────────────────────────────────────────────────────────────┘
```

| 物理量 | 代码位置 |
|--------|----------|
| 形式因子 $e^{-\pi(m^2+n^2)/(2N_\phi)} e^{-i\pi mn/N_\phi}/N_\phi$ | `matrix_coefficients.cpp` → `coeff_mn` |
| twist 相位 $e^{i m \theta / N_\phi}$ | `coeff_m_theta` |
| 基重叠 $e^{-i2\pi j m/N_\phi}$ | `coeff_jm` |
| 无序 $V_{mn}$ | `get_potential_coeff` / `potential_coeff.cu` |
| $H_{jk}(\theta_1,\theta_2)$ | `hamiltonian_rlm_cpu.cpp` / `hamiltonian_rlm_gpu.cu` |
| 对角化 | `mkl_diag.cpp` / `magma_diag.cu` / `cusolver_diag.cu` |
| Chern 数 | `chern.cpp` → `peer_cal_Chern` |
| HLM 子带投影 | `hamiltonian_hlm_cpu.cpp` → `get_wfs_ktheta` + `gemm3m` |

---

## 8. 精度与性能注意事项

| 问题 | 说明 |
|------|------|
| 单精度 Chern 误差 | $N_\phi > 1500$ 时 Chern 数严重失真；**标度分析必须用 `rlm_cpu_dp`（双精度）** |
| 单精度用途 | 生成波函数供机器学习（`rlm_gpu_magma_wfs` + `DwfsIO`） |
| 计算瓶颈 | 对角化占 ~97% 时间（README benchmark） |
| Chern 网格 | 论文取 $L_g = 50$（连续）或 $L_g \approx \sqrt{2/3}\,L$（晶格）；代码 `-g` 即 `n_mesh` |
| 总 Chern 检验 | $\sum_m C(m) = 1$；失败则 reject 整个 sample（论文 Sec. III C） |

---

## 9. 主要公式速查

**RLM 投影 Hamiltonian（核心）**：

$$
\boxed{
H_{jk}(\theta_1,\theta_2) = \sum_{m,n} V_{mn}\,
e^{-\frac{\pi}{2N_\phi}(m^2 + n^2)}\,
e^{-i\pi mn/N_\phi}\,
e^{i(m\theta_2 - n\theta_1)/N_\phi}\,
e^{-i2\pi m j/N_\phi}\,
\delta_{j,\, k-n \;(\text{mod } N_\phi)}
}
$$

（正方形样品 $L_1 = L_2$，$m,n$ 相对于截断中心 $m_0$ 定义。）

**Chern 数（离散 Fukui 方法）**：

$$
\boxed{
C(m) = \frac{1}{2\pi} \sum_{\text{plaquettes}} \arg \prod_{\ell=1}^{4} \langle \psi_m^{(\ell)} | \psi_m^{(\ell+1)} \rangle
}
$$

**标度律**：

$$
\boxed{
N_c \sim N_\phi^{\,1 - 1/(2\nu)}, \qquad \nu \approx 2.48
}
$$

---

## 10. 参考文献

1. Q. Zhu, P. Wu, R. N. Bhatt, and X. Wan, Phys. Rev. B **99**, 024205 (2019). [arXiv:1804.00398](https://arxiv.org/abs/1804.00398)
2. Y. Huo and R. N. Bhatt, Phys. Rev. Lett. **71**, 2376 (1993).
3. T. Fukui, Y. Hatsugai, and H. Suzuki, J. Phys. Soc. Jpn. **74**, 1674 (2005).
4. D. P. Arovas, R. N. Bhatt, F. D. M. Haldane, P. B. Littlewood, and R. R. P. Singh, Phys. Rev. Lett. **60**, 1695 (1988).
