# fSSE推导多采样周期线性混杂系统一般形式

## 1.数字采样控制系统状态量$\vec{x}(T)$的表达式

对于如下线性系统：

$$
\frac{d \vec{x}}{dT} = \mathcal{A} \vec{x} + \mathcal{B} \vec{u}
$$

其中$\vec{x}$是状态量，$\vec{x} \in \mathbb{R}^{n}$。

采用有时滞的数字采样控制，系统的输入与信号之间为线性关系，即

$$
\begin{aligned}
    \vec{u} &= K_{0}(\mathrm{Gain}) \vec{x}_{N} + K_{1}(\mathrm{Gain}) \vec{x}_{N - 1} + \cdots + K_{M}(\mathrm{Gain}) \vec{x}_{N -M} \\
    &= K(\mathrm{Gain}) \Theta_{N},\qquad T \in [N, N + 1), N \in \mathbb{Z}^{+}
\end{aligned}
$$

其中，$\mathrm{Gain}$是表示增益组合的元组，不是某一个矩阵。$\Theta_{N}$由一系列采样时刻的状态量拼接而来：

$$
\Theta_{N} = \begin{pmatrix}
    \vec{x}_{N} \\ \vec{x}_{N - 1} \\ \vdots \\ \vec{x}_{N - M}
\end{pmatrix}
$$

则状态空间方程可写为：

$$
\frac{d \vec{x}}{dT} = \mathcal{A} \vec{x} + \mathcal{B}_{0} \vec{x}_{N} + \mathcal{B}_{1} \vec{x}_{N - 1} + \cdots +\mathcal{B}_{M} \vec{x}_{N - M},
\qquad T \in [N, N + 1), N \in \mathbb{Z}^{+}
$$

其中，系统矩阵$\mathcal{A}$和输入矩阵$\mathcal{B}$与系统自身有关，$\mathcal{B}_{i}(i \in \mathbb{Z}^{+} \cap [1,M],M \in \mathbb{Z}^{+})$与输入矩阵$\mathcal{B}$以及增益的组合$\mathrm{Gain}$有关。

令

$$
\begin{aligned}
    A_{0}(T_{0},\mathrm{Gain}) &= \exp{(\mathcal{A} T_{0 })}, \\
    {}^{temp}B_{0}(T_{0},\mathrm{Gain}) &= \int_{0}^{T_{0}} \exp{(\mathcal{A} t_{0})} d t_{0} \mathcal{B}_{0} \\
    &= \mathcal{A}^{-1} \bigg(\exp{(\mathcal{A} T_0)} - I \bigg) \mathcal{B}_{0}, \\
    B_{0}(T_{0},\mathrm{Gain}) &= A_{0}(T_{0}) + {}^{temp}B_{0}(T_{0}) \\
    B_{1}(T_{0},\mathrm{Gain}) &= \mathcal{A}^{-1} \bigg(\exp{(\mathcal{A} T_{0})} - I \bigg) \mathcal{B}_{1}, \\
    &\vdots \\
    B_{M}(T_{0},\mathrm{Gain}) &= \mathcal{A}^{-1} \bigg(\exp{(\mathcal{A} T_{0})} - I \bigg) \mathcal{B}_{M}. \\
\end{aligned}
$$

显然，$B_{i}(i \in \{0,1,2,\cdots,M\})$是满秩的。

给定增益组合后，上述的矩阵都只是关于$T_{0}$的函数，可以计算得到

$$
\begin{aligned}
    \vec{x}_{N + 1} &= A_{0}(1) \vec{x}_{N} + {}^{temp}\!B_{0}(1) \vec{x}_{N} + B_{1}(1) \vec{x}_{N - 1} + \cdots + B_{M}(1) \vec{x}_{N - M} \\
    &= (A_{0}(1) + {}^{temp}\!B_{0}(1)) \vec{x}_{N} + B_{1}(1) \vec{x}_{N - 1} + \cdots + B_{M}(1) \vec{x}_{N - M} \\
    &= B_{0}(1) \vec{x}_{N} + B_{1}(1) \vec{x}_{N - 1} + \cdots + B_{M}(1) \vec{x}_{N - M}
\end{aligned}
$$

令

$$
B(T_{0}) = \begin{bmatrix}
        B_{0}(T_{0}) & B_{1}(T_{0}) & \cdots & B_{M}(T_{0})
    \end{bmatrix}
$$

显然，

$$
\begin{aligned}
    A_{0}(0) = I_{n \times n}, B_{0}(0) = A_{0}(0) = I, B_{1}(0) = \cdots = B_{M}(T_{0}) = \mathbf{0} \\
    \Rightarrow B(0) = \begin{bmatrix}
        I_{n \times n} & \mathbf{0} & \cdots & \mathbf{0}
    \end{bmatrix} = \begin{bmatrix}
        I_{n \times n} & \mathbf{0}_{n \times Mn}
    \end{bmatrix}
\end{aligned}
$$

有

$$
\begin{aligned}
    \Theta_{N + 1} &= \begin{pmatrix}
        \vec{x}_{N + 1} \\ \vec{x}_{N} \\ \vdots \\ \vec{x}_{N - M + 1}
    \end{pmatrix} \\
    &= \begin{bmatrix}
        B_{0}(1) & B_{1}(1) & \cdots & \cdots & B_{M}(1) \\
        I & \mathbf{0} & \cdots & \cdots & \mathbf{0} \\
        \mathbf{0} & I & \mathbf{0} & \cdots & \mathbf{0} \\
        \mathbf{0} & \mathbf{0} & \ddots & \ddots & \mathbf{0} \\
        \mathbf{0} & \mathbf{0} & \cdots & I & \mathbf{0}
    \end{bmatrix} \Theta_{N} \\
    &= \left [ \begin{array}{cc}
        \begin{matrix}
            B_{0}(1) & B_{1}(1) & \cdots & \cdots & B_{M}(1)
        \end{matrix} \\
        \hdashline
        \begin{array}{c:c}
            I_{Mn \times Mn} & \mathbf{0}_{Mn \times n}
        \end{array}
    \end{array} \right ]_{(M + 1)n \times (M + 1)n} \Theta_{N} \\
    &= \left [ \begin{array}{cc}
        B(1) \\
        \hdashline
        \begin{array}{c:c}
            I_{Mn \times Mn} & \mathbf{0}_{Mn \times n}
        \end{array}
    \end{array} \right ]_{(M + 1)n \times (M + 1)n} \Theta_{N} \\
    &= A \Theta_{N} \\
    &= A(\mathrm{Gain}) \Theta_{N}
\end{aligned}
$$

显然，$A$是满秩的，即$A$是$\mathbb{R}^{(M + 1)n} \mapsto \mathbb{R}^{(M + 1)n}$的满射。

根据$\Theta$的定义，有

$$
\vec{x}_{N} = \begin{bmatrix}
    I_{n \times n} & \mathbf{0}_{n \times n} & \cdots & \mathbf{0}_{n \times n}
\end{bmatrix} \Theta_{N} = \begin{bmatrix}
    I_{n \times n} & \mathbf{0}_{n \times Mn}
\end{bmatrix} \Theta_{N} = B(0) \Theta_{N}
$$

类似于之前的求解过程，对于$T = N + T_{0} (N \in \mathbb{N}, T_{0} \in (0,1))$，可以计算得到

$$
\begin{aligned}
    \vec{x}_{N + T_{0}} &= B_{0}(T_{0}) \vec{x}_{N} + B_{1}(T_{0}) \vec{x}_{N - 1} + \cdots + B_{M}(T_{0}) \vec{x}_{N - M} \\
    &= \begin{bmatrix}
        B_{0}(T_{0}) & B_{1}(T_{0}) & \cdots & B_{M}(T_{0})
    \end{bmatrix} \begin{pmatrix}
        \vec{x}_{N} \\ \vec{x}_{N - 1} \\ \vdots \\ \vec{x}_{N - M}
    \end{pmatrix} \\
    &= \begin{bmatrix}
        B_{0}(T_{0}) & B_{1}(T_{0}) & \cdots & B_{M}(T_{0})
    \end{bmatrix} \Theta_{N} \\
    &= B(T_{0}) \Theta_{N}
\end{aligned}
$$

可以计算得到，

$$
\begin{aligned}
    \vec{x}_{1} &= B_{0}(1) \vec{x}_{0} \\
    \vec{x}_{2} &= B_{0}(1) \vec{x}_{1} + B_{1}(1) \vec{x}_{0} \\
    &= (B_{0}^{2}(1) + B_{1}(1)) \vec{x}_{0} \\
    \vec{x}_{2} &= B_{0}(1) \vec{x}_{2} + B_{1}(1) \vec{x}_{1} + B_{2}(1) \vec{x}_{0} \\
    &= (B_{0}^{3}(1) + B_{0}^{2}(1)B_{1}(1) + B_{0}(1)B_{1}(1) +B_{1}^{2}(1) + B_{2}(1)) \vec{x}_{0} \\
    & \vdots
\end{aligned}
$$

可见，$\Theta_{M} = \begin{pmatrix}
    \vec{x}_{M} & \vec{x}_{M - 1} & \cdots & \vec{x}_{0}
\end{pmatrix}^{\top}$是$\vec{x}_{0}$在$\mathbb{R}^{(M + 1)n}$中的线性映射，该映射只与增益的组合$\mathrm{Gain}$有关，表示为

$$
\vec{\Theta}_{M} = C(\mathrm{Gain}) \vec{x}_{0}
$$

显然，$C = C(\mathrm{Gain})$是一个$\mathbb{R}^{n} \mapsto \mathbb{R}^{(M + 1)n}$的满射。

因此，对于$T = n + T_{0}(n \in \mathbb{Z}^{+},T_{0} \in [0,1))$，有

$$
\vec{x}_{T} = \vec{x}_{n + T_{0}} = B(T_{0}) A^{N - M} C \vec{x}_{0} = D(N, T_{0}, \mathrm{Gain}) \vec{x}_{0}
$$

此外，

$$
\begin{aligned}
    \frac{d \vec{x}_{T}}{d T} &= \mathcal{A} \vec{x}_{T} + \mathcal{B} \Theta_{N} \\
    &= \mathcal{A} B(T_{0}) A^{N - M} C \vec{x}_{0} + \mathcal{B} A^{N - M} C \vec{x}_{0} \\
    &= \bigg(\mathcal{A} B(T_{0}) + \mathcal{B} \bigg) A^{N - M} C \vec{x}_{0} \\
    &= E(N, T_{0}, \mathrm{Gain}) \vec{x}_{0}
\end{aligned}
$$

## 2.CLF函数指数稳定约束的引入

CLF函数的形式为

$$
V(\vec{x}_{T}) = \frac{1}{2} \vec{x}_{T}^{\top} H \vec{x}_{T}
$$

其中$H$是个正定的对称矩阵。要求做到指数收敛：

$$
\frac{dV}{dT} + \lambda V < 0
$$

相当于找到$\mathrm{Gain}$的取值范围，使得$\forall \vec{x}_{0} \in \mathbb{R}^{n}, N \in \mathbb{Z}^{+}, T_{0} \in [0,1)$，都有

$$
\begin{aligned}
    \frac{d \vec{x}_{T}^{\top}}{d T} \frac{H}{2} \vec{x}_{T} + \vec{x}_{T}^{\top} \frac{H}{2} \frac{d \vec{x}_{T}}{d T} + \frac{\lambda}{2} \vec{x}_{T}^{\top} H \vec{x}_{T} &< 0 \\
    \vec{x}_{0}^{\top} E^{\top} H D \vec{x}_{0} + \vec{x}_{0}^{\top} D^{\top} H E \vec{x}_{0} + \lambda \vec{x}_{0}^{\top} D^{\top} H D \vec{x}_{0} &< 0 \\
    \vec{x}_{0}^{\top} \bigg(E^{\top} H D + D^{\top} H E + \lambda D^{\top} H D \bigg) \vec{x}_{0} &< 0 \\
    \vec{x}_{0}^{\top} C^{\top} (A^{\top})^{N - M} \Psi(T_{0}, \mathrm{Gain}) A^{N - M} C \vec{x}_{0} &< 0
\end{aligned}
$$

其中，

$$
\begin{aligned}
    \Psi(T_{0}, \mathrm{Gain}) &=&& \bigg(\mathcal{A} B(T_{0}, \mathrm{Gain}) + \mathcal{B} \bigg)^{\top} H \bigg(\mathcal{A} B(T_{0}, \mathrm{Gain}) + \mathcal{B} \bigg) \\
    &&&+ B(T_{0}, \mathrm{Gain})^{\top} H \bigg(\mathcal{A} B(T_{0}, \mathrm{Gain}) + \mathcal{B} \bigg) \\
    &&&+ \lambda B(T_{0}, \mathrm{Gain})^{\top} H B(T_{0}, \mathrm{Gain})
\end{aligned}
$$

根据之前的推导，可知$(A)^{N - M} C$是一个$\mathbb{R}^{n} \mapsto \mathbb{R}^{(M + 1)n}$的满射。因此，等价于找到$\mathrm{Gain}$的取值范围，使得如下LMI约束在$T_{0} \in [0,1)$恒存在的问题：

$$
\begin{aligned}
    \Psi(T_{0}, \mathrm{Gain}) &=&& \bigg(\mathcal{A} B(T_{0}, \mathrm{Gain}) + \mathcal{B} \bigg)^{\top} H \bigg(\mathcal{A} B(T_{0}, \mathrm{Gain}) + \mathcal{B} \bigg) \\
    &&&+ B(T_{0}, \mathrm{Gain})^{\top} H \bigg(\mathcal{A} B(T_{0}, \mathrm{Gain}) + \mathcal{B} \bigg) \\
    &&&+ \lambda B(T_{0}, \mathrm{Gain})^{\top} H B(T_{0}, \mathrm{Gain}) < 0
\end{aligned}
$$

至此，得到了对于数字采样反馈控制系统，满足CLF指数稳定的条件。

## 3.下一步的工作

下一步有以下打算：

1. 找个线性系统，采用相同的CLF函数指数稳定约束，看引入数字采样和时滞后增益的可行域变化；
2. Koopman算子理论线性化某一个非线性系统（实验室二楼的单臂），结合实验和仿真，辨识高维的线性系统，验证该方法在非线性系统上的准确性。

## 4.创新点

这份工作有以下创新点：

1. 在数字采样控制系统这一混杂系统的稳定性问题中，引入了**CLF函数**，并（在后续）对比了数字采样以及时滞对增益的**可行域**（暂时这么叫吧）的影响；
2. 将CLF约束在时间尺度$T = N + T_{0}(N \in \mathbb{Z}^{+}, T_{0} \in [0,1))$上恒成立的问题，简化为$\Psi$在时间尺度$T_{0} \in [0,1)$上恒成立的问题，**缩小了需要求解的时间尺度**，消除了了矩阵**幂运算**，以及一系列幂运算带来的**数值稳定性**造成的误差；
3. 使用**Koopman算子理论线性化**非线性控制系统，验证了该方法在非线性系统上的准确性。
