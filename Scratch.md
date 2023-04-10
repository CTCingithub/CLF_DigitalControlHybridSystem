# 基于control Lyapunov function理论的数字采样反馈控制系统稳定性问题研究

## 1.连续系统离散化

对于如下线性系统：

$$
\dot{\vec{\theta}} = \mathcal{A} \vec{\theta} + \mathcal{B} \left. \vec{\theta} \right |_{t=(n-1)\tau},\qquad t\in[n\tau, (n+1)\tau),n \in \N,
$$

可求解常微分方程得到：

$$
\begin{aligned}
    \left. \vec{\theta} \right |_{t=(n+1)\tau} &= \exp{(\mathcal{A} \tau)} \left. \vec{\theta} \right |_{t=n\tau} + \int_{0}^{\tau} \exp{(\mathcal{A} t_0)} dt_0 \mathcal{B} \left. \vec{\theta} \right |_{t=n\tau} \\
    &= \Bigg (\exp{(\mathcal{A} \tau)} + \bigg (\int_{0}^{\tau} \exp{(\mathcal{A} t_0)} dt_0 \bigg) \mathcal{B} \Bigg) \left. \vec{\theta} \right |_{t=n\tau} \\
    &= A \left. \vec{\theta} \right |_{t=n \tau},
\end{aligned}
$$

将这种方式得到的离散映射用$\left. \overrightarrow{{}^{1}\theta} \right |_{(n+1)} = A_{1} \left. \overrightarrow{{}^{1}\theta} \right |_{n}$表示。另外，用$\left. \overrightarrow{{}^{2}\theta} \right |_{(n+1)} = A_{2} \left. \overrightarrow{{}^{2}\theta} \right |_{n}$表示通过前向欧拉法得到的离散映射。

### 例1

对于如下单自由度质量块施加PD控制，位置增益和速度增益分别为$k_p$和$k_d$。

![单自由度质量块](images/2023-04-10-15-14-41.png)

<!--Mathematica作图代码
x = 0;
y = 0;
(*绘制矩形质量块*)
rectangle = Rectangle[{x, y}, {x + 3, y + 2}];
(*指定矩形的边框颜色为黑色，填充颜色为白色*)
style = EdgeForm[Black];
(*绘制力箭头*)
arrow = Arrow[{{x + 3, y + 1}, {x + 6, y + 1}}];
(*绘制力箭头上的标注*)
label = Text[Style["\!\(\*
StyleBox[\"F\",\nFontSlant->\"Italic\"]\)", 30], {x + 5.8, y + 1.3}];
(*绘制质量块上的标注*)
massLabel = Text[Style["\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)", 30], {x + 1.5, y + 1}];
(*绘制地面*)
groundlines = 
  Plot[{-1, 0}, {t, -2.5, 5.5}, Filling -> Axis, 
   FillingStyle -> LightGray, PlotStyle -> {White, Black}];
(*组合图形*)
Show[{Graphics[{style, FaceForm[White], rectangle, arrow, label, 
    massLabel}], groundlines}, Axes -> False, 
 AspectRatio -> Automatic]
-->

物体的动力学方程为：

$$
\ddot{x} = - k_{p} \left. x \right |_{t=n \tau} - k_{d} \left. \dot{x} \right |_{t=n \tau}, \qquad t \in [n\tau, (n+1)\tau),n \in \N,
$$

若使用求解ODE的方法，则定义：

$$
\overrightarrow{{}^{1}\theta} = \begin{pmatrix}
    x \\ \dot{x}
\end{pmatrix}.
$$

则有，

$$
\dot{\overrightarrow{{}^{1}\theta}} = \begin{bmatrix}
    0 & 1 \\ 0 & 0
\end{bmatrix} \overrightarrow{{}^{1}\theta} + \begin{bmatrix}
    0 & 0 \\ -k_{p} & -k_{d}
\end{bmatrix} \left. \overrightarrow{{}^{1}\theta} \right |_{t=n\tau}, \qquad
t \in [n\tau, (n+1)\tau),n \in \N,
$$

可知

$$
\mathcal{A} = \begin{bmatrix}
    0 & 1 \\ 0 & 0
\end{bmatrix}, \qquad
\mathcal{B} = \begin{bmatrix}
    0 & 0 \\ -k_{p} & -k_{d}
\end{bmatrix}.
$$

可以计算得到：

$$
\begin{aligned}
    A_{1} &= \exp{(\mathcal{A} \tau)} + \bigg (\int_{0}^{\tau} \exp{(\mathcal{A} t_0)} dt_0 \bigg) \mathcal{B} \\
    &= \begin{bmatrix}
        1 & \tau \\
        0 & 1
    \end{bmatrix} + \begin{bmatrix}
        \tau & \frac{\tau^{2}}{2} \\
        0 & \tau
    \end{bmatrix} \begin{bmatrix}
        0 & 0 \\ -k_{p} & -k_{d}
    \end{bmatrix} \\
    &= \begin{bmatrix}
        1 - k_{p} \frac{\tau^{2}}{2} & \tau - k_{d} \frac{\tau^{2}}{2} \\
        -k_{p} \tau & 1 - k_{d} \tau
    \end{bmatrix}.
\end{aligned}
$$

若采用前向欧拉法，则定义：

$$
\overrightarrow{{^{2}}\theta} = \begin{pmatrix}
    x \\ \dot{x} \\ \ddot{x}
\end{pmatrix}.
$$

根据运动学关系：

$$
\begin{aligned}
    \left. x \right |_{t=(n+1)\tau} &= \left. x \right |_{t=n\tau} + \frac{\tau}{2} \left. \dot{x} \right |_{t=n\tau} + \frac{\tau^2}{2} \left. \ddot{x} \right |_{t=n\tau}, \\
    \left. \dot{x} \right |_{t=(n+1)\tau} &= \left. \dot{x} \right |_{t=n\tau} +\tau \left. \ddot{x} \right |_{t=n\tau}.
\end{aligned}
$$

再结合动力学方程，得到：

$$
A_{2} = \begin{bmatrix}
    1 & \frac{\tau}{2} & \frac{\tau^{2}}{2} \\
    0 & 1 & \frac{\tau}{2} \\
    -k_{p} & -k_{d} & 0
\end{bmatrix}.
$$

## 2.时间尺度的改变
