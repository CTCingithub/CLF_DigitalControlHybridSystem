# ����control Lyapunov function���۵����ֲ�����������ϵͳ�ȶ��������о�

## 1.����ϵͳ��ɢ��

������������ϵͳ��

$$
\dot{\vec{\theta}} = \mathcal{A} \vec{\theta} + \mathcal{B} \left. \vec{\theta} \right |_{t=(n-1)\tau},\qquad t\in[n\tau, (n+1)\tau),n \in \N,
$$

����ⳣ΢�ַ��̵õ���

$$
\begin{aligned}
    \left. \vec{\theta} \right |_{t=(n+1)\tau} &= \exp{(\mathcal{A} \tau)} \left. \vec{\theta} \right |_{t=n\tau} + \int_{0}^{\tau} \exp{(\mathcal{A} t_0)} dt_0 \mathcal{B} \left. \vec{\theta} \right |_{t=n\tau} \\
    &= \Bigg (\exp{(\mathcal{A} \tau)} + \bigg (\int_{0}^{\tau} \exp{(\mathcal{A} t_0)} dt_0 \bigg) \mathcal{B} \Bigg) \left. \vec{\theta} \right |_{t=n\tau} \\
    &= A \left. \vec{\theta} \right |_{t=n \tau},
\end{aligned}
$$

�����ַ�ʽ�õ�����ɢӳ����$\left. \overrightarrow{{}^{1}\theta} \right |_{(n+1)} = A_{1} \left. \overrightarrow{{}^{1}\theta} \right |_{n}$��ʾ�����⣬��$\left. \overrightarrow{{}^{2}\theta} \right |_{(n+1)} = A_{2} \left. \overrightarrow{{}^{2}\theta} \right |_{n}$��ʾͨ��ǰ��ŷ�����õ�����ɢӳ�䡣

### ��1

�������µ����ɶ�������ʩ��PD���ƣ�λ��������ٶ�����ֱ�Ϊ$k_p$��$k_d$��

![�����ɶ�������](images/2023-04-10-15-14-41.png)

<!--Mathematica��ͼ����
x = 0;
y = 0;
(*���ƾ���������*)
rectangle = Rectangle[{x, y}, {x + 3, y + 2}];
(*ָ�����εı߿���ɫΪ��ɫ�������ɫΪ��ɫ*)
style = EdgeForm[Black];
(*��������ͷ*)
arrow = Arrow[{{x + 3, y + 1}, {x + 6, y + 1}}];
(*��������ͷ�ϵı�ע*)
label = Text[Style["\!\(\*
StyleBox[\"F\",\nFontSlant->\"Italic\"]\)", 30], {x + 5.8, y + 1.3}];
(*�����������ϵı�ע*)
massLabel = Text[Style["\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)", 30], {x + 1.5, y + 1}];
(*���Ƶ���*)
groundlines = 
  Plot[{-1, 0}, {t, -2.5, 5.5}, Filling -> Axis, 
   FillingStyle -> LightGray, PlotStyle -> {White, Black}];
(*���ͼ��*)
Show[{Graphics[{style, FaceForm[White], rectangle, arrow, label, 
    massLabel}], groundlines}, Axes -> False, 
 AspectRatio -> Automatic]
-->

����Ķ���ѧ����Ϊ��

$$
\ddot{x} = - k_{p} \left. x \right |_{t=n \tau} - k_{d} \left. \dot{x} \right |_{t=n \tau}, \qquad t \in [n\tau, (n+1)\tau),n \in \N,
$$

��ʹ�����ODE�ķ��������壺

$$
\overrightarrow{{}^{1}\theta} = \begin{pmatrix}
    x \\ \dot{x}
\end{pmatrix}.
$$

���У�

$$
\dot{\overrightarrow{{}^{1}\theta}} = \begin{bmatrix}
    0 & 1 \\ 0 & 0
\end{bmatrix} \overrightarrow{{}^{1}\theta} + \begin{bmatrix}
    0 & 0 \\ -k_{p} & -k_{d}
\end{bmatrix} \left. \overrightarrow{{}^{1}\theta} \right |_{t=n\tau}, \qquad
t \in [n\tau, (n+1)\tau),n \in \N,
$$

��֪

$$
\mathcal{A} = \begin{bmatrix}
    0 & 1 \\ 0 & 0
\end{bmatrix}, \qquad
\mathcal{B} = \begin{bmatrix}
    0 & 0 \\ -k_{p} & -k_{d}
\end{bmatrix}.
$$

���Լ���õ���

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

������ǰ��ŷ���������壺

$$
\overrightarrow{{^{2}}\theta} = \begin{pmatrix}
    x \\ \dot{x} \\ \ddot{x}
\end{pmatrix}.
$$

�����˶�ѧ��ϵ��

$$
\begin{aligned}
    \left. x \right |_{t=(n+1)\tau} &= \left. x \right |_{t=n\tau} + \frac{\tau}{2} \left. \dot{x} \right |_{t=n\tau} + \frac{\tau^2}{2} \left. \ddot{x} \right |_{t=n\tau}, \\
    \left. \dot{x} \right |_{t=(n+1)\tau} &= \left. \dot{x} \right |_{t=n\tau} +\tau \left. \ddot{x} \right |_{t=n\tau}.
\end{aligned}
$$

�ٽ�϶���ѧ���̣��õ���

$$
A_{2} = \begin{bmatrix}
    1 & \frac{\tau}{2} & \frac{\tau^{2}}{2} \\
    0 & 1 & \frac{\tau}{2} \\
    -k_{p} & -k_{d} & 0
\end{bmatrix}.
$$

## 2.ʱ��߶ȵĸı�
