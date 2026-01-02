Full M-Set Reduction
---------------
The majority of this section is by following the math myself with a few spot checks. I don't fully understand it, so might be wrong.

This step reduces:

$$ g \rightarrow n+m $$

$$ \begin{bmatrix}
    M_{mm}       & M_{mn} \\
    M_{nm}       & M_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    \ddot u_{m} \\
    \ddot u_{n} \\
\end{Bmatrix} + 
\begin{bmatrix}
    K_{mm}       & K_{mn} \\
    K_{nm}       & K_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    u_{m} \\
    u_{n} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

Per basic dynamics, let:

  $$ u_m = [GM]{u_n} $$

Following:

$$ \begin{bmatrix}
    M_{mm}       & M_{mn} \\
    M_{nm}       & M_{nn} \\
\end{bmatrix} [GMI] \ddot u_{n} +
\begin{bmatrix}
    K_{mm}       & K_{mn} \\
    K_{nm}       & K_{nn} \\
\end{bmatrix} [GMI] u_n = \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

Per basic dynamics, let:

$$ [GMI] = \begin{bmatrix}
    [GM] & [0]  \\
    [0]  & [I]  \\
\end{bmatrix} $$

Following

$$ [GMI]^{-1} = \begin{Bmatrix}
    [GM^{-1}] & [0]  \\
    [0 ]      & [I]  \\
\end{Bmatrix} $$

$$ [GMI]^{-1} \begin{bmatrix}
    M_{mm}       & M_{mn} \\
    M_{nm}       & M_{nn} \\
\end{bmatrix} [GMI] \ddot u_{n} +
[GMI]^{-1} \begin{bmatrix}
    K_{mm}       & K_{mn} \\
    K_{nm}       & K_{nn} \\
\end{bmatrix} [GMI] u_n = [GMI]^{-1} \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

M-Set Reduction
---------------

$$ \begin{bmatrix}
    K_{mm}       & K_{mn} \\
    K_{nm}       & K_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    u_{m} \\
    u_{n} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

$$ K_{mm} u_m + K_{mn} u_n = F_m $$

$$ K_{mn} u_m + K_{nn} u_n = F_n $$

Let:

$$ u_m = G_{mn} u_n $$

$$ (K_{mm} G_{mn} + K_{mn}) u_n = F_m $$

$$ (K_{mn} G_{mn} + K_{nn}) u_n = F_n $$

MPC Definitinon
---------------
$$ \Sigma A_j u_j = 0 $$

RBE2 Definition
---------------
Independent (n) = Dependent (m)

$$ u_{m,2} = u_{n,1} $$

$$ G_mn = [1] $$

RBE3 Definition
---------------
```
       /2
#-k---1---3
       \4
```
$$ F_n1 = F_{d2} + F_{d3} + F_{d4} $$

$ K11 = 2
$ F1 = 10
K = [2 0 0 0]; F=[0]
    [0 0 0 0]    [Fd2]
    [0 0 0 0]    [Fd3]
    [0 0 0 0]    [Fd4]

$$ (K_{mn} G_{mn} + K_{nn}) u_1 = F_1 $$

$$ (K_{mn} G_{mn} + K_{nn}) u_n = F_n $$

(K_{mn} G_{mn} + K_{nn}) u_n

$$ K_{mm} u_m + K_{mn} u_n = F_m $$
$$ K_{mn} u_m + K_{nn} u_n = F_n $$

Per Elements
------------
$$ [RG_m, RG_n] [nm, xn]^T = 0 $$

$$ [RG_m] x_m + [RG_n] x_n = 0 $$

$$ x_m = -[RG_m]^{-1} [RG_n] x_n = [GM] x_n $$

So:

$$ x_m = [GM] x_n $$

$$ [GM] = -[RG_m]^{-1} [RG_n] $$

Prior to reduction:

$$ \begin{bmatrix}
    M_{mm}  & 0      \\
    0       & M_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    \ddot u_{m} \\
    \ddot u_{n} \\
\end{Bmatrix} + \begin{bmatrix}
    K_{mm}  & 0      \\
    0       & K_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    u_{m} \\
    u_{n} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

$$ M_{mm} \ddot x_m + K_{mm} x_m = F_m $$

$$ M_{nn} \ddot x_n + K_{nn} x_n = F_n $$

$$ M_{mm} [GM] \ddot x_n + K_{mm} [GM] x_n = F_m $$

Pre-multiplying:

$$ [GM]^T M_{mm} [GM] \ddot x_n + [GM]^T K_{mm} [GM] x_n = [GM]^T F_m $$

Adding:

$$ ( M_{nn} + [GM]^T M_{mm} [GM]) \ddot x_n + (K_{nn} + [GM]^T K_{mm} [GM]) x_n = F_n + [GM]^T F_m $$

S-Set Reduction
---------------
This step is straightforward, but might not be the best way to do it. It might be better to assume $u_s=0$ to simplify the math and handle $u_s != 0$ using constraints. This step reduces:

$$ n \rightarrow s+f $$

$$ \begin{bmatrix}
    K_{ss}  &  K_{sf} \\
    K_{fs}  &  K_{ff} \\
\end{bmatrix} \begin{Bmatrix}
    u_{s} \\
    u_{f} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{s} \\
    F_{f} \\
\end{Bmatrix} $$

$$ K_{ss} u_s + K_{sf} u_f = F_s $$

$$ K_{ss} u_s = F_s - K_{sf} u_f $$

$$ u_s = K_{ss}^{-1} (F_s - K_{sf} u_f) $$


$$ K_{sf} u_s + K_{ff} u_f = F_f $$

$$ K_{ff} u_f = F_f - K_{sf} u_s $$

$$ u_f = K_{ff}^{-1} (F_f - K_{sf} u_s) $$

O-Set Reduction
---------------
Per Mystran's docs (Reduction of the F-set to the A-set section), this step reduces:

$$ f \rightarrow a+o $$

$$ \begin{bmatrix}
    M_{aa}  & M_{ao}  \\
    M_{ao}^T  & M_{oo} \\
\end{bmatrix} \begin{Bmatrix}
    \ddot u_{a} \\
    \ddot u_{o} \\
\end{Bmatrix} + \begin{bmatrix}
    K_{aa}  & K_{ao} \\
    K_{ao}^T  & K_{oo} \\
\end{bmatrix} \begin{Bmatrix}
    u_{a} \\
    u_{o} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{a} \\
    F_{o} \\
\end{Bmatrix} $$

Let's do a Guyan Reduction ($\ddot u=0$):

$$ \begin{bmatrix}
    K_{aa}  & K_{ao} \\
    K_{ao}^T  & K_{oo} \\
\end{bmatrix} \begin{Bmatrix}
    u_{a} \\
    u_{o} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{a} \\
    F_{o} \\
\end{Bmatrix} $$

Take the 2nd equation and premultiply by $K_{oo}^{-1}$:

$$  K_{ao}^T u_a + K_{oo} u_o = F_o  $$
$$  K_{oo}^{-1} K_{ao}^T u_a + u_o = K_{oo}^{-1} F_o $$

Let's solve for $u_o$, so let:

$$ G_{oa} = -K_{oo}^{-1} K_{ao}^T $$
$$ u_o^0 = K_{oo}^{-1} F_o       $$

So:

$$  -G_{oa} u_a + u_o = u_o^0 $$
$$  u_o = u_o^0 + G_{oa} u_a $$

$$ \begin{Bmatrix}
    u_{a}  \\
    u_{o}  \\
\end{Bmatrix} = \begin{bmatrix}
    I_{aa}  \\
    G_{oa}  \\
\end{bmatrix} u_a + \begin{Bmatrix}
    0     \\
    u_o^0 \\
\end{Bmatrix} = [IG] u_a + \begin{Bmatrix}
    0     \\
    u_o^0 \\
\end{Bmatrix} $$

Dropping the $u_o^0$ part and pre-multiplying by [IG]^T:

$$ \begin{bmatrix}
    I_{aa}  &  G_{oa}^T  \\
\end{bmatrix} \begin{bmatrix}
    M_{aa}  & M_{ao}  \\
    M_{ao}^T  & M_{oo} \\
\end{bmatrix} \begin{bmatrix}
    I_{aa}  \\
    G_{oa}  \\
\end{bmatrix} \ddot u_a + \begin{bmatrix}
    I_{aa}  &  G_{oa}^T  \\
\end{bmatrix} \begin{bmatrix}
    K_{aa}  & K_{ao} \\
    K_{ao}^T  & K_{oo} \\
\end{bmatrix} \begin{bmatrix}
    I_{aa}  \\
    G_{oa}  \\
\end{bmatrix} u_a = \begin{bmatrix}
    I_{aa}  &  G_{oa}^T  \\
\end{bmatrix} \begin{Bmatrix}
    F_{a} \\
    F_{o} \\
\end{Bmatrix} $$

$$ \hat K_aa = \begin{bmatrix}
    I_{aa}  &  G_{oa}^T  \\
\end{bmatrix} \begin{bmatrix}
    K_{aa}    & K_{ao} \\
    K_{ao}^T  & K_{oo} \\
\end{bmatrix} \begin{bmatrix}
    I_{aa}  \\
    G_{oa}  \\
\end{bmatrix} = K_{aa} + K_{ao} G_{oa} + G_{oa}^T K_{ao}^T + G_{oa}^T K_{oo} G_{oa} $$

$$ \hat K_{aa} = K_{aa} + K_{ao} G_{oa} + (K_{ao} G_{oa})^T + G_{oa}^T K_{oo} G_{oa} $$

Plugging in $G_{oa}$ and simplifying:

$$ \hat K_{aa} = K_{aa} + K_{ao} G_{oa} $$

$$ \hat M_{aa} = M_{aa} + M_{ao} M_{oa} + (M_{ao} G_{oa})^T + G_{oa}^T M_{oo} G_{oa} $$
$$ \hat P_{a} = P_{a} + G_{oa} P_o $$

$$ [\hat M_{aa}] \ddot u_a + [\hat K_{aa}] u_a = [\hat P_a] $$
```
[I  GT] [aa  ao] [I] = [I GT] [aa*I + ao*G]  = aa+ao*G + GT*(aoT + oo*G) = aa+ao*G + GT*aoT + GT*oo*G
        [aoT oo] [G]          [aoT*I + oo*G]
```

        
A-Set Reduction
---------------
This step reduces:

$$ a \rightarrow t+q $$

q is typically 0, so a=t.

T-Set Reduction
---------------
The majority of this section is per basic dynamics (Theoretical Considerations for Using SUPORT). This step reduces:

$$ t \rightarrow l+r $$

Let:

$$ {u_l} = [D]{u_r} $$

$$ \begin{bmatrix}
    K_{ll}  &  K_{lr}  \\  
    K_{rl}  &  K_{rr}  \\
\end{bmatrix} \begin{Bmatrix}
    u_{l} \\
    u_{r} \\
\end{Bmatrix} = \begin{Bmatrix}
    0 \\
    F_{r} \\
\end{Bmatrix} $$

where

$$ [D] = [K_{ll}]^{-1} {K_{lr}} $$

$$ [\Phi_r] = \begin{Bmatrix}
    D \\
    I_{r} \\
\end{Bmatrix} $$
$$ [\Phi_r]^T = [D^T I_r] $$

$$ [\Phi_r]^T [M] [\Phi_r] > [0]  $$
$$ [\Phi_r]^T [K] [\Phi_r] = [0]  $$

Internal strain energy (work) is calculated as:

$$ [X] = [\Phi_r]^T \begin{bmatrix}
    K_{ll}  &  K_{lr}  \\  
    K_{rl}  &  K_{rr}  \\
\end{bmatrix} [\Phi_r] = [D]^T [K_{ll}] [D] + [K_{rr}] $$

Rigid body error ratio using the L2-norm is:

$$ e = \frac{| [Krr] + [K_{lr}^T][D] | } {|K_{rr}|} $$

The statics problem is:

$$ [K_{tt}] {u_t} = {F_t} $$ 

The SPC/SPCD Methods in Frequency Response Analysis
---------------------------------------------------
Redoing the sets for dynamics...

$$ ( -\omega^2 \begin{bmatrix}
    K_{ff}  &  K_{fs}  \\  
    K_{sf}  &  K_{ss}  \\
\end{bmatrix} +j \omega \begin{bmatrix}
    B_{ff}  &  B_{fs}  \\  
    B_{sf}  &  B_{ss}  \\
\end{bmatrix} + \begin{bmatrix}
    K_{ff}  &  K_{fs}  \\  
    K_{sf}  &  K_{ss}  \\
\end{bmatrix} + \begin{bmatrix}
    \hat K_{ff}  &  \hat K_{fs}  \\  
    \hat K_{sf}  &  \hat K_{ss}  \\
\end{bmatrix} ) \begin{pmatrix}
    U_f  \\  
    U_s  \\
\end{pmatrix} = \begin{Bmatrix}
    P_{f}  \\  
    P_{s}+q_s  \\
\end{Bmatrix} $$

where

$u_f$ is the displacement of the free DOFs

$u_s$ is the applied enforced motion

$q_s$ is the dynamic portion of the solution and doesn't make sense or statics.$M$ is the mass matrix

$B$ is the viscous damping matrix (BAAX?)

$K$ is the stiffness matrix (KAAX?)

$\hat K$ is the structural damping matrix (K4DD?)

$$ \hat K = GK + \sum G_E K_E $$

where

$G$   is the global/uniform structural damping

$G_E$ is the material structural damping

$K_E$ is the element stiffness matrix
