# Derivation of the Schwarzschild Equations of Motion (Hamiltonian Formalism)

This document derives the equations of motion for geodesics of test particles and photons in a static, spherically symmetric spacetime using Hamiltonian mechanics. It ends with a system of first-order ordinary differential equations that can be integrated numerically with the fourth-order Runge–Kutta method. The implementation is in the [`schwarzschild`](../src/schwarzschild) Python package of this repository; Section 16.5 maps the equations to the code.

A PDF version of this document (built from [`derivation.tex`](derivation.tex)) is attached to each release.

---

## 1. The Covariant Metric ($g_{\mu\nu}$)

For a static spacetime with spherical symmetry the line element $ds^2$ is

$$
ds^2 = -f(r)\,dt^2 + \frac{dr^2}{f(r)} + r^2 \left( d\theta^2 + \sin^2\theta\,d\phi^2 \right) \tag{1}
$$

From equation (1) the covariant metric $g_{\mu\nu}$ is diagonal with non-zero components

$$
\begin{align}
    g_{tt} &= -f(r) \tag{2a} \\
    g_{rr} &= \frac{1}{f(r)} \tag{2b} \\
    g_{\theta\theta} &= r^2 \tag{2c} \\
    g_{\phi\phi} &= r^2 \sin^2\theta \tag{2d}
\end{align}
$$

---

## 2. The Contravariant Metric ($g^{\mu\nu}$)

The contravariant metric $g^{\mu\nu}$ is the inverse of $g_{\mu\nu}$. Since $g_{\mu\nu}$ is diagonal, each contravariant component is simply the reciprocal of the corresponding covariant component in equation (2):

$$
\begin{align}
    g^{tt} &= \frac{1}{g_{tt}} = -\frac{1}{f(r)} \tag{3a} \\
    g^{rr} &= \frac{1}{g_{rr}} = f(r) \tag{3b} \\
    g^{\theta\theta} &= \frac{1}{g_{\theta\theta}} = \frac{1}{r^2} \tag{3c} \\
    g^{\phi\phi} &= \frac{1}{g_{\phi\phi}} = \frac{1}{r^2 \sin^2\theta} \tag{3d}
\end{align}
$$

---

## 3. The Hamiltonian

In General Relativity the Hamiltonian of a particle (or photon) is

$$
H = \frac{1}{2} g^{\mu\nu} p_\mu p_\nu \tag{4}
$$

Substituting the contravariant components from equation (3) into equation (4) gives

$$
H = \frac{1}{2} \left[ -\frac{p_t^2}{f(r)} + f(r) p_r^2 + \frac{p_\theta^2}{r^2} + \frac{p_\phi^2}{r^2 \sin^2\theta} \right] \tag{5}
$$

---

## 4. Hamilton's Equations and Conserved Quantities

The momenta conjugate to the coordinates obey the canonical equation, with $\lambda$ an affine parameter:

$$
\dot{p}_\mu = -\frac{\partial H}{\partial x^\mu} \tag{6}
$$

The Hamiltonian (5) does not depend explicitly on the time coordinate $t$ or on the azimuthal angle $\phi$ (cyclic coordinates). Consequently:

* **Conservation of energy ($E$):**
  $$
  \frac{\partial H}{\partial t} = 0 \implies \dot{p}_t = 0 \implies p_t = -E \quad (\text{constant}) \tag{7a}
  $$
  This constant of motion is defined as the negative of the specific energy of the particle.

* **Conservation of angular momentum ($L$):**
  $$
  \frac{\partial H}{\partial \phi} = 0 \implies \dot{p}_\phi = 0 \implies p_\phi = L \quad (\text{constant}) \tag{7b}
  $$
  This constant of motion is the specific angular momentum about the symmetry axis.

---

## 5. The Equatorial Plane

Because the spacetime is fully spherically symmetric, a particle's trajectory always stays in a single plane. Without loss of generality the coordinates can be oriented so that the motion takes place in the equatorial plane:

$$
\theta = \frac{\pi}{2}, \quad p_\theta = 0, \quad \sin\theta = 1
$$

Substituting this condition and the conserved quantities (7a) and (7b) into equation (5), the Hamiltonian reduces to

$$
\begin{align}
    H &= \frac{1}{2} \left[ -\frac{(-E)^2}{f(r)} + f(r) p_r^2 + 0 + \frac{L^2}{r^2(1)} \right] \nonumber \\
    H &= \frac{1}{2} \left[ -\frac{E^2}{f(r)} + f(r) p_r^2 + \frac{L^2}{r^2} \right] \tag{8}
\end{align}
$$

---

## 6. Momentum and Radial Velocity

The second canonical equation relates the coordinate velocity to the derivative of the Hamiltonian with respect to the conjugate momentum:

$$
\dot{x}^\mu = \frac{\partial H}{\partial p_\mu} \tag{9}
$$

For the radial component,

$$
\begin{align}
    \dot{r} &= \frac{\partial H}{\partial p_r} = \frac{1}{2} \left[ 2 f(r) p_r \right] \nonumber \\
    \dot{r} &= f(r) p_r \tag{10}
\end{align}
$$

The radial momentum can therefore be written in terms of the radial velocity:

$$
p_r = \frac{\dot{r}}{f(r)} \tag{11}
$$

---

## 7. The Radial Equation of Motion

From the normalization of the 4-momentum ($g^{\mu\nu} p_\mu p_\nu = -\epsilon$), the constant value of the Hamiltonian is

$$
H = -\frac{\epsilon}{2} \tag{12}
$$

where

* $\epsilon = 1$ for massive particles (timelike geodesics),
* $\epsilon = 0$ for photons and other massless particles (null geodesics).

Substitute equations (11) and (12) into equation (8):

$$
\begin{align}
    -\frac{\epsilon}{2} &= \frac{1}{2} \left[ -\frac{E^2}{f(r)} + f(r) \left( \frac{\dot{r}}{f(r)} \right)^2 + \frac{L^2}{r^2} \right] \nonumber \\
    -\epsilon &= -\frac{E^2}{f(r)} + f(r) \frac{\dot{r}^2}{f(r)^2} + \frac{L^2}{r^2} \nonumber \\
    -\epsilon &= -\frac{E^2}{f(r)} + \frac{\dot{r}^2}{f(r)} + \frac{L^2}{r^2}
\end{align}
$$

Multiply both sides by $f(r)$:

$$
-\epsilon f(r) = -E^2 + \dot{r}^2 + f(r) \frac{L^2}{r^2}
$$

Move terms to isolate $\dot{r}^2$:

$$
\begin{align}
    \dot{r}^2 &= E^2 - \epsilon f(r) - f(r) \frac{L^2}{r^2} \nonumber \\
    \dot{r}^2 &= E^2 - f(r) \left( \epsilon + \frac{L^2}{r^2} \right) \tag{13}
\end{align}
$$

---

## 8. The Effective Potential

The radial equation (13) can be written like the conservation of mechanical energy in classical mechanics:

$$
\dot{r}^2 + V_{\text{eff}}(r) = E^2 \tag{14}
$$

by defining the **effective potential** $V_{\text{eff}}(r)$ as

$$
V_{\text{eff}}(r) = f(r) \left( \epsilon + \frac{L^2}{r^2} \right) \tag{15}
$$

### Special Case: the Schwarzschild Metric

Outside a non-rotating body of mass $M$ (the Schwarzschild metric) the metric function is

$$
f(r) = 1 - \frac{2M}{r} \tag{16}
$$

in geometrized units where $G = c = 1$.

Substituting $f(r)$ into equation (15) gives the explicit effective potential

$$
\begin{align}
    V_{\text{eff}}(r) &= \left( 1 - \frac{2M}{r} \right) \left( \epsilon + \frac{L^2}{r^2} \right) \nonumber \\
    V_{\text{eff}}(r) &= \epsilon - \frac{2M\epsilon}{r} + \frac{L^2}{r^2} - \frac{2ML^2}{r^3} \tag{17}
\end{align}
$$

Physical interpretation of the terms:

1. $\epsilon$: the rest energy of the particle ($1$ for massive particles, $0$ for photons).
2. $-\frac{2M\epsilon}{r}$: the Newtonian gravitational potential.
3. $\frac{L^2}{r^2}$: the Newtonian centrifugal barrier.
4. $-\frac{2ML^2}{r^3}$: the General Relativity correction, responsible for the precession of the orbit's periapsis and for the strong attraction at short distances.

---

## 9. Why Use the Full Hamilton Equations?

The radial equation (13) has the form $\dot{r}^2 = E^2 - V_{\text{eff}}(r)$. For numerical work this form is **not ideal**:

1. Obtaining $\dot{r}$ requires a square root, $\dot{r} = \pm\sqrt{E^2 - V_{\text{eff}}}$, so the $\pm$ sign has to be switched by hand each time the particle passes a turning point ($\dot{r} = 0$), for example at periapsis or apoapsis.
2. Near a turning point $E^2 - V_{\text{eff}} \to 0$; a small numerical error can make the argument of the square root negative and the integration fails.

The solution is to use **all** of Hamilton's canonical equations (6) and (9). They form a system of **first-order** ordinary differential equations (ODEs) that has no square root and passes through turning points automatically:

$$
\boxed{\;\dot{x}^\mu = \frac{\partial H}{\partial p_\mu}, \qquad \dot{p}_\mu = -\frac{\partial H}{\partial x^\mu}\;} \tag{18}
$$

where a dot denotes the derivative with respect to the affine parameter $\lambda$ (for massive particles $\lambda = \tau$ is the proper time). The 4 coordinates $(t, r, \theta, \phi)$ and 4 momenta $(p_t, p_r, p_\theta, p_\phi)$ give 8 first-order ODEs.

As a starting point, the general Hamiltonian (5) is rewritten:

$$
H = \frac{1}{2} \left[ -\frac{p_t^2}{f(r)} + f(r)\, p_r^2 + \frac{p_\theta^2}{r^2} + \frac{p_\phi^2}{r^2 \sin^2\theta} \right]
$$

---

## 10. Velocity Equations ($\dot{x}^\mu = \partial H / \partial p_\mu$)

Each momentum $p_\mu$ appears in exactly one term inside the brackets of $H$, so its partial derivative follows immediately.

**The $t$ component.** Only the first term depends on $p_t$:

$$
\dot{t} = \frac{\partial H}{\partial p_t} = \frac{1}{2}\left( -\frac{2 p_t}{f(r)} \right) = -\frac{p_t}{f(r)} \tag{19a}
$$

**The $r$ component.** Only the second term depends on $p_r$ (same as equation (10)):

$$
\dot{r} = \frac{\partial H}{\partial p_r} = \frac{1}{2}\left( 2 f(r)\, p_r \right) = f(r)\, p_r \tag{19b}
$$

**The $\theta$ component.** Only the third term depends on $p_\theta$:

$$
\dot{\theta} = \frac{\partial H}{\partial p_\theta} = \frac{1}{2}\left( \frac{2 p_\theta}{r^2} \right) = \frac{p_\theta}{r^2} \tag{19c}
$$

**The $\phi$ component.** Only the fourth term depends on $p_\phi$:

$$
\dot{\phi} = \frac{\partial H}{\partial p_\phi} = \frac{1}{2}\left( \frac{2 p_\phi}{r^2 \sin^2\theta} \right) = \frac{p_\phi}{r^2 \sin^2\theta} \tag{19d}
$$

Note that equation (19) is just $\dot{x}^\mu = g^{\mu\nu} p_\nu = p^\mu$, the relation between covariant and contravariant momenta (raising an index with the metric).

---

## 11. Momentum Equations ($\dot{p}_\mu = -\partial H / \partial x^\mu$)

### 11.1 The $t$ and $\phi$ components (cyclic coordinates)

Because $H$ does not contain $t$ or $\phi$ explicitly,

$$
\begin{align}
    \dot{p}_t &= -\frac{\partial H}{\partial t} = 0 \tag{20a} \\
    \dot{p}_\phi &= -\frac{\partial H}{\partial \phi} = 0 \tag{20b}
\end{align}
$$

which recovers the conserved quantities $p_t = -E$ and $p_\phi = L$ of equation (7).

### 11.2 The $r$ component

The coordinate $r$ appears in all four terms. Write $f'(r) \equiv \dfrac{df}{dr}$ and differentiate term by term:

* First term (use $\frac{d}{dr}\left(\frac{1}{f}\right) = -\frac{f'}{f^2}$):
  $$
  \frac{\partial}{\partial r}\left( -\frac{p_t^2}{f} \right) = -p_t^2 \left( -\frac{f'}{f^2} \right) = \frac{f'\, p_t^2}{f^2}
  $$
* Second term:
  $$
  \frac{\partial}{\partial r}\left( f\, p_r^2 \right) = f'\, p_r^2
  $$
* Third term:
  $$
  \frac{\partial}{\partial r}\left( \frac{p_\theta^2}{r^2} \right) = -\frac{2 p_\theta^2}{r^3}
  $$
* Fourth term:
  $$
  \frac{\partial}{\partial r}\left( \frac{p_\phi^2}{r^2 \sin^2\theta} \right) = -\frac{2 p_\phi^2}{r^3 \sin^2\theta}
  $$

Add them and multiply by $\tfrac{1}{2}$:

$$
\frac{\partial H}{\partial r} = \frac{f'}{2}\left( \frac{p_t^2}{f^2} + p_r^2 \right) - \frac{p_\theta^2}{r^3} - \frac{p_\phi^2}{r^3 \sin^2\theta}
$$

Hence

$$
\dot{p}_r = -\frac{\partial H}{\partial r} = -\frac{f'(r)}{2}\left( \frac{p_t^2}{f(r)^2} + p_r^2 \right) + \frac{p_\theta^2}{r^3} + \frac{p_\phi^2}{r^3 \sin^2\theta} \tag{20c}
$$

### 11.3 The $\theta$ component

Only the fourth term depends on $\theta$. With $\frac{d}{d\theta}\left(\sin^{-2}\theta\right) = -2\sin^{-3}\theta \cos\theta$,

$$
\frac{\partial H}{\partial \theta} = \frac{1}{2}\, \frac{p_\phi^2}{r^2} \left( -\frac{2\cos\theta}{\sin^3\theta} \right) = -\frac{p_\phi^2 \cos\theta}{r^2 \sin^3\theta}
$$

Hence

$$
\dot{p}_\theta = -\frac{\partial H}{\partial \theta} = \frac{p_\phi^2 \cos\theta}{r^2 \sin^3\theta} \tag{20d}
$$

### 11.4 Explicit form for Schwarzschild

For $f(r) = 1 - \dfrac{2M}{r}$ (equation (16)) the derivative is

$$
f'(r) = \frac{2M}{r^2} \tag{21}
$$

and the combination appearing in (20c) simplifies:

$$
\frac{f'}{f^2} = \frac{2M/r^2}{\left(1 - \frac{2M}{r}\right)^2} = \frac{2M/r^2}{(r - 2M)^2 / r^2} = \frac{2M}{(r - 2M)^2}
$$

Substituting into (20c) gives

$$
\dot{p}_r = -\frac{M\, p_t^2}{(r - 2M)^2} - \frac{M\, p_r^2}{r^2} + \frac{p_\theta^2}{r^3} + \frac{p_\phi^2}{r^3 \sin^2\theta} \tag{22}
$$

### 11.5 Summary: the full 8-equation system (3D)

$$
\begin{align}
    \dot{t} &= -\frac{p_t}{1 - 2M/r} &
    \dot{p}_t &= 0 \nonumber \\
    \dot{r} &= \left(1 - \frac{2M}{r}\right) p_r &
    \dot{p}_r &= -\frac{M\, p_t^2}{(r - 2M)^2} - \frac{M\, p_r^2}{r^2} + \frac{p_\theta^2}{r^3} + \frac{p_\phi^2}{r^3 \sin^2\theta} \nonumber \\
    \dot{\theta} &= \frac{p_\theta}{r^2} &
    \dot{p}_\theta &= \frac{p_\phi^2 \cos\theta}{r^2 \sin^3\theta} \nonumber \\
    \dot{\phi} &= \frac{p_\phi}{r^2 \sin^2\theta} &
    \dot{p}_\phi &= 0 \tag{23}
\end{align}
$$

System (23) holds for motion in an arbitrary plane and can be used directly to simulate 3D motion (for example many photons with different orientations for ray tracing). It is singular at $\sin\theta = 0$ (the polar axis), however, so for a single particle it is more practical to use the reduction to the equatorial plane below.

---

## 12. Reduction to the Equatorial Plane: the System Used in Simulations

### 12.1 Consistency of the equatorial plane

Choose the initial conditions $\theta(0) = \pi/2$ and $p_\theta(0) = 0$. From (23),

$$
\dot{\theta}\big|_{\lambda=0} = \frac{0}{r^2} = 0, \qquad
\dot{p}_\theta\big|_{\lambda=0} = \frac{L^2 \cos(\pi/2)}{r^2 \sin^3(\pi/2)} = \frac{L^2 \cdot 0}{r^2} = 0
$$

Since neither $\theta$ nor $p_\theta$ changes, they **remain** $\pi/2$ and $0$ along the whole trajectory. The motion is indeed confined to the equatorial plane, which justifies the assumption of Section 5.

### 12.2 The reduced system

Substitute $\theta = \pi/2$, $p_\theta = 0$, $p_t = -E$ and $p_\phi = L$ into (23). The equations for $\theta$, $p_\theta$, $p_t$ and $p_\phi$ become trivial, leaving **four first-order ODEs**:

$$
\boxed{
\begin{aligned}
    \frac{dt}{d\lambda} &= \frac{E}{1 - \dfrac{2M}{r}} \\[6pt]
    \frac{dr}{d\lambda} &= \left(1 - \frac{2M}{r}\right) p_r \\[6pt]
    \frac{d\phi}{d\lambda} &= \frac{L}{r^2} \\[6pt]
    \frac{dp_r}{d\lambda} &= -\frac{M E^2}{(r - 2M)^2} - \frac{M\, p_r^2}{r^2} + \frac{L^2}{r^3}
\end{aligned}
} \tag{24}
$$

Equation (24) is the **final result, ready for numerical simulation**. Some important remarks:

* **The parameter $\epsilon$ does not appear explicitly in (24).** Massive particles ($\epsilon = 1$) and photons ($\epsilon = 0$) obey exactly the same equations; the difference enters only through the **initial condition** $p_r(0)$ (see Section 14).
* The right-hand side of (24) **does not depend** on $t$, $\phi$ or $\lambda$ (an autonomous system). Only $r$ and $p_r$ determine the rates of change; $t$ and $\phi$ are simply "accumulated".
* Physical meaning of the terms in $\dot{p}_r$:
  * $-\dfrac{ME^2}{(r-2M)^2}$: the gravitational pull, affected by the energy (it diverges at the horizon),
  * $-\dfrac{M p_r^2}{r^2}$: a correction from the curvature of the radial part of the metric,
  * $+\dfrac{L^2}{r^3}$: the centrifugal force.

### 12.3 General form (arbitrary $f(r)$)

To change the metric (for example Reissner–Nordström, de Sitter–Schwarzschild, etc.), use the general form before substituting $f$:

$$
\begin{aligned}
    \dot{t} &= \frac{E}{f(r)}, &
    \dot{r} &= f(r)\, p_r, &
    \dot{\phi} &= \frac{L}{r^2}, &
    \dot{p}_r &= -\frac{f'(r)}{2}\left( \frac{E^2}{f(r)^2} + p_r^2 \right) + \frac{L^2}{r^3}
\end{aligned} \tag{25}
$$

so the code only needs different functions `f(r)` and `df(r)` (in the package: a class implementing the `Metric` protocol of [`metric.py`](../src/schwarzschild/metric.py)).

---

## 13. Verification: Consistency with the Effective Potential

To check (24), derive the radial acceleration $\ddot{r}$ from Hamilton's equations and compare it with the result from equation (14).

**From equation (14).** Differentiate $\dot{r}^2 + V_{\text{eff}} = E^2$ with respect to $\lambda$:

$$
2\dot{r}\ddot{r} + \frac{dV_{\text{eff}}}{dr}\dot{r} = 0 \implies \ddot{r} = -\frac{1}{2}\frac{dV_{\text{eff}}}{dr}
$$

With $V_{\text{eff}} = f\left(\epsilon + \frac{L^2}{r^2}\right)$,

$$
\ddot{r} = -\frac{f'}{2}\left(\epsilon + \frac{L^2}{r^2}\right) + \frac{f L^2}{r^3} \tag{26}
$$

**From Hamilton's equations.** Differentiate (19b) with the chain rule, $\frac{d f}{d\lambda} = f'\dot{r} = f' f p_r$:

$$
\ddot{r} = \frac{d}{d\lambda}\left( f p_r \right) = f' \dot{r}\, p_r + f \dot{p}_r = f f' p_r^2 + f \left[ -\frac{f'}{2}\left( \frac{E^2}{f^2} + p_r^2 \right) + \frac{L^2}{r^3} \right]
$$

$$
\ddot{r} = \frac{f'}{2}\left( f p_r^2 \right) - \frac{f'}{2}\frac{E^2}{f} + \frac{f L^2}{r^3}
$$

The Hamiltonian constraint $H = -\epsilon/2$ (equations (8) and (12)) gives $f p_r^2 = \dfrac{E^2}{f} - \epsilon - \dfrac{L^2}{r^2}$. Substituting,

$$
\ddot{r} = \frac{f'}{2}\left( \frac{E^2}{f} - \epsilon - \frac{L^2}{r^2} \right) - \frac{f'}{2}\frac{E^2}{f} + \frac{f L^2}{r^3}
= -\frac{f'}{2}\left(\epsilon + \frac{L^2}{r^2}\right) + \frac{f L^2}{r^3}
$$

This result is **identical** to (26). ✓

For Schwarzschild ($f' = 2M/r^2$),

$$
\ddot{r} = -\frac{M\epsilon}{r^2} + \frac{L^2}{r^3} - \frac{3ML^2}{r^4} \tag{27}
$$

which is the Newtonian "force" plus the relativistic correction $-3ML^2/r^4$.

---

## 14. Initial Conditions, Control Quantities and Stopping Criteria

### 14.1 Units

Use geometrized units $G = c = 1$ and choose $M = 1$. Then $r$, $t$, $\lambda$ and $L$ are expressed in units of $M$, while $E$ is dimensionless. The event horizon is at $r_s = 2M = 2$.

### 14.2 Determining $p_r(0)$ from the Hamiltonian constraint

The initial state vector is $\left(t_0, r_0, \phi_0, p_{r,0}\right)$, usually with $t_0 = 0$ and $\phi_0 = 0$. The constants $E$ and $L$ can be chosen freely, but $p_{r,0}$ is **not** free: it must satisfy $H = -\epsilon/2$. From (13) and (11),

$$
p_{r,0} = \sigma\, \frac{\sqrt{E^2 - V_{\text{eff}}(r_0)}}{f(r_0)}, \qquad \sigma = \begin{cases} -1 & \text{ingoing} \\ +1 & \text{outgoing} \end{cases} \tag{28}
$$

provided $E^2 \geq V_{\text{eff}}(r_0)$. The square root is used **only once**, here; during the integration the sign of $p_r$ is handled automatically by (24).

### 14.3 Standard test cases

**(a) Circular orbit of a massive particle ($\epsilon = 1$).** The conditions are $V_{\text{eff}}'(r_c) = 0$ and $E^2 = V_{\text{eff}}(r_c)$. From $V_{\text{eff}}' = \frac{2M}{r^2}\left(1 + \frac{L^2}{r^2}\right) - \frac{2 f L^2}{r^3} = 0$, multiplying by $r^4/2$:

$$
M r^2 + M L^2 - L^2 r + 2ML^2 = 0 \implies L^2 (r - 3M) = M r^2
$$

$$
L = \sqrt{\frac{M r_c^2}{r_c - 3M}}, \qquad
E = \frac{1 - 2M/r_c}{\sqrt{1 - 3M/r_c}}, \qquad
p_{r,0} = 0 \tag{29}
$$

Circular orbits are stable only for $r_c \geq 6M$ (the *Innermost Stable Circular Orbit*, ISCO). The numerical result must give a constant $r(\lambda) = r_c$ and the angular frequency $d\phi/dt = \sqrt{M/r_c^3}$ (Kepler's law, which is exact in Schwarzschild coordinates).

**(b) Bound orbit with precession ($\epsilon = 1$).** Start at apoapsis: choose $r_0$ and $L$, set $p_{r,0} = 0$ and $E = \sqrt{V_{\text{eff}}(r_0)}$. For example $r_0 = 20M$, $L = 4.2M$ gives a "rosette" orbit with $r \in [10.33M,\; 20M]$.

**(c) Photons ($\epsilon = 0$).** Photon paths depend only on the impact parameter $b = L/E$, so $E = 1$ and $L = b$ can be used. The peak of the potential $V_{\text{eff}} = fL^2/r^2$ is at $r = 3M$ (the *photon sphere*) with value $L^2/(27M^2)$, so the critical impact parameter is

$$
b_c = 3\sqrt{3}\, M \approx 5.196\, M \tag{30}
$$

Photons with $b < b_c$ are captured by the black hole; those with $b > b_c$ are deflected and escape. For a photon coming from far away ($r_0 \gg M$):
$p_{r,0} = -\sqrt{1 - f(r_0) b^2 / r_0^2}\,/\,f(r_0)$.

### 14.4 The Hamiltonian constraint as an error check

Because $H$ is conserved exactly, the quantity

$$
\mathcal{C}(\lambda) = H + \frac{\epsilon}{2} = \frac{1}{2}\left[ -\frac{E^2}{f(r)} + f(r)\, p_r^2 + \frac{L^2}{r^2} \right] + \frac{\epsilon}{2} \tag{31}
$$

should be zero along the whole trajectory. In a simulation $|\mathcal{C}|$ is computed at every step as an accuracy indicator. Close to the horizon the term $E^2/f$ grows, so it is better to monitor the relative error $|\mathcal{C}| \,/\, \left(E^2 / 2f\right)$.

### 14.5 Stopping criteria

The integration stops when any of the following holds:

* $r \leq 2M(1 + \delta)$, with $\delta \sim 10^{-3}$: the particle falls into the horizon. Schwarzschild coordinates are singular at $r = 2M$ ($f \to 0$, so $\dot{t} \to \infty$ and $p_r \to \infty$), so the integration cannot continue through the horizon.
* $r \geq r_{\text{max}}$: the particle or photon has escaped to infinity.
* $\lambda \geq \lambda_{\text{max}}$ or the maximum number of steps is reached.

### 14.6 Conversion to Cartesian coordinates (visualization)

In the equatorial plane,

$$
x = r\cos\phi, \qquad y = r\sin\phi \tag{32}
$$

---

## 15. The Fourth-Order Runge–Kutta Method (RK4)

### 15.1 The initial value problem

RK4 solves a system of first-order ODEs of the form

$$
\frac{d\mathbf{y}}{d\lambda} = \mathbf{F}(\lambda, \mathbf{y}), \qquad \mathbf{y}(\lambda_0) = \mathbf{y}_0 \tag{33}
$$

where $\mathbf{y}$ is the state vector and $\mathbf{F}$ the right-hand side. The parameter $\lambda$ is discretized with step $h$: $\lambda_n = \lambda_0 + n h$ and $\mathbf{y}_n \approx \mathbf{y}(\lambda_n)$.

### 15.2 The basic idea

Exactly, $\mathbf{y}_{n+1} = \mathbf{y}_n + \int_{\lambda_n}^{\lambda_n + h} \mathbf{F}(\lambda, \mathbf{y}(\lambda))\, d\lambda$. Euler's method approximates this integral using only the slope at the start of the interval, $h\,\mathbf{F}(\lambda_n, \mathbf{y}_n)$, so its error is $\mathcal{O}(h)$. RK4 improves on this by evaluating the slope at **four points** in the interval (start, midpoint twice, end) and taking a weighted average, analogous to Simpson's rule $\frac{h}{6}\left[F_{\text{start}} + 4F_{\text{mid}} + F_{\text{end}}\right]$, where the midpoint weight $4$ is split into $2 + 2$.

### 15.3 The RK4 formulas

$$
\begin{align}
    \mathbf{k}_1 &= \mathbf{F}\left(\lambda_n,\; \mathbf{y}_n\right) \tag{34a} \\
    \mathbf{k}_2 &= \mathbf{F}\left(\lambda_n + \tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_1\right) \tag{34b} \\
    \mathbf{k}_3 &= \mathbf{F}\left(\lambda_n + \tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_2\right) \tag{34c} \\
    \mathbf{k}_4 &= \mathbf{F}\left(\lambda_n + h,\; \mathbf{y}_n + h\,\mathbf{k}_3\right) \tag{34d} \\[4pt]
    \mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{6}\left( \mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4 \right) \tag{34e}
\end{align}
$$

Meaning of each stage:

* $\mathbf{k}_1$: the slope at the start of the interval (same as Euler),
* $\mathbf{k}_2$: the slope at the midpoint, using a half-step Euler prediction with $\mathbf{k}_1$,
* $\mathbf{k}_3$: the slope at the midpoint again, using the improved prediction with $\mathbf{k}_2$,
* $\mathbf{k}_4$: the slope at the end of the interval, using a full-step prediction with $\mathbf{k}_3$.

### 15.4 Error

The local (per-step) error of RK4 is $\mathcal{O}(h^5)$ and the global error (accumulated over $N \propto 1/h$ steps) is $\mathcal{O}(h^4)$. Halving $h$ therefore reduces the global error by a factor of about $2^4 = 16$. This gives a simple convergence test: run the simulation with $h$ and $h/2$ and compare.

Note: RK4 is **not** a symplectic integrator, so $H$ drifts slowly over very long integrations (thousands of orbits). For simulations of tens to hundreds of orbits with a reasonable $h$ the drift is very small and can be monitored through $\mathcal{C}(\lambda)$ in (31).

---

## 16. Applying RK4 to the Schwarzschild Equations of Motion

### 16.1 State vector and right-hand side

From system (24), define

$$
\mathbf{y} = \begin{pmatrix} t \\ r \\ \phi \\ p_r \end{pmatrix}, \qquad
\mathbf{F}(\mathbf{y}) = \begin{pmatrix}
    F^t \\ F^r \\ F^\phi \\ F^{p}
\end{pmatrix} = \begin{pmatrix}
    \dfrac{E}{1 - 2M/r} \\[8pt]
    \left(1 - \dfrac{2M}{r}\right) p_r \\[8pt]
    \dfrac{L}{r^2} \\[8pt]
    -\dfrac{ME^2}{(r-2M)^2} - \dfrac{M p_r^2}{r^2} + \dfrac{L^2}{r^3}
\end{pmatrix} \tag{35}
$$

with $E$, $L$ and $M$ constant parameters. Since $\mathbf{F}$ does not depend on $\lambda$, the argument $\lambda$ in (34) can be ignored.

### 16.2 One RK4 step in explicit form

Let the state at step $n$ be $(t_n, r_n, \phi_n, p_{r,n})$. Since $\mathbf{F}$ depends only on $r$ and $p_r$, the only intermediate values that need to be computed are $r$ and $p_r$.

**Stage 1** (at $r_n,\, p_{r,n}$):

$$
\begin{aligned}
    k_1^t &= \frac{E}{f(r_n)}, &
    k_1^r &= f(r_n)\, p_{r,n}, &
    k_1^\phi &= \frac{L}{r_n^2}, &
    k_1^p &= G(r_n, p_{r,n})
\end{aligned}
$$

with $f(r) = 1 - 2M/r$ and $G(r, p_r) \equiv -\dfrac{ME^2}{(r-2M)^2} - \dfrac{M p_r^2}{r^2} + \dfrac{L^2}{r^3}$.

**Stage 2** (at the midpoint using $\mathbf{k}_1$):

$$
r^{(2)} = r_n + \tfrac{h}{2}k_1^r, \qquad p^{(2)} = p_{r,n} + \tfrac{h}{2}k_1^p
$$

$$
\begin{aligned}
    k_2^t &= \frac{E}{f(r^{(2)})}, &
    k_2^r &= f(r^{(2)})\, p^{(2)}, &
    k_2^\phi &= \frac{L}{\left(r^{(2)}\right)^2}, &
    k_2^p &= G\left(r^{(2)}, p^{(2)}\right)
\end{aligned}
$$

**Stage 3** (at the midpoint using $\mathbf{k}_2$):

$$
r^{(3)} = r_n + \tfrac{h}{2}k_2^r, \qquad p^{(3)} = p_{r,n} + \tfrac{h}{2}k_2^p
$$

$$
\begin{aligned}
    k_3^t &= \frac{E}{f(r^{(3)})}, &
    k_3^r &= f(r^{(3)})\, p^{(3)}, &
    k_3^\phi &= \frac{L}{\left(r^{(3)}\right)^2}, &
    k_3^p &= G\left(r^{(3)}, p^{(3)}\right)
\end{aligned}
$$

**Stage 4** (at the end of the interval using $\mathbf{k}_3$):

$$
r^{(4)} = r_n + h\,k_3^r, \qquad p^{(4)} = p_{r,n} + h\,k_3^p
$$

$$
\begin{aligned}
    k_4^t &= \frac{E}{f(r^{(4)})}, &
    k_4^r &= f(r^{(4)})\, p^{(4)}, &
    k_4^\phi &= \frac{L}{\left(r^{(4)}\right)^2}, &
    k_4^p &= G\left(r^{(4)}, p^{(4)}\right)
\end{aligned}
$$

**Update:**

$$
\begin{align}
    t_{n+1} &= t_n + \frac{h}{6}\left( k_1^t + 2k_2^t + 2k_3^t + k_4^t \right) \nonumber \\
    r_{n+1} &= r_n + \frac{h}{6}\left( k_1^r + 2k_2^r + 2k_3^r + k_4^r \right) \nonumber \\
    \phi_{n+1} &= \phi_n + \frac{h}{6}\left( k_1^\phi + 2k_2^\phi + 2k_3^\phi + k_4^\phi \right) \nonumber \\
    p_{r,n+1} &= p_{r,n} + \frac{h}{6}\left( k_1^p + 2k_2^p + 2k_3^p + k_4^p \right) \tag{36}
\end{align}
$$

### 16.3 Choosing the step size $h$

* Far from the black hole ($r \gtrsim 6M$) a fixed step $h \sim 0.01$–$0.1\,M$ is already very accurate.
* Close to the horizon $p_r \sim 1/(r - 2M)$ and $\dot{p}_r \sim 1/(r - 2M)^2$ grow sharply. With a fixed $h$ the RK4 stages can "jump" to $r < 2M$ and produce meaningless values. A simple remedy is to shrink the step in proportion to the distance from the horizon:

$$
h_n = \min\left( h_0,\; \alpha\,(r_n - 2M) \right), \qquad \alpha \sim 0.02 \tag{37}
$$

  Since $\mathbf{F}$ does not depend on $\lambda$, changing $h$ from step to step does not break the RK4 scheme (just accumulate $\lambda_{n+1} = \lambda_n + h_n$).

### 16.4 Algorithm

1. Set $M$, $\epsilon$, $E$, $L$, $r_0$, $\phi_0$, $t_0$, $h_0$, $\alpha$, $\delta$, $r_{\text{max}}$ and the maximum number of steps $N$.
2. Compute $p_{r,0}$ from (28) (or from (29)/(30) for the test cases).
3. Form $\mathbf{y}_0 = (t_0, r_0, \phi_0, p_{r,0})$.
4. For $n = 0, 1, \ldots, N-1$:
   1. Compute $h_n$ from (37).
   2. Compute $\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3, \mathbf{k}_4$ from (35) as in Section 16.2.
   3. Update $\mathbf{y}_{n+1}$ with (36) and $\lambda_{n+1} = \lambda_n + h_n$.
   4. Store $\mathbf{y}_{n+1}$ and compute $\mathcal{C}$ from (31).
   5. Stop if a stopping criterion (Section 14.5) is met.
5. Convert $(r, \phi)$ to $(x, y)$ with (32) for visualization.

### 16.5 Implementation

The algorithm is implemented in the [`schwarzschild`](../src/schwarzschild) package. This table maps the equations of this document to the code:

| Equations | Code |
|---|---|
| (15) | [`metric.effective_potential`](../src/schwarzschild/metric.py) |
| (16), (21) | [`metric.Schwarzschild.f`, `.df`](../src/schwarzschild/metric.py) |
| (25) (equals (24) for Schwarzschild) | [`equations.rhs`](../src/schwarzschild/equations.py) |
| (28) | [`initial_conditions.initial_radial_momentum`](../src/schwarzschild/initial_conditions.py) |
| (29) | [`analytic.circular_orbit_constants`](../src/schwarzschild/analytic.py) |
| (30) | [`analytic.critical_impact_parameter`](../src/schwarzschild/analytic.py) |
| (31) | [`equations.hamiltonian_constraint`, `equations.relative_constraint_error`](../src/schwarzschild/equations.py) |
| (32) | [`integrators.Trajectory.xy`](../src/schwarzschild/integrators.py) |
| (34), (36) | [`integrators.rk4_step`](../src/schwarzschild/integrators.py) |
| (37), Section 14.5 | [`integrators.integrate`](../src/schwarzschild/integrators.py) |

### 16.6 Validation

The test suite ([`tests/`](../tests)) and the [notebook](../notebooks/schwarzschild_simulation.ipynb) compare the numerical results with closed-form predictions:

* The circular orbit $r_c = 10M$ stays at $r = 10M$, and the numerical $d\phi/dt = 0.0316227766$ equals $\sqrt{M/r_c^3}$.
* The bound orbit oscillates in $r \in [10.3308M,\; 20M]$ with $|\mathcal{C}| \sim 10^{-16}$. Its periapsis advance per orbit, $2.127094008$ rad, agrees with the exact elliptic-integral result $2.127093985$ rad (below) to a relative difference of $10^{-8}$.
* The photon with $b = 5.19M$ is captured, while $b = 5.20M$ escapes after circling the photon sphere down to $r_{\min} \approx 3.07M$, consistent with $b_c = 3\sqrt{3}M \approx 5.196M$ in (30).
* For large $b$ the deflection angle follows the weak-field series $\alpha \approx 4M/b + \frac{15\pi}{4}(M/b)^2 + \frac{128}{3}(M/b)^3 + \frac{3465\pi}{64}(M/b)^4$; for $b = 40M$ the numerical value is $0.108104$ against $0.108096$ from the series.
* The measured global convergence order of RK4 is $3.996$ (theory: $4$).

**Exact periapsis advance.** With $u = 1/r$, the orbit equation of a massive particle is $(du/d\phi)^2 = 2M(u - u_1)(u - u_2)(u - u_3)$. The turning points $u_1 = 1/r_{\text{apo}}$ and $u_2 = 1/r_{\text{peri}}$ fix the third root through $u_1 + u_2 + u_3 = 1/(2M)$, and the angle swept in one radial period is

$$
\Delta\phi = \frac{4\,K(k^2)}{\sqrt{2M(u_3 - u_1)}}, \qquad k^2 = \frac{u_2 - u_1}{u_3 - u_1},
$$

where $K$ is the complete elliptic integral of the first kind. The periapsis advance per orbit is $\Delta\phi - 2\pi$ (function [`analytic.precession_per_orbit`](../src/schwarzschild/analytic.py)).
