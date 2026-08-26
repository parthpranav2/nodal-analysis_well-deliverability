# Formulas in Use — Nodal Analysis & Well Deliverability

> **Project:** Nodal Analysis of Petroleum Production Wells  
> **Notebook:** `Nodal_Analysis_Petroleum_Production.ipynb`  
> **Scope:** Black-Oil PVT, Inflow Performance Relationship (IPR), Vertical Lift Performance (VLP/TPC), Multi-Phase Flow Dynamics, Nodal Analysis, and Well Survey Geometry

---

## Table of Contents

1. [Well Survey Geometry (Minimum Curvature Method)](#1-well-survey-geometry-minimum-curvature-method)
2. [Black-Oil PVT — Oil Phase](#2-black-oil-pvt--oil-phase)
3. [Black-Oil PVT — Gas Phase](#3-black-oil-pvt--gas-phase)
4. [Black-Oil PVT — Water Phase](#4-black-oil-pvt--water-phase)
5. [Inflow Performance Relationship (IPR)](#5-inflow-performance-relationship-ipr)
6. [Vertical Lift Performance / Tubing Performance Curve (VLP/TPC)](#6-vertical-lift-performance--tubing-performance-curve-vlptpc)
7. [Multi-Phase Flow Correlations](#7-multi-phase-flow-correlations)
8. [Nodal Analysis & Operating Point](#8-nodal-analysis--operating-point)
9. [Absolute Open Flow (AOF)](#9-absolute-open-flow-aof)
10. [Sensitivity & Optimization](#10-sensitivity--optimization)
11. [Summary Table of Governing Correlations](#11-summary-table-of-governing-correlations)

---

## 1. Well Survey Geometry (Minimum Curvature Method)

The **Minimum Curvature Method** is the standard mathematical approach for computing 3D well trajectories (True Vertical Depth, Northing, Easting) from discrete directional survey measurements: Measured Depth ($MD$), Inclination ($\alpha$), and Azimuth ($\phi$).

### 1.1 Dogleg Angle ($\delta$)

Conventional fomula:
$$\delta = \arccos \Bigl( \cos(\phi_2 - \phi_1)\sin\alpha_1\sin\alpha_2 + \cos\alpha_1\cos\alpha_2\Bigr)$$

The total angular displacement $\delta$ between consecutive survey stations $1$ and $2$ is given by (used):

$$\delta = \arccos \Bigl( \cos(\alpha_2 - \alpha_1) - \sin\alpha_1 \sin\alpha_2 \bigl( 1 - \cos(\phi_2 - \phi_1) \bigr) \Bigr)$$

where:
- $\alpha_1, \alpha_2$ = inclination angles at stations 1 and 2 (radians)
- $\phi_1, \phi_2$ = azimuth angles at stations 1 and 2 (radians)

### 1.2 Ratio Factor ($RF$)

The ratio factor straightens the circular arc into coordinate displacement projections:

$$RF = \begin{cases} 
\dfrac{2}{\delta} \tan\left(\dfrac{\delta}{2}\right), & \delta > 10^{-9} \\[10pt]
1.0, & \delta \le 10^{-9} 
\end{cases}$$

### 1.3 Coordinate Increments

For an interval length $\Delta MD = MD_2 - MD_1$:

$$\Delta \text{TVD} = \frac{\Delta MD}{2} \, (\cos\alpha_1 + \cos\alpha_2) \cdot RF$$

$$\Delta \text{North} = \frac{\Delta MD}{2} \, (\sin\alpha_1 \cos\phi_1 + \sin\alpha_2 \cos\phi_2) \cdot RF$$

$$\Delta \text{East} = \frac{\Delta MD}{2} \, (\sin\alpha_1 \sin\phi_1 + \sin\alpha_2 \sin\phi_2) \cdot RF$$

### 1.4 Coordinate Accumulation

$$\text{TVD}_n = \sum_{i=1}^n \Delta \text{TVD}_i, \quad \text{North}_n = \sum_{i=1}^n \Delta \text{North}_i, \quad \text{East}_n = \sum_{i=1}^n \Delta \text{East}_i$$

### 1.5 Unit Conversion & Angle Conventions

$$\text{MD}_{\text{ft}} = 3.28084 \cdot \text{MD}_{\text{m}}$$

$$\theta_{\text{code}} = 90^\circ - \theta_{\text{survey}}$$

---

## 2. Black-Oil PVT — Oil Phase

### 2.1 API to Specific Gravity ($\gamma_o$)

$$\gamma_o = \frac{141.5}{131.5 + \text{API}}$$

### 2.2 Standing (1947) Solution Gas-Oil Ratio ($R_s$) — Saturated ($P \le P_b$)

$$a = 0.0125 \cdot \text{API} - 0.00091 \cdot T$$

$$R_s = \gamma_g \left[ \left( \frac{P}{18.2} + 1.4 \right) \cdot 10^a \right]^{1.2048} \quad [\text{scf/STB}]$$

where:
- $P$ = pressure ($\text{psia}$)
- $T$ = reservoir temperature ($^\circ\text{F}$)
- $\gamma_g$ = gas specific gravity ($\text{air} = 1.0$)
- $\text{API}$ = oil API gravity ($^\circ\text{API}$)

### 2.3 Standing (1947) Bubble-Point Pressure ($P_b$)

Inverting the solution GOR correlation for a given initial dissolved GOR $R_{sb}$:

$$P_b = 18.2 \left[ \left( \frac{R_{sb}}{\gamma_g} \right)^{\frac{1}{1.2048}} \cdot 10^{-a} - 1.4 \right] \quad [\text{psia}]$$

### 2.4 Standing (1947) Saturated Oil Formation Volume Factor ($B_o$)

$$F = R_s \left( \frac{\gamma_g}{\gamma_o} \right)^{0.5} + 1.25 \cdot T$$

$$B_o = 0.9759 + 1.2 \times 10^{-4} \cdot F^{1.2} \quad [\text{rb/STB}]$$

### 2.5 Vasquez & Beggs (1980) Undersaturated Oil Compressibility ($c_o$)

For $P > P_b$:

$$c_o = \frac{-1433 + 5 R_{sb} + 17.2 T - 1180 \gamma_g + 12.61 \text{API}}{10^5 \cdot P} \quad [\text{psi}^{-1}]$$

### 2.6 Undersaturated Oil Formation Volume Factor ($B_o$)

For $P > P_b$:

$$B_o = B_{ob} \cdot \exp \bigl[ c_o (P_b - P) \bigr] \quad [\text{rb/STB}]$$

where $B_{ob} = B_o(P_b, R_{sb})$.

### 2.7 Live Oil Density ($\rho_o$)

$$\rho_o = \frac{350.17 \cdot \gamma_o + 0.0764 \cdot R_s \cdot \gamma_g}{5.615 \cdot B_o} \quad [\text{lbm/ft}^3]$$

### 2.8 Beggs & Robinson (1975) Dead-Oil Viscosity ($\mu_{od}$)

$$Z = 3.0324 - 0.02023 \cdot \text{API}$$

$$Y = 10^Z$$

$$X = Y \cdot T^{-1.163}$$

$$\mu_{od} = 10^X - 1 \quad [\text{cp}]$$

### 2.9 Beggs & Robinson (1975) Saturated Live-Oil Viscosity ($\mu_{ob}$)

$$\mu_{ob} = a \cdot \mu_{od}^b \quad [\text{cp}]$$

where:

$$a = 10.715 \cdot (R_s + 100)^{-0.515}$$

$$b = 5.44 \cdot (R_s + 150)^{-0.338}$$

### 2.10 Beggs & Robinson (1975) Undersaturated Oil Viscosity ($\mu_o$)

For $P > P_b$:

$$m = 2.6 \cdot P^{1.187} \cdot \exp \bigl( -11.513 - 8.98 \times 10^{-5} \cdot P \bigr)$$

$$\mu_o = \max \left( \mu_{ob} \left( \frac{P}{P_b} \right)^m, \; 0.05 \right) \quad [\text{cp}]$$

### 2.11 Baker & Swerdloff (1956) Interfacial Tension ($\sigma$)

$$\sigma_{68} = 39.0 - 0.2571 \cdot \text{API} \quad [\text{dynes/cm}]$$

$$\sigma_{100} = 37.5 - 0.2571 \cdot \text{API} \quad [\text{dynes/cm}]$$

$$\sigma_{\text{dead}} = \sigma_{68} - (T_c - 68) \cdot \frac{\sigma_{68} - \sigma_{100}}{32}, \quad T_c = \min \bigl( \max(T, 68), 100 \bigr)$$

Live-oil pressure correction:

$$\sigma_{\text{live}} = \max \left( \sigma_{\text{dead}} \cdot \bigl( 1 - 0.024 \cdot P^{0.45} \bigr), \; 1.0 \right) \quad [\text{dynes/cm}]$$

---

## 3. Black-Oil PVT — Gas Phase

### 3.1 Sutton (1985) Pseudo-Critical Properties

$$T_{pc} = 169.2 + 349.5 \cdot \gamma_g - 74.0 \cdot \gamma_g^2 \quad [^\circ\text{R}]$$

$$P_{pc} = 756.8 - 131.0 \cdot \gamma_g - 3.6 \cdot \gamma_g^2 \quad [\text{psia}]$$

### 3.2 Pseudo-Reduced Properties

$$T_{pr} = \frac{T + 460}{T_{pc}}, \quad P_{pr} = \frac{P}{P_{pc}}$$

### 3.3 Dranchuk & Abou-Kassem (1975) Real-Gas Z-Factor (DAK)

The EOS relation is solved implicitly for $Z$:

$$Z = 1 + c_1 \rho_r + c_2 \rho_r^2 - c_3 \rho_r^5 + A_{10} \bigl(1 + A_{11} \rho_r^2\bigr) \left( \frac{\rho_r^2}{T_{pr}^3} \right) \exp \bigl( -A_{11} \rho_r^2 \bigr)$$

where:

$$\rho_r = \frac{0.27 \cdot P_{pr}}{Z \cdot T_{pr}}$$

$$c_1 = A_1 + \frac{A_2}{T_{pr}} + \frac{A_3}{T_{pr}^3} + \frac{A_4}{T_{pr}^4} + \frac{A_5}{T_{pr}^5}$$

$$c_2 = A_6 + \frac{A_7}{T_{pr}} + \frac{A_8}{T_{pr}^2}$$

$$c_3 = A_9 \left( \frac{A_7}{T_{pr}} + \frac{A_8}{T_{pr}^2} \right)$$

#### DAK Coefficients Table

| Coefficient | Value | Coefficient | Value |
| :--- | :--- | :--- | :--- |
| $A_1$ | $+0.3265$ | $A_7$ | $-0.7361$ |
| $A_2$ | $-1.0700$ | $A_8$ | $+0.1844$ |
| $A_3$ | $-0.5339$ | $A_9$ | $+0.1056$ |
| $A_4$ | $+0.01569$ | $A_{10}$ | $+0.6134$ |
| $A_5$ | $-0.05165$ | $A_{11}$ | $+0.7210$ |
| $A_6$ | $+0.5475$ | | |

### 3.4 Gas Formation Volume Factor ($B_g$)

$$B_g = \frac{0.02827 \cdot Z \cdot (T + 460)}{P} \quad [\text{rcf/scf}]$$

### 3.5 Gas Density ($\rho_g$)

$$\rho_g = \frac{2.7 \cdot \gamma_g \cdot P}{Z \cdot (T + 460)} \quad [\text{lbm/ft}^3]$$

### 3.6 Lee, Gonzalez & Eakin (1966) Gas Viscosity ($\mu_g$)

$$M_a = 28.97 \cdot \gamma_g$$

$$\rho_{g,\text{gcc}} = \frac{\rho_g}{62.4} \quad [\text{g/cm}^3]$$

$$K = \frac{(9.4 + 0.02 M_a) \cdot (T + 460)^{1.5}}{209 + 19 M_a + (T + 460)}$$

$$X = 3.5 + \frac{986}{T + 460} + 0.01 M_a$$

$$Y = 2.4 - 0.2 X$$

$$\mu_g = 10^{-4} \cdot K \cdot \exp \left( X \cdot \rho_{g,\text{gcc}}^Y \right) \quad [\text{cp}]$$

---

## 4. Black-Oil PVT — Water Phase

$$B_w = 1.0 \quad [\text{rb/STB}]$$

$$\rho_w = 62.4 \times 1.07 = 66.768 \quad [\text{lbm/ft}^3]$$

$$\mu_w = \max \bigl( 1.0 - 0.012 \cdot (T - 60), \; 0.30 \bigr) \quad [\text{cp}]$$

---

## 5. Inflow Performance Relationship (IPR)

### 5.1 Straight-Line (Darcy / Productivity Index) IPR — $\bar{P} \ge P_b$

$$q = J \cdot (\bar{P} - P_{wf})$$

$$J = \frac{q_{\text{test}}}{\bar{P} - P_{wf,\text{test}}} \quad [\text{STB/day/psi}]$$

Semi-steady state radial inflow:

$$J = \frac{0.00708 \cdot k \cdot h}{\mu_o \cdot B_o \cdot \left[ \ln\left( \dfrac{r_e}{r_w} \right) - 0.75 + S \right]}$$

Absolute Open Flow ($P_{wf} = 0$):

$$q_{\text{AOF}} = J \cdot \bar{P}$$

### 5.2 Vogel (1968) IPR — Two-Phase Reservoir ($\bar{P} \le P_b$)

$$\frac{q}{q_{\max}} = 1 - 0.2 \left( \frac{P_{wf}}{\bar{P}} \right) - 0.8 \left( \frac{P_{wf}}{\bar{P}} \right)^2$$

$$q_{\max} = \frac{q_{\text{test}}}{1 - 0.2 \left( \dfrac{P_{wf,\text{test}}}{\bar{P}} \right) - 0.8 \left( \dfrac{P_{wf,\text{test}}}{\bar{P}} \right)^2}$$

### 5.3 Composite (Klins & Clark, 1993) IPR — Saturated/Undersaturated Mixed ($\bar{P} > P_b$)

$$\text{For } P_{wf} \ge P_b: \quad q = J \cdot (\bar{P} - P_{wf})$$

$$\text{For } P_{wf} < P_b: \quad q = q_b + \frac{J \cdot P_b}{1.8} \left[ 1 - 0.2 \left( \frac{P_{wf}}{P_b} \right) - 0.8 \left( \frac{P_{wf}}{P_b} \right)^2 \right]$$

where:

$$q_b = J \cdot (\bar{P} - P_b)$$

### 5.4 Jones, Blount & Glaze (1976) Non-Darcy Turbulent IPR

$$\bar{P} - P_{wf} = a \cdot q + b \cdot q^2$$

where:
- $a$ = laminar flow resistance coefficient ($\text{psi} / (\text{STB/day})$)
- $b$ = turbulence/non-Darcy coefficient ($\text{psi} / (\text{STB/day})^2$)

Absolute Open Flow:

$$q_{\text{AOF}} = \frac{-a + \sqrt{a^2 + 4 b \bar{P}}}{2 b}$$

### 5.5 Skin Factor ($S$) & Flow Efficiency ($FE$)

$$S = \left( \frac{J_{\text{ideal}}}{J_{\text{actual}}} - 1 \right) \cdot \left[ \ln\left( \frac{r_e}{r_w} \right) - 0.75 \right]$$

$$FE = \frac{J_{\text{actual}}}{J_{\text{ideal}}} = \frac{\bar{P} - P_{wf,\text{ideal}}}{\bar{P} - P_{wf,\text{actual}}}$$

---

## 6. Vertical Lift Performance / Tubing Performance Curve (VLP/TPC)

### 6.1 Mechanical Energy Conservation Gradient Breakdown

The steady-state mechanical energy conservation equation along measured depth $z$ ($MD$) accounts for three simultaneous pressure dissipation mechanisms:

$$\frac{dP}{dz} = \left( \frac{dP}{dz} \right)_{\text{gravity}} + \left( \frac{dP}{dz} \right)_{\text{friction}} + \left( \frac{dP}{dz} \right)_{\text{kinetic}}$$

1. **Gravitational / Hydrostatic Gradient:**
   $$\left( \frac{dP}{dz} \right)_{\text{gravity}} = \frac{\rho_m \cdot g \cdot \sin\theta}{144 \cdot g_c} \quad [\text{psi/ft}]$$
   where $\theta$ is the well inclination angle relative to horizontal ($90^\circ$ for vertical, $0^\circ$ for horizontal in Beggs-Brill notation).

2. **Frictional Gradient:**
   $$\left( \frac{dP}{dz} \right)_{\text{friction}} = \frac{f_{tp} \cdot \rho_{ns} \cdot v_m^2}{2 \cdot 144 \cdot g_c \cdot d} \quad [\text{psi/ft}]$$

3. **Kinetic / Acceleration Gradient:**
   $$\left( \frac{dP}{dz} \right)_{\text{kinetic}} = E_k \cdot \left( \frac{dP}{dz} \right)_{\text{total}} = \frac{E_k}{1 - E_k} \left[ \left( \frac{dP}{dz} \right)_{\text{gravity}} + \left( \frac{dP}{dz} \right)_{\text{friction}} \right] \quad [\text{psi/ft}]$$
   where $E_k$ is the dimensionless kinetic acceleration parameter:
   $$E_k = \min\left( \frac{v_m \cdot v_{sg} \cdot \rho_{ns}}{g_c \cdot P \cdot 144}, \; 0.6 \right)$$

Total coupled gradient:
$$\left( \frac{dP}{dz} \right)_{\text{total}} = \frac{\left( \dfrac{dP}{dz} \right)_{\text{gravity}} + \left( \dfrac{dP}{dz} \right)_{\text{friction}}}{1 - E_k} \quad [\text{psi/ft}]$$

where:
- $g = 32.174 \;\text{ft/s}^2$, $g_c = 32.174 \;\text{lbm}\cdot\text{ft}/(\text{lbf}\cdot\text{s}^2)$
- $d$ = internal diameter of tubing ($\text{ft}$)
- $P$ = in-situ pressure ($\text{psia}$)
- $\rho_{ns}, \rho_m$ = no-slip and slip mixture densities ($\text{lbm/ft}^3$)
- $v_m, v_{sg}$ = mixture and superficial gas velocities ($\text{ft/s}$)

### 6.2 Multi-Phase In-Situ Volumetric Rates

$$q_l = 5.615 \cdot \frac{q_o B_o + q_w B_w}{86400} \quad [\text{ft}^3/\text{s}]$$

$$q_g = \frac{\max(GOR - R_s, \; 0) \cdot q_o \cdot B_g}{86400} \quad [\text{ft}^3/\text{s}]$$

### 6.3 Superficial and Mixture Velocities

$$v_{sl} = \frac{q_l}{A_t} = \frac{4 q_l}{\pi d^2} \quad [\text{ft/s}], \qquad v_{sg} = \frac{q_g}{A_t} = \frac{4 q_g}{\pi d^2} \quad [\text{ft/s}]$$

$$v_m = v_{sl} + v_{sg} \quad [\text{ft/s}]$$

### 6.4 Liquid Holdup & Mixture Density

$$\rho_m = H_l \rho_l + (1 - H_l) \rho_g \quad [\text{lbm/ft}^3]$$

$$\rho_l = f_w \rho_w + (1 - f_w) \rho_o, \quad f_w = \frac{\text{WOR}}{1 + \text{WOR}} = \text{WC}$$

### 6.5 Reynolds Number ($Re$) & Colebrook-White Friction Factor ($f_m$)

$$Re = \frac{1488 \cdot \rho_m \cdot v_m \cdot d}{\mu_m}$$

$$\mu_m = H_l \mu_l + (1 - H_l) \mu_g$$

Colebrook-White relationship:

$$\frac{1}{\sqrt{f_m}} = -2 \log_{10} \left( \frac{\varepsilon / d}{3.7} + \frac{2.51}{Re \sqrt{f_m}} \right)$$

### 6.6 Bottomhole Pressure Traversal Integration

$$P_{wf} = P_{wh} + \int_0^{MD} \left( \frac{dP}{dz} \right) dz$$

---

## 7. Multi-Phase Flow Correlations

### 7.1 Hagedorn & Brown (1965) Dimensionless Groups

1. **Liquid Velocity Number ($N_{vl}$):**
   $$N_{vl} = 1.938 \cdot v_{sl} \left( \frac{\rho_l}{\sigma} \right)^{0.25}$$

2. **Gas Velocity Number ($N_{vg}$):**
   $$N_{vg} = 1.938 \cdot v_{sg} \left( \frac{\rho_l}{\sigma} \right)^{0.25}$$

3. **Pipe Diameter Number ($N_d$):**
   $$N_d = 120.872 \cdot d \left( \frac{\rho_l}{\sigma} \right)^{0.5}$$

4. **Liquid Viscosity Number ($N_l$):**
   $$N_l = 0.15726 \cdot \mu_l \left( \frac{1}{\rho_l \cdot \sigma^3} \right)^{0.25}$$

### 7.2 Beggs & Brill (1973) Multi-Angle Correlation

- **No-Slip Holdup ($\lambda_l$):**
  $$\lambda_l = \frac{v_{sl}}{v_{sl} + v_{sg}} = \frac{v_{sl}}{v_m}$$

- **Froude Number ($N_{Fr}$):**
  $$N_{Fr} = \frac{v_m^2}{g \cdot d}$$

- **Flow Pattern Boundaries:**
  $$L_1 = 316 \cdot \lambda_l^{0.302}, \quad L_2 = 0.0009252 \cdot \lambda_l^{-2.4684}$$
  $$L_3 = 0.10 \cdot \lambda_l^{-1.4516}, \quad L_4 = 0.5 \cdot \lambda_l^{-6.738}$$

- **Inclination Correction Factor ($\psi$):**
  $$H_l(\theta) = H_l(0^\circ) \cdot \psi(\theta)$$
  $$\psi(\theta) = 1 + C \left[ \sin(1.8 \theta) - \frac{1}{3} \sin^3(1.8 \theta) \right]$$

### 7.3 Griffith & Wallis (1961) Drift-Flux Slug Flow Model

$$v_d = 0.35 \sqrt{g \cdot d}$$

$$H_l = 1 - \frac{v_{sg}}{v_{sl} + v_{sg} + v_d}$$

---

## 8. Nodal Analysis & Operating Point

### 8.1 Bottomhole Node Balance

$$\Delta P_{\text{node}}(q) = P_{wf,\text{IPR}}(q) - P_{wf,\text{VLP}}(q) = 0$$

Operating flow rate $q_{\text{op}}$ and pressure $P_{wf,\text{op}}$ are solved using Brent's method:

$$q_{\text{op}} = \text{root} \bigl\{ P_{wf,\text{IPR}}(q) - P_{wf,\text{VLP}}(q) \bigr\}$$

$$P_{wf,\text{op}} = P_{wf,\text{IPR}}(q_{\text{op}}) = P_{wf,\text{VLP}}(q_{\text{op}})$$

### 8.2 Total Energy Pressure Drop Decomposition

$$\bar{P} - P_{\text{sep}} = \Delta P_{\text{reservoir}} + \Delta P_{\text{completion}} + \Delta P_{\text{tubing}} + \Delta P_{\text{choke}} + \Delta P_{\text{surface line}}$$

---

## 9. Absolute Open Flow (AOF)

### 9.1 Summary of AOF by IPR Model

1. **Darcy PI Model:**
   $$q_{\text{AOF}} = J \cdot \bar{P}$$

2. **Vogel (1968) Model:**
   $$q_{\text{AOF}} = q_{\max} = \frac{q_{\text{test}}}{1 - 0.2 \left( \dfrac{P_{wf,\text{test}}}{\bar{P}} \right) - 0.8 \left( \dfrac{P_{wf,\text{test}}}{\bar{P}} \right)^2}$$

3. **Composite Klins & Clark Model:**
   $$q_{\text{AOF}} = J(\bar{P} - P_b) + \frac{J \cdot P_b}{1.8}$$

4. **Jones Non-Darcy Model:**
   $$q_{\text{AOF}} = \frac{-a + \sqrt{a^2 + 4 b \bar{P}}}{2 b}$$

---

## 10. Sensitivity & Optimization

### 10.1 Water-Cut Effect on Mixture Properties

$$\text{WC} = \frac{\text{WOR}}{1 + \text{WOR}}$$

$$\rho_l = (1 - \text{WC}) \cdot \rho_o + \text{WC} \cdot \rho_w$$

$$\mu_l = \mu_o^{(1 - \text{WC})} \cdot \mu_w^{\text{WC}}$$

### 10.2 Choke Performance Relations

- **Critical (Sonic) Flow:**
  $$q_o = C_d \cdot A_c \cdot \sqrt{\frac{\Delta P_{\text{choke}}}{\rho_l}}$$

- **Subcritical Flow:**
  $$q_o = C_d \cdot A_c \cdot \sqrt{\frac{2 \Delta P_{\text{choke}}}{\rho_l}}$$

---

## 11. Summary Table of Governing Correlations

| Phenomenon / Variable | Correlation / Law | Primary Author & Year |
| :--- | :--- | :--- |
| **Wellbore Trajectory** | Minimum Curvature Method | API RP 11L / Standard (1985) |
| **Solution GOR ($R_s$)** | Empirical Flash Dissolution | Standing (1947) |
| **Bubble Point ($P_b$)** | Inverted Empirical GOR | Standing (1947) |
| **Oil FVF ($B_o$)** | Saturated Black-Oil | Standing (1947) |
| **Oil Compressibility ($c_o$)** | Undersaturated Regime | Vasquez & Beggs (1980) |
| **Dead-Oil Viscosity ($\mu_{od}$)** | Low-Pressure Empirical | Beggs & Robinson (1975) |
| **Live-Oil Viscosity ($\mu_{ob}, \mu_o$)** | Live/Undersaturated Multiplier | Beggs & Robinson (1975) |
| **Interfacial Tension ($\sigma$)** | Dead/Live Surface Tension | Baker & Swerdloff (1956) |
| **Gas Pseudo-Criticals** | Polynomial Natural Gas | Sutton (1985) |
| **Gas $Z$-Factor** | Implicit 11-Constant EOS | Dranchuk & Abou-Kassem (1975) |
| **Gas Viscosity ($\mu_g$)** | Semi-Empirical Kinetic Model | Lee, Gonzalez & Eakin (1966) |
| **Linear IPR** | Radial Semi-Steady Darcy Flow | Darcy (1856) |
| **Two-Phase IPR** | Numerical Inflow Simulation | Vogel (1968) |
| **Mixed IPR** | Composite Continuity Formulation | Klins & Clark (1993) |
| **Non-Darcy IPR** | Quadratic Forchheimer Inflow | Jones, Blount & Glaze (1976) |
| **Vertical VLP** | Dimensionless Multi-Phase Holdup | Hagedorn & Brown (1965) |
| **Deviated VLP** | Multi-Angle Flow Regime Map | Beggs & Brill (1973) |
| **Slug Holdup** | Drift-Flux Kinematics | Griffith & Wallis (1961) |
| **Pipe Friction Factor** | Implicit Rough Pipe Relation | Colebrook & White (1937) |
| **Nodal Solution** | Derivative-Free Root Finding | Brent (1973) |

---
*All standard field units apply: Pressure in $\text{psia}$, Temperature in $^\circ\text{F}$ / $^\circ\text{R}$, Rates in $\text{STB/day}$ and $\text{Mscf/day}$, Lengths/Depths in $\text{ft}$, Viscosity in $\text{cp}$, and Interfacial Tension in $\text{dynes/cm}$.*
