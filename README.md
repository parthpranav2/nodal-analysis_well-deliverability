# Nodal Analysis & Well Deliverability Simulation System

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Jupyter](https://img.shields.io/badge/Notebook-Jupyter-orange.svg)](https://jupyter.org/)
[![Domain](https://img.shields.io/badge/Domain-Petroleum%20Production%20Engineering-red.svg)]()
[![Formulas](https://img.shields.io/badge/Documentation-FORMULAS__IN__USE.md-brightgreen.svg)](FORMULAS_IN_USE.md)

An end-to-end, physics-driven **Nodal Analysis & Well Deliverability** modelling suite for vertical, horizontal, directional, and custom trajectory oil and gas production wells. 

The system couples reservoir inflow deliverability (**IPR**), fluid thermophysical properties (**Black-Oil PVT**), and wellbore hydraulics (**VLP / TPC**) via mechanistic and empirical multi-phase flow correlations, resolving the bottomhole nodal operating point and providing comprehensive optimization diagnostics (tubing sizing, water-cut loading, GOR variation, reservoir depletion, artificial gas-lift, and critical liquid loading).

---

## Table of Contents

- [1. Executive Summary & Physics Architecture](#1-executive-summary--physics-architecture)
- [2. Mathematical & Engineering Foundations](#2-mathematical--engineering-foundations)
  - [2.1 Black-Oil PVT Fluid Characterization](#21-black-oil-pvt-fluid-characterization)
  - [2.2 Inflow Performance Relationship (IPR)](#22-inflow-performance-relationship-ipr)
  - [2.3 Vertical Lift Performance (VLP / TPC) & Multi-Phase Flow](#23-vertical-lift-performance-vlp--tpc--multi-phase-flow)
  - [2.4 Well Trajectory & Minimum Curvature Integration](#24-well-trajectory--minimum-curvature-integration)
  - [2.5 Liquid Loading & Critical Unloading Velocity](#25-liquid-loading--critical-unloading-velocity)
- [3. Test Run Results & Benchmark Analysis](#3-test-run-results--benchmark-analysis)
  - [3.1 Test Case 1: Vertical Well (Baseline)](#31-test-case-1-vertical-well-baseline)
  - [3.2 Test Case 2: Horizontal Lateral Well](#32-test-case-2-horizontal-lateral-well)
  - [3.3 Test Case 3: Directional Well (Schoonebeek-1301 Survey Benchmark)](#33-test-case-3-directional-well-schoonebeek-1301-survey-benchmark)
  - [3.4 Cross-Case Comparative Synthesis](#34-cross-case-comparative-synthesis)
- [4. Project Structure & Key Files](#4-project-structure--key-files)
- [5. Installation & Execution Guide](#5-installation--execution-guide)
- [6. References & Standards](#6-references--standards)

---

## 1. Executive Summary & Physics Architecture

Nodal Analysis partitions the entire production system (from the drainage radius in the reservoir to the surface separator) into distinct hydraulic nodes to balance fluid deliverability against energy dissipation.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                               NODAL ANALYSIS WORKFLOW                                  │
│                                                                                        │
│   RESERVOIR INFLOW (IPR)                     WELLBORE HYDRAULICS (VLP)                 │
│   • Standing / Vogel / Composite             • Beggs & Brill (1973) Multi-Angle        │
│   • Jones Non-Darcy Quadratic                • Hagedorn & Brown (1965) Vertical        │
│   • Skin Factor & Flow Efficiency            • Colebrook-White Pipe Friction           │
│                    │                                         │                         │
│                    └─────────────────┬───────────────────────┘                         │
│                                      ▼                                                 │
│                        BOTTOMHOLE NODE PRESSURE BALANCE                                │
│                     Pwf_IPR(q*) = Pwf_VLP(q*)  →  Operating Point (q*, Pwf*)           │
│                                      │                                                 │
│                                      ▼                                                 │
│                             SENSITIVITY & OPTIMIZATION                                 │
│   • Tubing Sizing (API 5CT)       • Water Cut & GOR Shifts    • Gas-Lift Optimization  │
│   • Depletion Trajectory          • Wellhead Back-Pressure    • Liquid Loading Limits  │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

For full analytical derivations, coefficient tables, and governing equations, refer to [`FORMULAS_IN_USE.md`](FORMULAS_IN_USE.md).

---

## 2. Mathematical & Engineering Foundations

### 2.1 Black-Oil PVT Fluid Characterization

The fluid engine models live oil, dissolved gas, and water across both saturated ($P \le P_b$) and undersaturated ($P > P_b$) states:

- **Bubble-Point Pressure ($P_b$) & Solution GOR ($R_s$):** Calculated using Standing (1947) correlations.
- **Oil Formation Volume Factor ($B_o$):** Saturated expansion via Standing (1947); undersaturated compressibility $c_o$ via Vasquez & Beggs (1980).
- **Live Oil Viscosity ($\mu_o$):** Dead-oil base via Beggs & Robinson (1975), live saturated reduction via dissolved $R_s$, and exponential scaling above $P_b$.
- **Real Gas $Z$-Factor:** Implicit 11-constant Dranchuk & Abou-Kassem (1975) equation of state solved numerically via Brent's method with Sutton (1985) pseudo-criticals ($T_{pc}, P_{pc}$).
- **Gas Viscosity ($\mu_g$):** Lee, Gonzalez & Eakin (1966) semi-empirical relation.
- **Interfacial Tension ($\sigma$):** Baker & Swerdloff (1956) with Standing live-oil pressure correction.

### 2.2 Inflow Performance Relationship (IPR)

- **Linear Darcy IPR:** $q = J \cdot (\bar{P} - P_{wf})$ when $\bar{P} \ge P_b$ and $P_{wf} \ge P_b$.
- **Vogel (1968) Two-Phase IPR:** $q / q_{\max} = 1 - 0.2(P_{wf}/\bar{P}) - 0.8(P_{wf}/\bar{P})^2$ for dissolved gas breakout.
- **Composite Klins & Clark (1993) IPR:** Enforces slope continuity across the bubble-point kink for saturated/undersaturated transitions.
- **Jones, Blount & Glaze (1976) Non-Darcy IPR:** Quadratic Forchheimer flow $\bar{P} - P_{wf} = a q + b q^2$ accounting for near-wellbore turbulence.

### 2.3 Vertical Lift Performance (VLP / TPC) & Multi-Phase Flow

The pressure gradient along measured depth $z$ balances hydrostatic head and wall shear friction:

$$\frac{dP}{dz} = \frac{\rho_m g \cos\theta}{144 g_c} + \frac{f_m \rho_m v_m^2}{2 \cdot 144 g_c d}$$

- **Beggs & Brill (1973):** Evaluates horizontal, vertical, and inclined flow regimes (Segregated, Intermittent, Distributed, Transition) using Froude number $N_{Fr}$ and no-slip holdup $\lambda_l$, modified by inclination factor $\psi(\theta)$.
- **Hagedorn & Brown (1965):** Uses dimensionless velocity, diameter, and viscosity numbers for vertical lift.
- **Colebrook-White:** Friction factor $f_m$ solved iteratively across laminar, transitional, and fully rough turbulent regimes.

### 2.4 Well Trajectory & Minimum Curvature Integration

Directional trajectories from survey stations ($MD_i, \alpha_i, \phi_i$) are calculated using the 3D **Minimum Curvature Method**:

$$\delta = \arccos \bigl[ \cos(\alpha_2 - \alpha_1) - \sin\alpha_1 \sin\alpha_2 (1 - \cos(\phi_2 - \phi_1)) \bigr]$$

$$RF = \frac{2}{\delta} \tan\left(\frac{\delta}{2}\right)$$

Coordinate increments $\Delta\text{TVD}, \Delta\text{North}, \Delta\text{East}$ are accumulated station-by-station to map local inclination angles $\theta(z)$ for accurate gravity-gradient resolution.

### 2.5 Liquid Loading & Critical Unloading Velocity

As reservoir pressure depletes or flow velocity drops, the gas phase can no longer lift the liquid film, causing liquid accumulation and well loading. The system determines:
- **VLP Minimum ($q_{\min}, P_{wf,\min}$):** The inflection point separating friction-dominated flow from hydrostatic-dominated loading.
- **Critical Reservoir Pressure ($\bar{P}_{\text{crit}}$):** The threshold pressure below which stable natural flow ceases.

---

## 3. Test Run Results & Benchmark Analysis

The `results/` directory stores three comprehensive test runs evaluating different completion geometries and survey trajectories.

---

### 3.1 Test Case 1: Vertical Well (Baseline)

- **Directory:** [`results/vertical/`](results/vertical/)
- **Configuration:** Vertical well ($\theta = 90^\circ$, $\text{TVD} = 8,500\text{ ft}$), Default API 5CT Tubing 2-7/8" ($2.441\text{ in ID}$), $P_r = 4,200\text{ psia}$, $P_b = 2,463\text{ psia}$, $\text{WC} = 20\%$, $\text{GOR} = 600\text{ scf/STB}$, $P_{wh} = 300\text{ psia}$.
- **Baseline Operating Point:** **$q^* = 2,267\text{ STB/d}$** (1,814 STB/d oil, 453 STB/d water) at **$P_{wf}^* = 2,805\text{ psia}$**.
- **AOF Potential:** $5,046\text{ STB/d}$ ($J = 1.625\text{ STB/d/psi}$).

#### Diagnostic Results & Plots

| Analysis | Output Artifact | Key Finding |
| :--- | :--- | :--- |
| **PVT Characterization** | [`results/vertical/689f15ca-02d2-46d7-949f-11b64387a848.png`](results/vertical/689f15ca-02d2-46d7-949f-11b64387a848.png) | Bubble point $P_b = 2,463\text{ psia}$, $B_{ob} = 1.34\text{ rb/STB}$, $\mu_o = 0.98\text{ cp}$. |
| **Inflow Performance (IPR)** | [`results/vertical/00fece15-f5b4-45c8-b49d-3987b919b966.png`](results/vertical/00fece15-f5b4-45c8-b49d-3987b919b966.png) | Single test calibration at $(1,300\text{ STB/d}, 3,400\text{ psia})$ above $P_b$. |
| **Pressure Traverse & VLP** | [`results/vertical/5d722eea-6fd7-4c2f-ba46-4a46da2fb25f.png`](results/vertical/5d722eea-6fd7-4c2f-ba46-4a46da2fb25f.png) | Continuous gradient traversal from $300\text{ psia}$ wellhead to $2,805\text{ psia}$ BHP. |
| **Base-Case Nodal Solution** | [`results/vertical/4e6efc02-875a-4822-b3cd-263ceee686f7.png`](results/vertical/4e6efc02-875a-4822-b3cd-263ceee686f7.png) | Stable operating intersection at $2,267\text{ STB/d}$. |
| **Tubing Size Sensitivity** | [`results/vertical/6a0e6a28-f2e4-4a58-8aa6-b7ef779a9a26.png`](results/vertical/6a0e6a28-f2e4-4a58-8aa6-b7ef779a9a26.png) | Rate increases from $1,328\text{ STB/d}$ ($1.900"$) up to $2,828\text{ STB/d}$ ($4-1/2"$). |
| **Water Cut Sensitivity** | [`results/vertical/8962eead-3a97-48fa-8967-8f370844a896.png`](results/vertical/8962eead-3a97-48fa-8967-8f370844a896.png) | Rate drops from $2,485\text{ STB/d}$ ($\text{WC}=0\%$) to $877\text{ STB/d}$ ($\text{WC}=80\%$). |
| **GOR Sensitivity** | [`results/vertical/22e79380-9685-4473-b7fd-ddaa4147af3a.png`](results/vertical/22e79380-9685-4473-b7fd-ddaa4147af3a.png) | Higher GOR lowers mixture density, increasing rate from $1,759$ to $2,580\text{ STB/d}$. |
| **Reservoir Depletion** | [`results/vertical/fe364d12-e20c-472c-853f-91c974b6ffbc.png`](results/vertical/fe364d12-e20c-472c-853f-91c974b6ffbc.png)<br>[`results/vertical/8cf0910d-6407-411c-a094-1236584d2355.png`](results/vertical/8cf0910d-6407-411c-a094-1236584d2355.png) | Deliverability tracks decline; ceases natural flow below $P_r \approx 2,939\text{ psia}$. |
| **Wellhead Pressure ($P_{wh}$)** | [`results/vertical/8653b3c2-5ff6-4144-adaa-d7c4428517ad.png`](results/vertical/8653b3c2-5ff6-4144-adaa-d7c4428517ad.png) | Choke backpressure sensitivity: $P_{wh} = 100 \to 1200\text{ psia}$ ($q^* = 2539 \to 369\text{ STB/d}$). |
| **Artificial Gas-Lift** | [`results/vertical/6b620a72-09d7-468b-9dd4-2ccf10706af3.png`](results/vertical/6b620a72-09d7-468b-9dd4-2ccf10706af3.png) | Optimum injection rate: $+1.0\text{ MMscf/d}$ at $5,950\text{ ft}$ yields $q^* = 2,460\text{ STB/d}$. |
| **Liquid Loading Diagnostics** | [`results/vertical/c3030d99-af10-41c3-97f2-9edb5a12be44.png`](results/vertical/c3030d99-af10-41c3-97f2-9edb5a12be44.png) | $q_{\min} = 726\text{ STB/d}$, critical depletion limit $\bar{P}_{\text{crit}} = 2,939\text{ psia}$. |
| **Executive 4-Panel Dashboard** | [`results/vertical/52312abf-2b2a-44b5-9560-803d37bcbfc9.png`](results/vertical/52312abf-2b2a-44b5-9560-803d37bcbfc9.png) | Multi-panel executive summary dashboard. |
| **Execution State Screenshot** | [`results/vertical/Screenshot 2026-08-21 at 12.43.28.png`](results/vertical/Screenshot%202026-08-21%20at%2012.43.28.png) | Interactive notebook UI session state. |

---

### 3.2 Test Case 2: Horizontal Lateral Well

- **Directory:** [`results/horizontal/`](results/horizontal/)
- **Configuration:** Horizontal wellbore ($\theta = 0^\circ$ code angle, lateral section at $\text{TVD} = 8,500\text{ ft}$), Default API 5CT Tubing 3-1/2" ($2.992\text{ in ID}$), matching reservoir & fluid baseline.
- **Physics Mechanism:** Zero gravity pressure loss in the horizontal lateral ($\cos 90^\circ = 0$ in the lateral section) coupled with higher frictional dissipation.

#### Diagnostic Results & Plots

| Analysis | Output Artifact | Key Finding |
| :--- | :--- | :--- |
| **PVT Characterization** | [`results/horizontal/6329b5b9-a8b3-4595-8e24-f4ea20c301dd.png`](results/horizontal/6329b5b9-a8b3-4595-8e24-f4ea20c301dd.png) | Identical baseline fluid system calibration. |
| **Inflow Performance (IPR)** | [`results/horizontal/a1e6ca9e-43ed-4825-bc38-7115b143a36a.png`](results/horizontal/a1e6ca9e-43ed-4825-bc38-7115b143a36a.png) | Standard reservoir deliverability relationship. |
| **Pressure Traverse & VLP** | [`results/horizontal/375b1050-d907-43f7-9f33-7c063b8153e7.png`](results/horizontal/375b1050-d907-43f7-9f33-7c063b8153e7.png) | Lateral friction gradient without hydrostatic head loss. |
| **Base-Case Nodal Solution** | [`results/horizontal/703a4fa7-76f0-498d-98b7-5f24e2907a97.png`](results/horizontal/703a4fa7-76f0-498d-98b7-5f24e2907a97.png) | Operating deliverability in 3-1/2" horizontal lateral. |
| **Tubing Size Sensitivity** | [`results/horizontal/2af5ddc3-4d87-49d1-8a1f-03807a08b8d1.png`](results/horizontal/2af5ddc3-4d87-49d1-8a1f-03807a08b8d1.png) | Large-bore tubing ($4-1/2"$) strongly favored to mitigate lateral friction. |
| **Water Cut Sensitivity** | [`results/horizontal/4eb466b9-5696-45d1-a163-363c09b2a711.png`](results/horizontal/4eb466b9-5696-45d1-a163-363c09b2a711.png) | Water-cut impact on horizontal multi-phase friction factor. |
| **GOR Sensitivity** | [`results/horizontal/40f707bc-de97-4d61-a92b-6c91a96744fb.png`](results/horizontal/40f707bc-de97-4d61-a92b-6c91a96744fb.png) | Gas slippage and free-gas velocity in horizontal lateral. |
| **Reservoir Depletion** | [`results/horizontal/73f493bf-5930-4f31-b3db-0d26c0e387bf.png`](results/horizontal/73f493bf-5930-4f31-b3db-0d26c0e387bf.png)<br>[`results/horizontal/5eba661f-78fb-4604-b683-53f1c6bb03c4.png`](results/horizontal/5eba661f-78fb-4604-b683-53f1c6bb03c4.png) | Production rate trajectories across reservoir depletion stages. |
| **Wellhead Pressure ($P_{wh}$)** | [`results/horizontal/b7bbdb84-4171-4a5c-ab6e-72ee50431b76.png`](results/horizontal/b7bbdb84-4171-4a5c-ab6e-72ee50431b76.png) | Impact of manifold/separator backpressure on horizontal well deliverability. |
| **Artificial Gas-Lift** | [`results/horizontal/24dd89ab-3065-413b-a238-58fcfe958e7a.png`](results/horizontal/24dd89ab-3065-413b-a238-58fcfe958e7a.png) | Gas-lift injection response curve for horizontal configuration. |
| **Liquid Loading Diagnostics** | [`results/horizontal/362e0054-6607-4717-bbc4-0b1dad589f73.png`](results/horizontal/362e0054-6607-4717-bbc4-0b1dad589f73.png) | Critical liquid holdup transition in horizontal flow. |
| **Executive 4-Panel Dashboard** | [`results/horizontal/48191d7f-0fe1-48ab-a636-dcfec518c051.png`](results/horizontal/48191d7f-0fe1-48ab-a636-dcfec518c051.png) | Comprehensive horizontal well operational dashboard. |
| **Execution State Screenshot** | [`results/horizontal/Screenshot 2026-08-21 at 12.45.54.png`](results/horizontal/Screenshot%202026-08-21%20at%2012.45.54.png) | Notebook execution capture for horizontal case. |

---

### 3.3 Test Case 3: Directional Well (Schoonebeek-1301 Survey Benchmark)

- **Directory:** [`results/directional_1/`](results/directional_1/)
- **Survey Dataset:** [`NLOG_GS_PUB_SCH-1301.xlsx`](results/directional_1/NLOG_GS_PUB_SCH-1301.xlsx) (Dutch onshore field benchmark well `SCHOONEBEEK-1301`, CDS-FINAL survey).
- **Trajectory Profile:** 64 survey stations, $MD = 0 \to 1,620\text{ m}$ ($5,315\text{ ft}$), max inclination $\alpha_{\max} = 88.83^\circ$ (building from vertical to near-horizontal S-curve).

#### Trajectory & 3D Animation

| Trajectory View | Output Artifact | Description |
| :--- | :--- | :--- |
| **2D Plan & Section** | [`results/directional_1/347a3eb9-7092-4056-83a8-ef29aaad7de9.png`](results/directional_1/347a3eb9-7092-4056-83a8-ef29aaad7de9.png) | Top-down East-North plan view and Vertical Section vs. TVD. |
| **3D Rotating Trajectory** | [`results/directional_1/download.gif`](results/directional_1/download.gif) | Full $360^\circ$ 3D orbital animation of the Schoonebeek wellpath. |

#### Diagnostic Results & Plots

| Analysis | Output Artifact | Key Finding |
| :--- | :--- | :--- |
| **PVT Characterization** | [`results/directional_1/757192b5-ea2f-43a2-a53e-921b9ebb3f04.png`](results/directional_1/757192b5-ea2f-43a2-a53e-921b9ebb3f04.png) | High-resolution 4-panel Black-Oil PVT curves. |
| **Inflow Performance (IPR)** | [`results/directional_1/5147f18d-33a8-49aa-840e-b33a750f1b4c.png`](results/directional_1/5147f18d-33a8-49aa-840e-b33a750f1b4c.png) | Calibrated composite IPR curve. |
| **Pressure Traverse & VLP** | [`results/directional_1/3643c4dc-f7b2-40f4-9498-a9d8a449626a.png`](results/directional_1/3643c4dc-f7b2-40f4-9498-a9d8a449626a.png) | Depth traversal accounting for varying station-by-station inclination. |
| **Base-Case Nodal Solution** | [`results/directional_1/8bfd1c5d-3781-4079-b9ab-4e19fd7cf0a9.png`](results/directional_1/8bfd1c5d-3781-4079-b9ab-4e19fd7cf0a9.png) | Nodal operating point matching directional survey geometry. |
| **Tubing Size Sensitivity** | [`results/directional_1/8f4f29e5-c27d-46c5-a58e-6e044c183b62.png`](results/directional_1/8f4f29e5-c27d-46c5-a58e-6e044c183b62.png) | Sizing comparison along build-and-hold directional profile. |
| **Water Cut Sensitivity** | [`results/directional_1/2d6e2b00-9360-417a-b0f8-d0b93b3ee91a.png`](results/directional_1/2d6e2b00-9360-417a-b0f8-d0b93b3ee91a.png) | Water loading behavior along deviated well profile. |
| **GOR Sensitivity** | [`results/directional_1/18da03be-c71e-4f71-be4a-ca6a5d622f9a.png`](results/directional_1/18da03be-c71e-4f71-be4a-ca6a5d622f9a.png) | GOR impact on hydrostatic lift in deviated section. |
| **Reservoir Depletion** | [`results/directional_1/9609e7af-c5dc-486e-8903-bf04917f6700.png`](results/directional_1/9609e7af-c5dc-486e-8903-bf04917f6700.png)<br>[`results/directional_1/1504f716-66e6-42a9-bd4c-cbd6177c140d.png`](results/directional_1/1504f716-66e6-42a9-bd4c-cbd6177c140d.png) | Production rate and bottomhole pressure depletion response. |
| **Wellhead Pressure ($P_{wh}$)** | [`results/directional_1/d25420d6-731a-4e8d-bfd0-a54836162676.png`](results/directional_1/d25420d6-731a-4e8d-bfd0-a54836162676.png) | Choke sensitivity in directional wellbore. |
| **Artificial Gas-Lift** | [`results/directional_1/6c3075d9-a815-4ad5-a741-88a09ee6d461.png`](results/directional_1/6c3075d9-a815-4ad5-a741-88a09ee6d461.png) | Gas-lift injection optimization along deviated geometry. |
| **Liquid Loading Diagnostics** | [`results/directional_1/11a9faec-1b83-4445-993d-ce50f5f6f319.png`](results/directional_1/11a9faec-1b83-4445-993d-ce50f5f6f319.png) | Critical liquid loading curve for Schoonebeek well. |
| **Executive 4-Panel Dashboard** | [`results/directional_1/73137dae-9aa3-42da-8040-c87200f086e8.png`](results/directional_1/73137dae-9aa3-42da-8040-c87200f086e8.png) | Complete directional well summary dashboard. |
| **Execution State Screenshot** | [`results/directional_1/Screenshot 2026-08-21 at 12.49.17.png`](results/directional_1/Screenshot%202026-08-21%20at%2012.49.17.png) | Execution screenshot showing survey table loading and metrics. |

---

### 3.4 Cross-Case Comparative Synthesis

| Production Scenario | Vertical Well (Case 1) | Horizontal Well (Case 2) | Directional Well (Case 3) | Primary Governing Mechanism |
| :--- | :---: | :---: | :---: | :--- |
| **Base Operating Rate ($q^*$)** | $2,267\text{ STB/d}$ | $2,595\text{ STB/d}$ | $2,180\text{ STB/d}$ | Gravity vs friction trade-off |
| **Flowing BHP ($P_{wf}^*$)** | $2,805\text{ psia}$ | $2,603\text{ psia}$ | $2,858\text{ psia}$ | Required drawdown at node |
| **Default Tubing Size** | $2\text{-}7/8"$ ($2.441"$) | $3\text{-}1/2"$ ($2.992"$) | $2\text{-}7/8"$ ($2.441"$) | Standard completion catalog |
| **Max Production with $4\text{-}1/2"$ Tubing** | $2,828\text{ STB/d}$ | $3,140\text{ STB/d}$ | $2,760\text{ STB/d}$ | Frictional backpressure relief |
| **Rate at $80\%$ Water Cut** | $877\text{ STB/d}$ | $990\text{ STB/d}$ | $810\text{ STB/d}$ | Hydrostatic liquid loading |
| **Optimum Gas-Lift Gain** | $+193\text{ STB/d}$ ($+8.5\%$) | $+220\text{ STB/d}$ ($+8.5\%$) | $+205\text{ STB/d}$ ($+9.4\%$) | Hydrostatic density reduction |
| **Critical Depletion Pressure ($\bar{P}_{\text{crit}}$)** | $2,939\text{ psia}$ | $2,810\text{ psia}$ | $2,980\text{ psia}$ | Minimum VLP point threshold |

---

## 4. Project Structure & Key Files

```
nodal-analysis_well-deliverability/
├── Nodal_Analysis_Petroleum_Production.ipynb  # Primary interactive computational notebook
├── FORMULAS_IN_USE.md                         # Exhaustive mathematical & engineering equations
├── README.md                                  # Complete project documentation & benchmark analysis
├── .gitignore                                 # Git rules (excludes binaries and PPTX artifacts)
└── results/                                   # Benchmark execution runs & diagnostic outputs
    ├── vertical/                              # Case 1: Vertical well outputs & dashboard
    ├── horizontal/                            # Case 2: Horizontal well outputs & dashboard
    └── directional_1/                         # Case 3: Directional well survey & 3D animation
        ├── NLOG_GS_PUB_SCH-1301.xlsx          # Directional survey data (Schoonebeek-1301)
        └── download.gif                       # 3D rotating trajectory animation
```

---

## 5. Installation & Execution Guide

### 5.1 Prerequisites

Ensure you have Python 3.9+ and the required scientific packages installed:

```bash
pip install numpy pandas scipy matplotlib ipywidgets openpyxl pillow
```

To enable interactive ipywidgets in JupyterLab / Jupyter Notebook:

```bash
jupyter nbextension enable --py widgetsnbextension
```

### 5.2 Running the Simulation

1. Launch Jupyter Notebook:
   ```bash
   jupyter notebook Nodal_Analysis_Petroleum_Production.ipynb
   ```
2. Run **Cell 1** (Well Type Selector):
   - Select **Vertical**, **Horizontal**, **Directional**, or **Custom**.
   - For **Directional**, upload an Excel (`.xlsx`) or CSV survey file (e.g. `results/directional_1/NLOG_GS_PUB_SCH-1301.xlsx`), assign MD / INC / AZI columns, and click **Load & Compute TVD**.
3. Run all subsequent cells (`Cell 2` through `Cell 19`) to execute Black-Oil PVT, IPR calibration, VLP multi-phase integration, Nodal solution, sensitivity analysis, liquid loading diagnostics, and the summary executive dashboard.

---

## 6. References & Standards

1. **Beggs, H.D. & Brill, J.P. (1973):** *"A Study of Two-Phase Flow in Inclined Pipes"*, Journal of Petroleum Technology, 25(5), 607–617.
2. **Hagedorn, A.R. & Brown, K.E. (1965):** *"Experimental Study of Pressure Gradients in Two-Phase Vertical Flow"*, JPT, 17(4), 475–484.
3. **Vogel, J.V. (1968):** *"Inflow Performance Relationships for Solution-Gas Drive Wells"*, JPT, 20(1), 83–92.
4. **Klins, M.A. & Clark, C.A. (1993):** *"An Improved Method to Predict Inflow Performance Relationships in Saturated and Saturated/Undersaturated Reservoirs"*, SPE Reservoir Engineering, 8(4), 249–254.
5. **Jones, L.G., Blount, E.M. & Glaze, O.H. (1976):** *"Use of Short-Term Multiple-Rate Flow Tests to Predict Performance of Wells with Turbulence"*, SPE 6133.
6. **Standing, M.B. (1947):** *"A General Pressure-Volume-Temperature Correlation for Constrained Condensates and Crude Oil Systems"*, API Drilling and Production Practice.
7. **Vasquez, M. & Beggs, H.D. (1980):** *"Correlations for Fluid Physical Property Prediction"*, JPT, 32(6), 968–970.
8. **Beggs, H.D. & Robinson, J.R. (1975):** *"Estimating the Viscosity of Crude Oil Systems"*, JPT, 27(9), 1140–1141.
9. **Dranchuk, P.M. & Abou-Kassem, J.H. (1975):** *"Calculation of Z Factors for Natural Gases Using Equations of State"*, Journal of Canadian Petroleum Technology, 14(3), 34–40.
10. **Lee, A.L., Gonzalez, M.H. & Eakin, B.E. (1966):** *"The Viscosity of Natural Gases"*, JPT, 18(8), 997–1000.
11. **API Recommended Practice 11L / Specification 5CT:** *"Specification for Casing and Tubing"*, American Petroleum Institute.

---
*Created for petroleum engineering simulation, completion design, and production optimization research.*
