# 🧊 ME605 - Project 2: Computer Solution of Parabolic PDEs

## 📘 Overview
This repository contains the implementation and report for **Project 2** of the course **ME 605 | Computational Fluid Dynamics**.

The objective of this project is to numerically solve a **2D unsteady heat conduction equation** (a **parabolic Partial Differential Equation**) using the **Finite Difference Method (FDM)** and to analyze its behavior under **Explicit, Implicit, and Crank–Nicolson** time integration schemes.

---

## 🧩 Problem Statement: 2D Transient Heat Conduction
We aim to solve the **2D transient heat conduction equation** (also known as the **Heat Equation**) given by:

$$
\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right)
$$

The equation is solved on a square domain of **$80\ \text{mm} \times 80\ \text{mm}$**, subject to the following conditions:

### **Initial Condition (IC)**
The plate is initially at a uniform temperature:

$$
T(x, y, 0) = 300~\text{K}
$$

### **Boundary Conditions (BCs)**

The boundaries are subjected to Dirichlet (fixed temperature) and Neumann (convective) conditions:

$$
\begin{aligned}
\text{Bottom Edge:} \quad & T(x, 0) = 323~\text{K} \\\\
\text{Top Edge:} \quad & T(x, 80) = 423~\text{K} \\\\
\text{Right Edge:} \quad & T(80, y) = 473~\text{K} \\\\
\text{Left Edge (Convection):} \quad & h \big(T_\infty - T(0, y)\big) = -\lambda \left.\frac{\partial T}{\partial x}\right|_{x=0}
\end{aligned}
$$

### **Given Parameters**
| Parameter | Symbol | Value | Unit |
| :--- | :--- | :--- | :--- |
| Convective heat transfer coefficient | $h$ | 250 | $\text{W/m}^2\text{K}$ |
| Thermal conductivity | $\lambda$ | 5 | $\text{W/mK}$ |
| Free stream temperature | $T_\infty$ | 573 | $\text{K}$ |
| Thermal diffusivity | $\alpha$ | $5 \times 10^{-6}$ | $\text{m}^2/\text{s}$ |

---

## ⚙️ Numerical Implementation
### **Discretization and Domain**
* **Method:** **Finite Difference Method (FDM)** is used for both spatial and temporal discretization.
* **Domain Size:** $80\ \text{mm} \times 80\ \text{mm}$
* **Grid Resolution:** $81 \times 81$ nodes ($\Delta x = \Delta y = 1\ \text{mm}$)
* **Convection BC:** A **first-order approximation** for the Neumann boundary condition is used to preserve the **tridiagonal structure** of the coefficient matrix, which is crucial for efficient implicit solutions.

### **Time Integration Schemes**
The solution is computed using three primary schemes to compare their stability and accuracy:
* **Explicit Method**
* **Implicit Method**
* **Crank–Nicolson Method** (Semi-implicit)

---

## 🧮 Tasks and Analysis

### (1) Explicit Scheme
* Solved the transient heat conduction equation using the **Explicit scheme**.
* **Stability Analysis:** Determined the **critical time step ($\Delta t_{\text{crit}}$)** necessary for ensuring numerical stability.
* Computed the full **temperature evolution** and the final **steady-state** profile.

### (2) Implicit Scheme
* Executed the simulation using a time step $\mathbf{\Delta t_{\text{CFD}} = 2 \Delta t_{\text{crit}}}$ to demonstrate the **unconditional stability** of the implicit approach.
* Followed the same procedure for temperature evolution and steady-state analysis.

### (3) Crank–Nicolson Scheme
* Applied the **semi-implicit time integration** (Crank–Nicolson method).
* Compared transient and steady-state results for accuracy and convergence rate against the other two methods.

### **Solution Approach Details**
* The discretized PDEs are converted into a large **system of linear algebraic equations** at each time step.
* For the implicit and Crank–Nicolson formulations, **Iterative Solvers** (specifically the **Gauss–Seidel** method) were employed to solve the resulting matrix system.

---

## 📊 Results and Visualization
The final report and repository include comprehensive plots illustrating the key findings:

* **Transient Profile Plots:** Temperature variation plotted along the **Centerline $x = 40\ \text{mm}$ (varying $y$)** and the **Centerline $y = 40\ \text{mm}$ (varying $x$)**.
    * Plots include **six equally spaced time instances** from $t = 0$ to $t = t_{\text{end}}$.
* **Contour Plots:** Visualization of the **Steady-State Temperature Distribution** across the entire $80\ \text{mm} \times 80\ \text{mm}$ domain.
* **Method Comparison:** Validation of stability criteria and comparison of **convergence** speed and **steady-state results** among the Explicit, Implicit, and Crank–Nicolson schemes.
