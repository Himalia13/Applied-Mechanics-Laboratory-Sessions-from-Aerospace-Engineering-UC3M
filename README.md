# 🔬 Dynamic Systems Simulations

This repository contains simulations and analyses of various dynamic systems in the field of physics, with MATLAB implementations and detailed documentation in PDF format.

## 📁 Repository Structure

```
.
|-- Compound Double Pendulum
|   |-- code
|   |   |-- main_exp.m
|   |   |-- main_num.m
|   |   |-- poncaire.mat
|   |   |-- tracked_positions.mat
|   |-- Compound Double Pendulum.pdf
|
|-- Particle connected to a spool
|   |-- code
|   |-- Particle connected to a spool.pdf
|
|-- Particle on oscillating loop
|   |-- code
|   |-- Particle on oscillating loop.pdf
```

## 📂 Project Descriptions

### 1. Compound Double Pendulum
**Description:** Simulation of a compound double pendulum, exploring chaotic trajectories and nonlinear dynamics.

- **Main Files:**
  - `main_exp.m` - Code for double pendulum experiments.
  - `main_num.m` - Numerical simulation of the pendulum.
  - `poncaire.mat` - Poincaré map data.
  - `tracked_positions.mat` - Tracked positions during the simulation.
- **Document:** `Compound Double Pendulum.pdf` - Theoretical explanation and results analysis.

### 2. Particle connected to a spool
**Description:** Simulation of a particle connected to a spool, analyzing the system's dynamics. A heavy particle \(P\) of mass \(m\) is connected by a string of total length \(l\) to a cylindrical spool of radius \(a\) and center at \(O\). The spool rotates at a known constant angular rate \(\omega\) counterclockwise. The string is massless, inelastic, thin, and remains in tension during the motion of the particle. The inertial reference frame \(Oxy\) is defined with the origin at the spool center. 

- **Laboratory Session 1 Tasks:**
  1. Verify that the system has one degree of freedom, parametrized by \(\phi\).
  2. Express \(\xi\) as a function of time \(t\) and angle \(\phi\). Derive position, velocity, and acceleration vectors.
  3. Formulate the forces acting on \(P\).
  4. Write the second-order differential equation for \(\phi\) and derive tension \(T\) as a function of \(\phi\) and its derivative.
  5. Implement the differential equation in MATLAB as a function for integration with `ode45`.
  6. Create a stop function to detect when \(\xi = 0\) or \(T = 0\).
  7. Integrate the equations using `ode45` and plot the evolution of \(\phi\) and \(\xi\) over time.
  8. Plot the trajectory of \(P\) in the \(Oxy\) plane, including a reference circle for the spool.
  9. Compute and plot the mechanical energy \(E\) of \(P\) over time.
 10. Analyze cases of particle interaction with the spool and string tension loss.
 11. Discuss potential code modifications for post-tension loss integration.

- **Parameters:**
  - \(l = 8 \, m\)
  - \(g = 9.81 \, m/s^2\)
  - \(m = 0.1 \, kg\)
  - \(a = 0.2 \, m\)
  - \(\phi(0) = 0 \, rad\)
  - \(\xi(0) = 1 \, m\)

- **Cases:**
  - **Case 1:** \(\omega = 0.1 \, rad/s, \, x'_P(0) = 0 \, m/s\)
  - **Case 2:** \(\omega = 0 \, rad/s, \, x'_P(0) = -1 \, m/s\)
  - **Case 3:** \(\omega = 0 \, rad/s, \, x'_P(0) = -10 \, m/s\)
  - **Case 4:** \(\omega = 0 \, rad/s, \, x'_P(0) = 4 \, m/s\)
  - **Case 5:** \(\omega = 0.1 \, rad/s, \, x'_P(0) = 1 \, m/s\)

- **Document:** `Particle connected to a spool.pdf` - Theoretical description and simulations.
- **Code:** To be added.

### 3. Particle on oscillating loop
**Description:** Analysis of a particle's movement on an oscillating loop, investigating trajectories and dynamic behavior.

- **Document:** `Particle on oscillating loop.pdf` - System analysis and graphics.
- **Code:** To be added.

## 🚀 How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/user/repo
   ```
2. Navigate to the project of interest and open the .m files in MATLAB to run the simulations.

3. Review the PDFs to understand the theoretical framework and the results of each simulation.

## 🛠️ Requirements
- MATLAB (Version 2020 or higher recommended)
- Signal processing and dynamic systems toolbox

## 🤝 Contribution
Contributions are welcome. If you want to add or improve a simulation, create a pull request or open an issue to discuss it.

## 📄 License
This project is licensed under the MIT License. See the LICENSE file for more details.

