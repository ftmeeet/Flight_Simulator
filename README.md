# Beaver DHC-2 Flight Simulator in MATLAB & Simulink

This project is a high-fidelity flight simulator for the **De Havilland Canada DHC-2 Beaver** aircraft, developed entirely within the MATLAB and Simulink environment. It implements a full 12-Degree-of-Freedom (12-DOF) non-linear dynamic model to simulate the aircraft's motion and provides real-time visualization through an integration with the open-source flight simulator, **FlightGear**.

The core of this project was developed as a mini-project for the Aeronautical Engineering program at SVIT, Gujarat Technological University.

## Key Features

-   **High-Fidelity 12-DOF Model:** Implements the complete set of 12 non-linear ordinary differential equations (ODEs) for rigid-body aircraft dynamics.
-   **Data-Driven Simulation:** Utilizes detailed aerodynamic coefficients, mass properties, and engine data specific to the DHC-2 Beaver aircraft, based on established thesis research.
-   **Real-Time Visualization:** Connects seamlessly with FlightGear to provide an immersive, real-time 3D visualization of the aircraft's state (position, orientation, and motion).
-   **Interactive Controls:** The Simulink model allows for real-time input from the user (via sliders for elevator, aileron, and rudder) to control the aircraft during simulation.
-   **Stable Trim Analysis:** The simulation is initialized around a calculated trim state for steady, level flight, and the results demonstrate stable aircraft dynamics.

## Technical Architecture

The simulator is built on a modular architecture that separates the physics engine from the visualization layer:

1.  **MATLAB (`odefun5.m`):** The core physics engine is a MATLAB function that defines the 12 state derivatives. It takes the current state of the aircraft and control inputs, calculates the forces and moments using the DHC-2's aerodynamic data, and returns the time rate of change for each state variable.
2.  **Simulink (`beaver_simulator.slx`):** Simulink acts as the integration hub. It uses an ODE solver (`ode45`) to integrate the state derivatives from the MATLAB function over time. It also houses the user interface blocks for pilot commands and sends the final state data to FlightGear.
3.  **FlightGear:** This open-source simulator acts as the graphical front-end. It receives the aircraft's position (latitude, longitude, altitude) and orientation (phi, theta, psi) from the Simulink Aerospace Blockset and renders the DHC-2 model in a 3D environment.


## The Physics: 12-DOF Equations of Motion

The simulation accurately models the aircraft's flight by solving 12 coupled non-linear ODEs. The state vector `X` represents the complete status of the aircraft at any given time:

`X = [p, q, r, u, v, w, φ, θ, ψ, x, y, z]`

Where:
-   `p, q, r`: Angular rates (Roll, Pitch, Yaw)
-   `u, v, w`: Velocities in the body-fixed frame
-   `φ, θ, ψ`: Euler angles (Roll, Pitch, Yaw)
-   `x, y, z`: Position in the Earth-fixed inertial frame

The simulation uses the Runge-Kutta-based `ode45` solver in MATLAB to numerically integrate these equations.

## Simulation Results

The model was simulated for a 100-second flight duration starting from a stable trim condition. The results show the expected oscillatory but stable behavior of a typical aircraft responding to initial conditions.

## How to Run the Simulation

To run this project, you will need the following dependencies:
-   MATLAB (R2019b or newer)
-   Simulink
-   Simulink Aerospace Blockset
-   FlightGear (version compatible with the Aerospace Blockset)

**Steps:**
1.  Ensure FlightGear is installed and configured to receive data from Simulink.
2.  Clone this repository to your local machine.
3.  Load the aircraft parameters by running `load('dhc2_vars.mat')` in the MATLAB command window.
4.  Open the `beaver_simulator.slx` file in Simulink.
5.  Run the simulation from the Simulink toolbar.
6.  The FlightGear window should launch and display the aircraft simulation in real-time.

## Acknowledgments

-   The aerodynamic data and aircraft parameters for the DHC-2 Beaver were adapted from the MSc thesis: *A SIMULINK environment for flight dynamics and control analysis: Application to the DHC-2 Beaver*.
-   The theoretical foundation for the equations of motion is based on principles from Robert C. Nelson's *Flight Stability and Automatic Control*.
