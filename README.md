# UAV Drone Simulator Project

## Overview

This project encompasses a UAV drone simulator that moves exclusively along two axes, the longitudinal x-axis, and features pitch capability. Developed as part of a course project, the primary objective was to derive the parameters of the matrix describing the dynamic system with minimal variances using a grey-box system. The simulator in Simulink incorporates input buffers and noise to emulate drone noise.

In a subsequent step, a local experiment was conducted, and non-public data, not disclosed here, was extracted. In the second task, the goal is to create a customized signal to excite the system and obtain the covariance matrix with the smallest possible magnitude (Fisher matrix) as an estimator of the error. Remarkably, results indicate that a pseudobinary random signal yields excellent outcomes. Additionally, a robustness analysis was conducted through Monte Carlo simulations.

## Project Structure
Task 1 is the model identification
Task 2 is the best signal design
- `Task 1/common/`: Contains all the pre implemented function to run the simulation on the quadrotor UAV.
- `Task 1/Main.m`: Main script to identify the model parameter of the UAV with error estimator and Bode Plots.
- `Task 1/Simulator_Single_Axis_excitation.slx`: Simulink model of the UAV.
- `Task 2/RBS/Main_RBS.m`: The code simulates a UAV drone's behavior using a Random Binary Sequence (RBS) as an input signal. It conducts simulations for various RBS sequences, estimates model parameters through a grey-box approach, and assesses the sensitivity of the results to variations in input parameters. The code then generates plots illustrating the covariance matrix trace, simulation trends, and relative errors compared to the actual system parameters.
- `Task 2/RBS/Main_RBS.m`:
- 
## Acknowledgments
Thanks to my professor Marco Lovera and other authors who have contributed to the simulator (Salvatore Meraglia, Mattia Giurato, Paolo Gattazzo).
A special thanks to my colleagues Alessandro Boldrini e Gianluca Napoletano for developing the project with me.
