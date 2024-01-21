# UAV Drone Simulator Project

## Overview

This project encompasses a UAV drone simulator that moves exclusively along two axes, the longitudinal x-axis, and features pitch capability. Developed as part of a course project, the primary objective was to derive the parameters of the matrix describing the dynamic system with minimal variances using a grey-box system. The simulator in Simulink incorporates input buffers and noise to emulate drone noise.

In a subsequent step, a local experiment was conducted, and non-public data, not disclosed here, was extracted. In the second task, the goal is to create a customized signal to excite the system and obtain the covariance matrix with the smallest possible magnitude (Fisher matrix) as an estimator of the error. Remarkably, results indicate that a pseudobinary random signal yields excellent outcomes. Additionally, a robustness analysis was conducted through Monte Carlo simulations.

## Project Structure

- `Simulator/`: Contains the Simulink files for the UAV drone simulator.
- `ExperimentalData/`: Placeholder for non-public experimental data (not included in this public repository).
- `SignalGeneration/`: Code and files related to the creation of a custom signal for system excitation.
- `RobustnessAnalysis/`: Monte Carlo analysis scripts and results.
