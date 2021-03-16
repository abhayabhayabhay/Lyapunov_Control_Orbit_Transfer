# Lyapunov_Control_Orbit_Transfer
Final project for ECH267 - Lyapunov Based Controller for Low-Thrust Orbit Transfer 

This code was developed to simulate a Lyapunov Based Controller for orbit transfer of a low-thrust spacecraft. Two cases are analyzed - an orbit altitude change case and an orbit inclination change case. The user is prompted to input which case to run. The simulation takes a few seconds to 3 minutes to run for the given cases. Different cases may take longer but if the solver takes more than 15 minutes for a simulation on the order of 10,000 time steps, then the solver might not be able to solve the given problem. 

To run the code please make sure to add the Kepler2Carts and PlotEarth directories to the MATLAB path. The directories are given in this repository and can also be found here: 

<https://www.mathworks.com/matlabcentral/fileexchange/80632-kepler2carts>

<https://www.mathworks.com/matlabcentral/fileexchange/25048-plot-earth>

Please see the report within the repository to learn more about the methodology used. 

Thank you for viewing my project! 

## Bibliography

[1] H. Leeghim, D.-H. Cho, S.-J. Jo, and D. Kim, “Generalized Guidance Scheme for Low-Thrust Orbit Transfer,” Mathematical Problems in Engineering, vol. 2014, pp. 1–9, 2014.

[2] D. E. Chang, D. F. Chichka, and J. E. Marsden, “Lyapunov-based transfer between elliptic Keplerian orbits,” Discrete \& Continuous Dynamical Systems - B, vol. 2, no. 1, pp. 57–67, 2002.

[3] Mehrez, M. (2021) 
MPC and MHE implementation in MATLAB using Casadi  <https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi>.
