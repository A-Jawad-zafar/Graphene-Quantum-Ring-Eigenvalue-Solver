# Graphene-Quantum-Ring-Eigenvalue-Solver
This project contains a MATLAB script that numerically solves the eigenvalue equation of a graphene quantum ring using the Dirac effective model.
# Graphene Quantum Ring Eigenvalue Solver

This repository contains MATLAB scripts that numerically solve the eigenvalue equation of a graphene quantum ring using the Dirac effective model. It focuses on solutions exhibiting physical significance by filtering out those with asymptotic behavior.

## Description

The main script, `GrapheneQuantumRingSolver.m`, employs the bisection method to find solutions where the function changes its sign and can handle different radii values. The purpose is to retain only the physical solutions of the equation and provide insights into the quantum properties of graphene rings.

## Repository Structure

- `GrapheneQuantumRingSolver.m`: The main MATLAB script that initializes the parameters and runs the simulation.
- `Solution.m`: Contains function definitions and performs calculations to find the energy levels near the K-Valley of the Graphene Quantum Ring.
- `asymptotes.m`: A helper function to find the asymptotes in the positive and negative regions.
- `zero_crossings.m`: Computes all values of zero crossing of the Eigenvalue Equation at given input parameters.

## Usage

1. **Clone or Download this Repository**
   - Click on the 'Code' button on the repository page and select "Download ZIP", or clone the repository using Git.
   
2. **Navigate to the Directory**
   - Open MATLAB and navigate to the directory containing the downloaded or cloned files.
   
3. **Run the Main Script**
   - Execute the main script `GrapheneQuantumRingSolver.m` in MATLAB.

## Parameters

- `h`: Reduced Planck constant (in eV·fs)
- `gamma`: The nearest neighbor hopping integral
- `Rin`: Inner radius (in nm)
- `Rout`: Outer radius (in nm)
- `m`: Magnetic quantum number
- `T`: Total Number of input points for Eigen function
- `S`: Number of solutions the code will go through

## References

This code is associated with the following publication:


-  Ahmal Jawad Zafar,, Aranyo Mitra, and Vadym Apalkov. "Ultrafast valley polarization of graphene nanorings." Physical Review B 106.15 (2022): 155147. https://doi.org/10.1103/PhysRevB.106.155147


## Acknowledgements

Major funding was provided by Grant No. DE-FG02- 01ER15213 from the Chemical Sciences, Biosciences, and Geosciences Division, Office of Basic Energy Sciences, Of- fice of Science, U.S. Department of Energy. Numerical simulations were performed using support from Grant No. DE-SC0007043 from the Materials Sciences and Engineering Division of the Office of the Basic Energy Sciences, Office of Science, U.S. Department of Energy.

