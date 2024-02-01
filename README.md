# DFT implementation for the He atom

This is a Julia implementation of a Density Functional Theory (DFT) solver for the He atom. In this specific implementation, the code computes the total ground-state energy for the He atom.

# Parameters

| Parameter | Description | Default value |
| --- | --- | --- |
| `N` | Number of grid points | 4096 |
| `r_min, r_max` | Minimum and maximum values for the radial coordinate | 1e-4, 50 |
| `prec` | Convergence criterion for iteration | 1e-6 |
| `max_iter` | Maximum number of allowed iterations | 10 |

## `dft_He`

At the core of the program is the function `dft_He`, here are its initialization and methods of operation.

### Initialization of the function

1. Define constants such as π and parameters specific to the He atom (e.g. nuclear charge, potential parameters).
2. Set the default values for the numerical grid (`N`), radial coordinate range (`r_min`, `r_max`), and convergence criteria (`prec`).
2. Create arrays for radial coordinates (`R`), electron density (`rho`), and various potentials (`V_nuc`, `V_hart`, `V_exch`, `V_corr`, `V`, `U`).
4. Set up parameters for the Hartree potential calculation (`U_hart`, `hart_rng`, `alpha`).
5. Set initial energy bounds (`E_min`, `E_max`) and initialize iteration variables (`prev_E_tot`, `E_tot`, `E_n`, `niter`).

### Self-consistent loop

1. Check for convergence by comparing the difference between the previous and current total energies with the specified precision (`prec`).
2. Calculate the nuclear potential contribution (`V_nuc`).
3. Perform a numerical solution for the Hartree potential (`U_hart`) using finite differences.
4. Compute the exchange potential contribution (`V_exch`).
5. Evaluate the correlation potential contribution (`V_corr`) based on the electron density.
6. Compute the total potential (`V`) as the sum of nuclear, Hartree, exchange, and correlation potentials.
7. Solve for electron energy levels using a binary search method to determine the energy eigenvalue (`E_n`).
8. Normalize the electron density (`rho`) and update the potential energy contributions.
9. Calculate the contributions of Hartree, exchange, and correlation potentials to the total energy.
10. Update the total energy (`E_tot`) based on the electron energy levels and potential energy contributions.
11. Append the total energy to the results array (`ens`) for later analysis.
12. Increment the iteration counter (`niter`) and repeat the iteration until convergence or reaching the maximum iteration count.

### Output

Return the array of total energies `ens`.

## Usage
```julia
res = @time dft_He(N=4096, r_min=1e-4, r_max=50.0, prec=1e-6, max_iter=10)
```

## Final considerations

1. This implementation is specific to the Helium atom and may need modifications for other atoms or systems.
2. The code uses numerical methods for solving the Schrödinger equation and self-consistent iteration to converge to the ground state.

Feel free to modify and adapt the code for your specific needs.