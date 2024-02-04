# DFT implementation for the He atom

This Julia implementation is a numerical solver for the DFT, focusing on calculating the total ground-state energy for the He atom. The included [report](dft_HeReport.pdf) explains more in detail what the calculations are. 

## Parameters

| Parameter | Description | Default value |
| --- | --- | --- |
| `N` | Number of grid points | 4096 |
| `r_min, r_max` | Minimum and maximum values for the radial coordinate | 1e-4, 50 |
| `ε` | Convergence criterion for iteration | 1e-6 |
| `n_max` | Maximum number of allowed iterations | 30 |

## `dft_He` function

The `dft_He` function is at the core of the program, encompassing initialization and operational methods for the DFT calculation. Here is its mode of operation.

### Initialization

1. Define constants such as π and parameters specific to the He atom (e.g., nuclear charge, potential parameters).
2. Set default values for the numerical grid (`N`), radial coordinate range (`r_min`, `r_max`), and convergence criteria (`ε`).
3. Create arrays for radial coordinates (`r`), electron density (`rho_gs`), and various potentials (`V_nuc`, `V_H`, `V_x`, `V_c`, `V_eff`, `phi_0`).
4. Set up parameters for the correlation potential (`a`, `b`, `c`, `d`, `γ`, `β_1`, `β_2`).
5. Initialize iteration variables (`E_nm1`, `E_n`, `n`), and set initial values for effective potential (`E_eff0`).

### Self-consistent loop

1. Check for convergence by comparing the difference between the previous and current total energies with the specified precision (`ε`).
2. Calculate contributions from nuclear potential (`V_nuc`), Hartree potential (`V_H`), exchange potential (`V_x`), and correlation potential (`V_c`).
3. Compute the total potential (`V_eff`) as the sum of nuclear, Hartree, exchange, and correlation potentials.
4. Solve for electron energy levels using a binary search method to determine the energy eigenvalue (`E_n`).
5. Normalize the electron density (`rho_gs`) and update potential energy contributions.
6. Calculate contributions of Hartree, exchange, and correlation potentials to the total energy.
7. Update the total energy (`E_n`) based on the electron energy levels and potential energy contributions.
8. Append the total energy to the results array (`E_ns`) for later analysis.
9. Continue the iteration until convergence or reaching the maximum iteration count.

### Output

Return the array of total energies `E_ns`.

## Example of usage
```julia
E_ns = @time dft_He(N=4096, r_min=1e-4, r_max=50.0, ε=1e-6, n_max=30)
```

## Final remarks

This implementation is specific to the He atom and may require adjustments for other atoms or systems, such as number of electrons changes.

Feel free to modify and adapt the code for your specific needs.