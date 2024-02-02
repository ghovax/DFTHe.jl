function dft_He(; N::Int64, r_min::Float64, r_max::Float64, prec::Float64, max_iter::Int64)
    # Array to store total energies during iterations
    ens = Float64[]

    # Constants
    step = (r_max - r_min) / (N - 1)
    N_elec = 2
    q_nucl = 2

    # Arrays for radial coordinates, electron density, and potentials
    R = (0:N-1) .* step .+ r_min
    rho, V_nuc, V_hart, V_exch, V_corr, V, U = zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N)

    # Parameters for the correlation potential
    a, b, c, d, gamma, beta1, beta2 = 0.0311, -0.048, 0.002, -0.0116, -0.1423, 1.0529, 0.3334

    # Energy bounds
    E_min, E_max = -20.0, 0.0

    # Iteration variables
    prev_E_tot, E_tot, E_n = 1.0, 0.0, 0.0
    niter = 0

    # Main loop for self-consistent iteration
    while abs(prev_E_tot - E_tot) > prec
        if niter == max_iter
            break
        end
        niter += 1
        prev_E_tot = E_tot

        # Nuclear potential
        V_nuc .= -q_nucl ./ R

        # Hartree potential calculation using finite differences
        U_hart = zeros(N)
        U_hart[1:2] .= 0
        for i in 2:N-1
            U_hart[i+1] = 2.0 * U_hart[i] - U_hart[i-1] - step^2 * rho[i] / R[i]
        end
        alpha = (N_elec - U_hart[N]) / R[N]
        U_hart .+= alpha .* R
        V_hart .= U_hart ./ R

        # Exchange potential
        V_exch .= -((3 * rho) ./ (4 * π^2 * R .^ 2)) .^ (1.0 / 3.0)

        # Correlation potential calculation
        for i in 1:N
            if rho[i] < 1e-10
                V_corr[i] = 0.0
            else
                rs = (3 * R[i]^2 / rho[i])^(1 / 3)
                if rs < 1
                    V_corr[i] = a * log(rs) + b - a / 3 + c * 2 / 3 * rs * log(rs) + (2 * d - c) * rs / 3
                elseif rs < 1e10
                    e_c = gamma / (1 + beta1 * sqrt(rs) + beta2 * rs)
                    V_corr[i] = e_c * (1 + beta1 * 7 / 6 * sqrt(rs) + beta2 * 4 / 3 * rs) / (1 + beta1 * sqrt(rs) + beta2 * rs)
                else
                    V_corr[i] = 0.0
                end
            end
        end

        # Total potential
        V .= V_nuc + V_hart + V_exch + V_corr

        # Electron energy levels calculation using binary search
        n, l, hse_prec = 1, 0, 1.0e-9
        E_max, E_min = 0.0, -20.0
        E_n = 0.0

        while abs(E_max - E_min) > hse_prec
            E_n = (E_min + E_max) / 2.0

            U[N-1:N] .= R[N-1:N] .* exp.(-R[N-1:N])
            for i in N-1:-1:2
                U[i-1] = 2 * U[i] - U[i+1] + step^2 * (-2.0 * E_n + 2 * V[i]) * U[i]
            end

            nodes = sum(U[1:N-1] .* U[2:N] .< 0)
            nodes > n - l - 1 ? (E_max = E_n) : (E_min = E_n)
        end

        # Normalize electron density
        norm = (U[1]^2 + U[N]^2) / 2.0 + sum(U[2:N-1] .^ 2.0)
        U ./= sqrt(norm * step)
        rho .= 2 .* U .^ 2

        # Energy contributions calculation
        E_hart = sum(V_hart .* rho) / 2.0 * step
        E_exch = sum(V_exch .* rho) / 2.0 * step
        E_corr = sum(V_corr .* rho) / 2.0 * step

        # Total energy calculation
        E_tot = 2.0 * E_n - E_hart - (E_exch - E_corr) / 2.0

        # Append total energy to the results array
        append!(ens, E_tot)
    end

    return ens
end

# Example usage
res = @time dft_He(N=4096, r_min=1e-4, r_max=50.0, prec=1e-6, max_iter=10)
println(last(res))