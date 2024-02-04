function dft_He(; N::Int64, r_min::Float64, r_max::Float64, ε::Float64, n_max::Int64)
    # Array to store total energies during iterations
    E_ns = Float64[]

    # Constants
    step = (r_max - r_min) / (N - 1)

    # Arrays for radial coordinates, electron density, and potentials
    r = (0:N-1) .* step .+ r_min
    rho_gs = zeros(N)
    V_nuc, V_H, V_x, V_c, V_eff = zeros(N), zeros(N), zeros(N), zeros(N), zeros(N)
    phi_0 = zeros(N)

    # Parameters for the correlation potential
    a, b, c, d, γ, β_1, β_2 = 0.0311, -0.048, 0.002, -0.0116, -0.1423, 1.0529, 0.3334

    # Iteration variables
    E_nm1, E_n = 1.0, 0.0
    n = 0

    E_eff0 = 0.0

    # Main loop for self-consistent iteration
    while abs(E_nm1 - E_n) > ε
        if n == n_max
            break
        end
        n += 1
        E_nm1 = E_n

        # Nuclear potential
        V_nuc .= -2 ./ r

        # Hartree potential calculation using finite differences
        U_H = zeros(N)
        U_H[1:2] .= 0
        for i in 2:N-1
            U_H[i+1] = 2.0 * U_H[i] - U_H[i-1] - step^2 * rho_gs[i] / r[i]
        end
        α = (2 - U_H[N]) / r[N]
        V_H .= U_H ./ r .+ α

        # Exchange potential
        V_x .= -((3 * rho_gs) ./ (4 * π^2 * r .^ 2)) .^ (1.0 / 3.0)

        # Correlation potential calculation
        for i in 1:N
            if rho_gs[i] < 1e-10
                V_c[i] = 0.0
            else
                r_s = (3 * r[i]^2 / rho_gs[i])^(1 / 3)
                if r_s < 1
                    V_c[i] = a * log(r_s) + b - a / 3 + c * 2 / 3 * r_s * log(r_s) + (2 * d - c) * r_s / 3
                elseif r_s < 1e10
                    e_c = γ / (1 + β_1 * sqrt(r_s) + β_2 * r_s)
                    V_c[i] = e_c * (1 + β_1 * 7 / 6 * sqrt(r_s) + β_2 * 4 / 3 * r_s) / (1 + β_1 * sqrt(r_s) + β_2 * r_s)
                else
                    V_c[i] = 0.0
                end
            end
        end

        # Total potential
        V_eff .= V_nuc + V_H + V_x + V_c

        # Electron energy levels calculation using binary search
        n, l = 1, 0

        hse_ε = 10e-9
        E_max, E_min = 0.0, -20.0
        E_eff0 = 0.0

        while abs(E_max - E_min) > hse_ε
            E_eff0 = (E_min + E_max) / 2.0

            phi_0[N-1:N] .= r[N-1:N] .* exp.(-r[N-1:N])
            for i in N-1:-1:2
                phi_0[i-1] = 2 * phi_0[i] - phi_0[i+1] + step^2 * (-2.0 * E_eff0 + 2 * V_eff[i]) * phi_0[i]
            end

            nodes = sum(phi_0[1:N-1] .* phi_0[2:N] .< 0)
            nodes > n - l - 1 ? (E_max = E_eff0) : (E_min = E_eff0)
        end

        # Normalize electron density
        normphi_0 = (phi_0[1]^2 + phi_0[N]^2) / 2.0 + sum(phi_0[2:N-1] .^ 2.0)
        phi_0 ./= sqrt(normphi_0 * step)
        rho_gs .= 2 .* phi_0 .^ 2

        # Energy contributions calculation
        E_H = sum(V_H .* rho_gs) / 2.0 * step
        E_x = sum(V_x .* rho_gs) / 2.0 * step
        E_c = sum(V_c .* rho_gs) / 2.0 * step

        # Total energy calculation
        E_n = 2.0 * E_eff0 - E_H - (E_x - E_c) / 2.0

        # Append total energy to the results array
        append!(E_ns, E_n)
    end

    return E_ns
end

# Example usage
E_ns = @time dft_He(N=4096, r_min=1e-4, r_max=50.0, ε=1e-6, n_max=30)
println("E_gs = ", last(E_ns), " E_h")