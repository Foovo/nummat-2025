module dn03

using DifferentialEquations
using Plots

export solve_de_dopri5, matematicno_nihalo, harmonicno_nihalo, plot_time_based_angle, compare_harmonic_and_math_plots

# DE system of first order
function matematicno_nihalo(dY, Y, p, t)
    g =  9.81
    l =  10

    theta, z = Y
    dY[1] = z
    dY[2] = -g/l * sin(theta)
end

function harmonicno_nihalo(dY, Y, p, t)
    g =  9.81
    l =  10

    theta, z = Y
    dY[1] = z
    dY[2] = -g/l * theta
end

# Solve differential equation using dopri5 method
function solve_de_dopri5(de, initial_angle, iters)
    Y0 = [initial_angle, 0.0] #initial velocity is 0
    time = (0.0, iters)

    prob = ODEProblem(de, Y0, time)

    return solve(prob, DP5(), saveat=0.01)
end

# Find time, where angle is closest to zero
function find_theta_at_almost_zero(times, theta_values)
    smallest = 0
    smallest_theta = 0
    for (i, theta) in enumerate(theta_values)
        smallest =  i
        smallest_theta = theta
        if theta < 0 
            if abs(theta) < smallest_theta
                return times[i]
            else
                return times[smallest]
            end
        end
    end
    return NaN
end

# Plot graph of pendulum iteration time based on different starting angles
function plot_time_based_angle(de)

    angles = 1:0.001:pi
    time = []

    for angle in angles

        solution = solve_de_dopri5(de, angle, 10)

        times = solution.t
        angles_theta = solution[1,:]

        nihalni_cas = find_theta_at_almost_zero(times, angles_theta)*4
     
        push!(time, nihalni_cas)

    end

    return angles, time
    
end

# Compare plots of the harmonic and mathematical plots based on the starting angle
function compare_harmonic_and_math_plots(angle)
    
    solution = solve_de_dopri5(matematicno_nihalo, angle, 50)

    times = solution.t
    angles_y = solution[1,:]

    plot(times, angles_y, xlabel="t", ylabel="theta(t)", label= "Matematično nihalo", title="Kot nihala čez čas z začetnim kotom $(round(rad2deg(angle), digits=2))°")

    solution = solve_de_dopri5(harmonicno_nihalo, angle, 50)

    times = solution.t
    angles_y = solution[1,:]
    p = plot!(times, angles_y, xlabel="t", ylabel="theta(t)", label= "Harmonično nihalo")

    savefig(p, "graph_$(round(rad2deg(angle), digits=2)).pdf")
    
end


end # module dn03
