using dn03

initial_theta = pi/4
solution = solve_de_dopri5(matematicno_nihalo, pi/4, 10)

times = solution.t
angles_y = solution[1,:]

using Plots

plot_time_based_angle()

#graf odvisnosti nihajnega časa matematičnega nihala od energije nihala.

compare_harmonic_and_math_plots(pi/1.5)
compare_harmonic_and_math_plots(pi/2)
compare_harmonic_and_math_plots(pi/4)
compare_harmonic_and_math_plots(pi/8)