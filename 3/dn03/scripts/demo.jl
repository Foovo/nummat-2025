using dn03

initial_theta = pi/4
solution = solve_de_dopri5(matematicno_nihalo, initial_theta, 10)

times = solution.t
angles_y = solution[1,:]


using Plots

angles, time = plot_time_based_angle(matematicno_nihalo)
plot(angles, time, label="Matematično nihalo", xlabel="Nihajni čas", ylabel="theta(t)", title="Nihajni čas glede na začetni kot")

angles, time = plot_time_based_angle(harmonicno_nihalo)
p = plot!(angles, time, label="Harmonično nihalo")
savefig(p, "nihajni_cas_glede_na_kot.pdf")




compare_harmonic_and_math_plots(pi/1.5)
compare_harmonic_and_math_plots(pi/2)
compare_harmonic_and_math_plots(pi/4)
compare_harmonic_and_math_plots(pi/8)