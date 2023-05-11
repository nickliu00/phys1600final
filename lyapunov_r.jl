using Plots 
using JLD
using Statistics 
using LinearAlgebra 
using LinearRegression
using Distances


loaded_data1 = load("planet_h_rs.jld")
loaded_data2 = load("planet_h_rs2.jld")

planet_h_rs = loaded_data1["planet_h_rs"]
planet_h_rs2 = loaded_data2["planet_h_rs2"]

display(planet_h_rs)
display(planet_h_rs2)

# plot(planet_h_thetas)
# plot!(planet_h_thetas2)
# println(size(planet_h_thetas))
# println(size(planet_h_thetas2))

minsize = min(size(planet_h_rs)[1], size(planet_h_rs2)[1])

delta_rs = [euclidean(planet_h_rs[i], planet_h_rs2[i]) for i in 1:minsize]
times = range(1, stop=100, length=minsize)

plot(times, delta_rs, legend=false, title="Δr vs. Time for Planet h from Initial Perturbation", xlabel = "Years", ylabel = "Δθ")

#linear regression to calculate Lyapunov exponent

lr = linregress(times, log.(delta_rs.+0.00001))
println("Slope: "*string(LinearRegression.slope(lr)))
println("Intercept: "*string(LinearRegression.bias(lr)))

Δθ_plot = Plots.plot(times, (delta_rs), xlabel = "t", ylabel = "Δr", legend = false, title = "Δr vs. t")
Δθ_plot_ln = Plots.plot(times, log.(delta_rs.+0.0001), xlabel = "t", ylabel = "ln(Δr)", legend = false, title = "ln(Δr) vs. t")
#Plots.plot(θt, Δθ_plot)
Plots.plot(Δθ_plot_ln)

# p1 = [0.06141494801123557, 0.0, 0.0002090199523157984]
# p2 = [0.06147683765279758, 0.0, 0.00020923058784260478]
# dist = norm(p1-p2)
# println(dist)

