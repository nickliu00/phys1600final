using Plots 
using JLD
using Statistics 
using LinearAlgebra 


loaded_data1 = load("planet_h_thetas.jld")
loaded_data2 = load("planet_h_thetas2.jld")

planet_h_thetas = loaded_data1["planet_h_thetas"]
planet_h_thetas2 = loaded_data2["planet_h_thetas2"]

# plot(planet_h_thetas)
# plot!(planet_h_thetas2)
println(size(planet_h_thetas))
println(size(planet_h_thetas2))

delta_thetas = abs.(planet_h_thetas[1:793509] - planet_h_thetas2[1:793509])
times = range(1, stop=100, length=793509)

plot(times, delta_thetas, legend=false, title="Δθ vs. Time for Planet h from Initial Perturbation", xlabel = "Years", ylabel = "Δθ")

# p1 = [0.06141494801123557, 0.0, 0.0002090199523157984]
# p2 = [0.06147683765279758, 0.0, 0.00020923058784260478]
# dist = norm(p1-p2)
# println(dist)

