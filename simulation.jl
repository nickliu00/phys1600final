#!/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia

# Simulate a solar system

using Plots # for plotting trajectory
using DifferentialEquations # for solving ODEs
using LinearAlgebra # for dot and cross products
using BenchmarkTools
using JLD

G = 4.0*pi^2 # time scale = year and length scale = AU and mass scale = solar mass

mutable struct body
   name::String # name of star or planet
   m::Float64 # mass
   r::Vector{Float64} # position vector
   v::Vector{Float64} # velocity vector
end 

function star(name="Sun", m = 1.0, r = zeros(3), v = zeros(3))
   return body(name, m, r, v)
end

function planet(name="Earth", m=3.0e-6, a=1.0, ϵ=0.017, i=0.0)
   perihelion = (1.0 - ϵ) * a
	aphelion = (1.0 + ϵ) * a
	speed = sqrt(G * 0.089 * (1.0 + ϵ)^2 / (a * (1.0 - ϵ^2)))
   println(speed)

   phi = 0.0
	r = [perihelion * cos(i*pi/180.0) * cos(phi), perihelion * cos(i*pi/180.0) * sin(phi), perihelion * sin(i*pi/180.0)]
	v = [-speed * sin(phi), speed * cos(phi), 0.0]
   return body(name, m, r, v)
end

function momentum(b::body)
   return b.m * b.v
end

function angularMomentum(b::body)
   return b.m * cross(b.r, b.v)
end

function kineticEnergy(b::body)
   v = b.v
   m = b.m
   return 0.5 * m * dot(v, v)
end

function potentialEnergy(body1::body, body2::body)
   r = body1.r - body2.r
   rSquared = dot(r, r)
   return -G * body1.m * body2.m / sqrt(rSquared)
end

function rv(b::body)
   return [b.r; b.v]
end

function force(body1::body, body2::body)
   r = body1.r - body2.r
   rSquared = dot(r, r)
   return -G * r * rSquared^(-1.5)
   #return -G * r / (rSquared * sqrt(rSquared))
end

mutable struct SolarSystem
   bodies::Vector{body}
   numberOfBodies::Int64
   phaseSpace::Matrix{Float64} # 6N dimensional phase space
end

function TotalMass(s::SolarSystem)
   M = 0.0
   for b in s.bodies
      M += b.m
   end
   return M
end

function TotalLinearMomentum(s::SolarSystem)
   P = zeros(3)
   for b in s.bodies
      P += momentum(b)
   end
   return P
end

function TotalAngularMomentum(s::SolarSystem)
   L = zeros(3)
   for b in s.bodies
      L += angularMomentum(b)
   end

   return L
end

function TotalEnergy(s::SolarSystem)
   ke = 0.0
   pe = 0.0

   for body1 in s.bodies
      ke += kineticEnergy(body1)
      for body2 in s.bodies
         if (body1 != body2)
            pe += 0.5 * potentialEnergy(body1, body2)
         end
      end
   end

   return pe + ke
end

function CenterOfMassFrame!(s::SolarSystem)

   M = TotalMass(s)
   P = TotalLinearMomentum(s)
   V = P / M
   
   s.phaseSpace = zeros(6, 0)
   for body in s.bodies
      body.v -= V #boost to COM frame
      s.phaseSpace = [s.phaseSpace rv(body)]
   end

   return nothing
end

function tendency!(dps, ps, s::SolarSystem, t)

   i = 1 # update phase space of individual bodies
   for b in s.bodies
      b.r = ps[1:3, i]
      b.v = ps[4:6, i]
      i += 1
   end

   # find velocities of bodies and forces on them. O(N^2) computational cost
   N = s.numberOfBodies
   for i in 1:N
      b1 = s.bodies[i]
      dps[1:3, i] = b1.v  #dr/dt
      dps[4:6, i] = zeros(3)
      for j in 1:i-1
            b2 = s.bodies[j]
            f = force(b1, b2) # call only once
            dps[4:6, i] += b2.m * f
            dps[4:6, j] -= b1.m * f # Newton's 3rd law
      end
   end

   return nothing
end

function parallel_tendency!(dps, ps, s::SolarSystem, t)

   i = 1 # update phase space of individual bodies
   for b in s.bodies
      b.r = ps[1:3, i]
      b.v = ps[4:6, i]
      i += 1
   end

   # find velocities of bodies and forces on them. O(N^2) computational cost
   N = s.numberOfBodies
   Threads.@threads for i in 1:N
      b1 = s.bodies[i]
      dps[1:3, i] = b1.v  #dr/dt
      dps[4:6, i] = zeros(3)
      for j in 1:N
            if (j != i)
               b2 = s.bodies[j]
               f = force(b1, b2)
               dps[4:6, i] += b2.m * f
            end
      end
   end

   return nothing
end

function eccentricity(body::Int64, solution)
   n = size(solution[1, 1, :])[1]
   samples = 1000
   interval = floor(Int, n/samples)
   ϵ = zeros(samples-1)
   time = zeros(samples-1)
   for i in 1:samples-1
      r2 = sol[1, body, i*interval:(i+1)*interval].^2 + sol[2, body, i*interval:(i+1)*interval].^2 + sol[3, body, i*interval:(i+1)*interval].^2
      aphelion = sqrt(maximum(r2))
      perihelion = sqrt(minimum(r2))
      ϵ[i] = (aphelion - perihelion) / (aphelion + perihelion)
      time[i] = solution.t[i*interval]
   end
   return (time, ϵ)
end

function obliquity(body::Int64, solution)
   n = size(solution[1, 1, :])[1]
   samples = 1000
   interval = floor(Int, n/samples)
   ob = zeros(samples)
   time = zeros(samples)
   for i in 1:samples
      r = solution[1:3, body, i*interval]
      v = solution[4:6, body, i*interval]
      ell = cross(r, v)
      norm = sqrt(dot(ell, ell))
      ob[i] = (180.0 / pi) * acos(ell[3]/norm)
      time[i] = solution.t[i*interval]
   end
   return (time, ob)
end

function SolarSystem()

   bodies = Vector{body}()
   push!(bodies, star("TRAPPIST-1", 1*0.089, zeros(3), zeros(3)))
   push!(bodies, planet("planet_b", 5*4.1354313e-6, 0.01154, 0.00622, 90-89.728))
   push!(bodies, planet("planet_c", 5*3.9354315e-6, 0.01580, 0.00654, 90-89.778))
   push!(bodies, planet("planet_d", 5*1.1666655e-6, 0.02227, 0.00837, 90-89.896))
   push!(bodies, planet("planet_e", 5*2.0816796e-6, 0.02925, 0.00510, 90-89.793))
   push!(bodies, planet("planet_f", 5*3.1264233e-6, 0.03849, 0.01007, 90-89.740))
   push!(bodies, planet("planet_g", 5*3.9753714e-6, 0.04683, 0.00208, 90-89.742))
   push!(bodies, planet("planet_h", 5*9.792783e-7, 0.06189, 0.00667, 90-89.805))

   numberOfBodies = size(bodies)[1]

   phaseSpace = zeros(6, 0)
   for b in bodies
      phaseSpace = [phaseSpace rv(b)]
   end
   
   return SolarSystem(bodies, numberOfBodies, phaseSpace)

end

s = SolarSystem()
CenterOfMassFrame!(s)
println(typeof(s))
println("Initial total energy = ", TotalEnergy(s))
println("Initial total linear momentum = ", TotalLinearMomentum(s))
println("Initial total angular momentum = ", TotalAngularMomentum(s))
println("Number of bodies = ", s.numberOfBodies)
for b in s.bodies
   println("body name = ", b.name)
end

t_final = 100.0 # final time of simulation
tspan = (0.0, t_final) # span of time to simulate
prob = ODEProblem(parallel_tendency!, s.phaseSpace, tspan, s) # specify ODE
sol = solve(prob, progress=true, maxiters=1e8, reltol=1e-6, abstol=1e-6) # solve using Tsit5 algorithm to specified accuracy
println("\n\t Results")
println("Final time  = ", sol.t[end])
println("Final total energy = ", TotalEnergy(s))
println("Final total linear momentum = ", TotalLinearMomentum(s))
println("Final total angular momentum  = ", TotalAngularMomentum(s))

trappist1 = 1
planet_b = 2
planet_c = 3
planet_d = 4
planet_e = 5
planet_f = 6
planet_g = 7
planet_h = 8

planet_h_thetas = atan.(sol[1, planet_h, :], sol[2, planet_h, :])

for (index, theta) in enumerate(planet_h_thetas)
   if index > 1 && planet_h_thetas[index]-planet_h_thetas[index-1] < 0 
      planet_h_thetas[index:end].+=2*pi
   end 
end

# plot(sol.t, planet_h_thetas)
save("planet_h_thetas2.jld", "planet_h_thetas2", planet_h_thetas)
# loaded_data = load("planet_h_thetas.jld")
# loaded_array = loaded_data["planet_h_thetas"]
# plot(sol.t, loaded_array)

# Plot of orbit
xy = plot(
[(sol[1, trappist1, :], sol[2, trappist1, :]), 
 (sol[1, planet_b, :], sol[2, planet_b, :]), 
 (sol[1, planet_c, :], sol[2, planet_c, :]),
 (sol[1, planet_d, :], sol[2, planet_d, :]), 
 (sol[1, planet_e, :], sol[2, planet_e, :]),
 (sol[1, planet_f, :], sol[2, planet_f, :]),
 (sol[1, planet_g, :], sol[2, planet_g, :]),
 (sol[1, planet_h, :], sol[2, planet_h, :])],
 xlabel = "x (AU)", ylabel = "y (AU)", legend = true, title = "TRAPPIST-1 Planetary System", 
 label = permutedims([body.name for body in s.bodies]), aspect_ratio=1, linewidth=0.5)

plot(xy, xlims=(-0.2, 0.2), ylims=(-0.10, 0.10))

# # Plot of obliquity
# tilt = obliquity(planet_b, sol)
# obliquityPlot = plot(tilt, xlabel = "t", ylabel = "tilt", legend = false, title = "obliquity")

# # Plot of eccentricity
# eccentric = eccentricity(planet_b, sol)
# eccentricityPlot = plot(eccentric, xlabel = "t", ylabel = "ϵ", legend = false, title = "eccentricity")

# plot(obliquityPlot, eccentricityPlot)


