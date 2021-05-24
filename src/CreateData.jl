using OrdinaryDiffEq, Plots, CSV, DataFrames

#Create ODE Object
function dCdt(du, u, p, t)
    Ka, Ke = p
    du[1] = -Ka*u[1] #Where u[1] = C_GI & u[2] = C_P
    du[2] = Ka*u[1] - Ke*u[2] 
end

plot()
for Ka in 0:0.05:1
    for Ke in 0:0.05:1
        for T in 1:1:12
            dosetimes = 0:T:168 # we will dose every T hours 
            condition(u,t,integrator) = t ∈ dosetimes
            affect!(integrator) = integrator.u[1] += 150 #add this concentration to the GI compartment
            cb = DiscreteCallback(condition,affect!)
            tspan = (0, 168) #Time range
            u0 = [150, 0] #initial concentrations in GI and Plasma compartment respectively
            p = [Ka, Ke]
            prob = ODEProblem(dCdt, u0, tspan, p)
            sol = solve(prob,Rosenbrock23(),callback=cb,tstops=dosetimes)
            plot!(sol.t, sol[2,:]; xlabel = "T", ylabel = "C", legend = false, ticks = false, color = :black, width = 3, xlim = tspan)
        end
    end
end

colNames = []
data = []
for Ka in 0:0.05:1
    for Ke in 0:0.05:1
        for T in 1:1:12
            dosetimes = 0:T:168 # we will dose every T hours 
            condition(u,t,integrator) = t ∈ dosetimes
            affect!(integrator) = integrator.u[1] += 150 #add this concentration to the GI compartment
            cb = DiscreteCallback(condition,affect!)
            tspan = (0, 168) #Time range
            u0 = [150, 0] #initial concentrations in GI and Plasma compartment respectively
            p = [Ka, Ke]
            prob = ODEProblem(dCdt, u0, tspan, p)
            sol = solve(prob,Rosenbrock23(),callback=cb,tstops = 0:1:168)
            push!(colNames, "Ka"*string(Ka)*"_Ke"*string(Ke)*"_T"*string(T))
            push!(data, sol[2,[findall(x -> x == i, sol.t)[end] for i in 0:1:168]]) #select values for 0:1:168 so that the solution sizes are the same. Keep in mind that we save at callbacks twice.
        end
    end
end
data2 = hcat(data...)
df = DataFrame(data2, Symbol.(colNames))
CSV.write("./LinearModelData.csv", df)



