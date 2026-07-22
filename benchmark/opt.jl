using OpenSOAP
using JuMP
using LinearAlgebra
# using Ipopt, MadNLP
using BenchmarkTools

using DataFrames

using CDDLib, COSMO, Clp, GLPK, HiGHS, Ipopt, JuMP, MadNLP, Pajarito, SCIP, SCS, SDPLR, Tulip

# KNITRO, Ipopt, MadNLP, CDDLib, Clp, COSMO, GLPK, HiGHS, Optim, Pajarito

function run()
    # compare solve time for an LP with different solvers under JuMP
    sim = load_mission("config/mission.yaml")

    P_t = sim.mission.spacecraft.attitude.wheels.torque_env
    P_m = sim.mission.spacecraft.attitude.wheels.momentum_env

    npts = 10^3  # number of points to solve (project onto hull)
    sources = rand(3,npts)
    pts = zeros(3,npts)
    
    # rotate each vector into random orientation (to cover the sphere)
    for i in 1:npts
        sources[:,i] = r_random()*sources[:,i]/LinearAlgebra.norm(sources[:,i]) .* min(P_t.distances...) * 2.0
    end

    optimizers = [CDDLib, COSMO, Clp, GLPK, HiGHS, Ipopt, MadNLP, Pajarito, SCIP, SCS, Tulip]#SDPLR, Tulip]

    modes = [:scalar, :vector]
    
    # results = zeros(length(optimizers),length(modes))

    results = DataFrame(optimizer=String[], scalartime=Float64[], vectortime=Float64[], scalarmem=Float64[], vectormem=Float64[])
    print("\nrunning")
    
    for (o, opt) in enumerate(optimizers)
        # print(opt)
        innerdata = zeros(4)
        k = 1
        for (m, mode) in enumerate(modes)
            print(".")
            # println("running with ", opt)
            # @btime project(x, $P_t, model=JuMP.Model(lopt.Optimizer)) setup=(x=rand(3); lopt=$opt)
            try
                f = @timed for i in 1:npts
                    pts[:,i] = project(sources[:,i], P_t, model=JuMP.Model(opt.Optimizer), mode=mode)
                end
                # println(f)
                # results[o,m] = f.time
                innerdata[k] = f.time
                innerdata[k+2] = f.bytes/1e6 # in MB
            catch
                
                # results[o,m] = Inf
                # println(opt, " failed!")
                innerdata[k] = Inf
                innerdata[k+2] = Inf
            end
            k += 1
        end
        push!(results, (string(nameof(opt)), innerdata...))
    end
    println()
    println(results)
end