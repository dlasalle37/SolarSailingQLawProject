"""
QLawIntegrator: zero-order hold method integrator for QLaw.

Description: Algorithm is essentially as follows:
    1.) Take inputs
    2.) Calculate control angle (alpha/beta)
    3.) integrate for short time
    4.) check stopping criteria
    5.) if not done, go to 2.) and repeat, if done, go to 6.)
    6.) Save output data (States, angles)

    Inputs:
    ps: structure containing everything needed for the propagation (See QLawParams for detailed info)

    Outputs:
    X0: Final states
    exitcode: Symbol to tell the user what caused the exit (success, timed out, etc.)

    data written to files (if ps.writeData=true)
"""
function QLawIntegrator(ps::QLawParams)
    #Pull initial states
    X0 = ps.oe0
    
    #Pull targets

    # Pull other parameters
    integStep = ps.step_size

    # Storage setup
    if ps.writeData
        kep        = Vector{Vector{Float64}}(undef, 0) # 
        sail_angles    = Vector{Vector{Float64}}(undef, 0) # alpha/beta storage
    end
    # begin integration
    done = false
    maxTime = ps.max_sim_time
    abstol = ps.abstol
    reltol = ps.reltol
    idx = 0
    integTime = 0 # time from eph.t0
    exitcode = :None
    while !done
        x = X0;
        # Compute control
        tf = integTime+integStep
        tspan = (integTime, tf)
        α, β, ~ = compute_control(x, ps)
        ps.alpha = α
        ps.beta = β
        u = @SArray [α; β]
        prob = ODEProblem(gauss_variational_eqn!, x, tspan, ps, abstol=abstol, reltol=reltol)
        sol = solve(prob)
        
        # Update variables for next loop
        integTime = tspan[2]
        X0 = sol.u[end]

        # Error Checking
        if sol.retcode != :Success
            println(sol.retcode)
            println(sol.u)
            error(":(")
        end
        
        idx += 1
        ps.current_time += integStep  # stepping forward in time for time-based computations (eg. earth true anom, etc.)

        # add current data for writing
        if ps.writeData
            push!(kep, zeros(6))
            kep[end][1:6] = view(X0, 1:6);
            
            push!(sail_angles, zeros(2))
            sail_angles[end][1:2] = view(u, 1:2)
        end

        # Stopping criteria
        # Compute targeting error
        aerr                = ps.Woe[1]*abs(X0[1] - ps.oet[1]) - ps.oeTols[1]
        eerr                = ps.Woe[2]*abs(X0[2] - ps.oet[2]) - ps.oeTols[2]
        ierr                = ps.Woe[3]*abs(X0[3] - ps.oet[3]) - ps.oeTols[3]
        ωerr                = ps.Woe[4]*abs(acos(cos(X0[4] - ps.oet[4]))) - ps.oeTols[4]
        λerr                = ps.Woe[5]*abs(acos(cos(X0[5] - ps.oet[5]))) - ps.oeTols[5]
        targError           = @SVector [aerr, eerr, ierr, ωerr, λerr]
        if tf >= maxTime
            done=true
            exitcode = :MaxTimeReached
        elseif maximum(targError) <= 0
            done=true
            exitcode=:success
        end
    end

    # Writing
    if ps.writeData
        open(datadir("kep.txt"),   "w") do io; writedlm(io,   kep); end
        open(datadir("angles.txt"),   "w") do io; writedlm(io,   sail_angles); end
    end

    return X0, exitcode
end
