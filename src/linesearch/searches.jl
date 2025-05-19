###     handles the step and stopping

export searchStandard, searchExtSteps, searchExtDist
export SearchStandard, SearchExtendedSteps, SearchExtendedDist, SearchTrueNewton


"""
    searchStandard(p,α,booster,hist,freqs,objFunction,(ϵls,kmax); showtrace=false)

Perform basic linesearch on the booster along `p` with step length `α`. 

# args
- `ϵls`: Minimum required decrease in objective value per step.
- `kmax`: Maximum steps for search.
"""
function searchStandard(p,α,booster,hist,freqs,objFunction,(ϵls,kmax); showtrace=false)
    updateHist!(booster,hist,freqs,objFunction)
    
    k = 0

    while k < kmax
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue > hist[2].objvalue-ϵls
            move(booster,-α*p; additive=true)
            updateHist!(booster,hist,freqs,objFunction)

            break
        end
        
        k += 1
    end

    return k
end

"""
    SearchStandard(ϵls,kmax)

Callback for basic linesearch option in [`linesearch`](@ref). See 
[`searchStandard`](@ref).
"""
const SearchStandard(ϵls,kmax) = Callback(searchStandard,(ϵls,kmax))



"""
    searchExtSteps(p,α,booster,hist,freqs,objFunction,(kmax,); showtrace=false)

Perform forced linesearch on the booster along `p` with step length `α` for `kmax` steps.
Move to best solution on finish.
"""
function searchExtSteps(p,α,booster,hist,freqs,objFunction,(kmax,); showtrace=false)
    updateHist!(booster,hist,freqs,objFunction)
    
    k0 = 0
    obj0 = hist[1].objvalue
    pos0 = copy(booster.pos)

    for k in 1:kmax
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue < obj0
            k0 = k
            obj0 = hist[1].objvalue
            pos0 = copy(booster.pos)
        end
    end

    move(booster,pos0; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return k0
end

"""
    SearchExtendedSteps(kmax)

Callback for extended linesearch option by steps in [`linesearch`](@ref). See 
[`searchExtSteps`](@ref).
"""
const SearchExtendedSteps(kmax) = Callback(searchExtSteps,(kmax,))




"""
    searchExtDist(p,α,booster,hist,freqs,objFunction,(d,); showtrace=false)

Perform forced linesearch on the booster along `p` with step length `α` for distance `d`.
Move to best solution on finish.
"""
function searchExtDist(p,α,booster,hist,freqs,objFunction,(d,); showtrace=false)
    kmax = d/(α*norm(p))

    updateHist!(booster,hist,freqs,objFunction)
    
    k0 = 0
    obj0 = hist[1].objvalue
    pos0 = copy(booster.pos)

    for i in 1:kmax
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue < obj0
            k0 = k
            obj0 = hist[1].objvalue
            pos0 = copy(booster.pos)
        end
    end

    move(booster,pos0; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return k0
end

"""
    SearchExtendedDist(d)

Callback for extended linesearch option by distance in [`linesearch`](@ref). See 
[`searchExtDist`](@ref).
"""
SearchExtendedDist(d) = Callback(searchExtDist,(d,))


"""
    searchTrueNewton(p,α,booster,hist,freqs,objFunction,args; showtrace=false)

Perform a true Newton step (that's probably a bad idea though, but feel free to try
anyways, it's not like I have plenty of experience with these kind of this and that the
booster hast to move anyways, so you could do the linesearch regardless, idiot).
"""
function searchTrueNewton(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    move(booster,α*p; additive=true)
end

"""
    SearchTrueNewton

Callback for true Newton step linesearch option by distance in [`linesearch`](@ref). See 
[`searchTrueNewton`](@ref).
"""
SearchTrueNewton = Callback(searchTrueNewton)
