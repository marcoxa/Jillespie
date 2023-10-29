### -*- Mode: Julia -*-

module Jillespie

using StatsBase
using Distributions


export greet, Model, simulate

"The Jillespie module greeting."
greet() = print("Jillespie: a simple Gillespie implementation.")

macro lambda(args, body)
    :($args -> $body)
end


### Model
### =====

"
Model

A collection of species and reactions.
"
struct Model
    species
    state_changes
    reactions
    function Model(species, sc, reactions)
        new(Dict{Symbol, Int64}(map(Pair,
                                    species,
                                    zeros(Int64, length(species)))),
            sc,
            reactions)
    end
end


Model(species :: Vector{Symbol}) = Model(species, [], [])

Model(initial :: Dict{Symbol, Int64}) = Model(keys(initial), [initial], [])

function Model(initial :: Vector{Tuple{Symbol, Int64}})
    Model(map(first, intial), [Dict(intial)], [])
end


### Reaction
### ========

struct Reaction
    model :: Model
    reactants                   # The reaction "variables"
    propensity                  # A function
    stoichiometry               # A vector of changes

    # Inner constructor prevents default ones from being defined.
    function Reaction(m, r, p, s)
        
        @assert length(r) == length(s)
        
        r = new(m, r, p, s)
        push!(m.reactions, r)
        r
    end
end


function Reaction(model, variables, stoich, form :: Expr)
    Reaction(model, variables, (@lambda variables form), stoich)
end


function happen(reaction)
    m =  reaction.model
    vs = reaction.reactants
    s =  reaction.stoichiometry

    # foreach(
    
end
    


function selectvars(react, state)
    
end


### Simulation engine
### =================

"""
    simulate(model; start, finish)

Simulate a model from `start` to `end`, intended as time.
"""
function simulate(model; start = 0, finish = 1024)
    local times = [0.0]
    local counts = model.initials
    local state = []
    local rates = []
    local reactions = model.reactions

    while times[end] <= finish
        state = counts[end]

        rates = [r.propensity(selectvars(r, state)... for r in reactions] 

        if all(rate == 0 for rate in rates)
            break
        end

        transition = rand(reactions)
        next_state = advance_state(transition, state)

        dt = rand(Exponential(sum(rates)))

        push!(times, times[end] + dt)
        append!(counts, next_state)
    end

    return times, counts
end

end # module Jillespie

### end of file Jillespie.jl

