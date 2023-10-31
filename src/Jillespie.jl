### -*- Mode: Julia -*-

"
module Jillespie

An implementation of Gillespie's Stochastic/Montecarlo simulation algorithm.
"
module Jillespie

using StatsBase
using Distributions


export greet, Model, simulate

"The Jillespie module greeting."
greet() = print("Jillespie: a Gillespie simulation algorithm  implementation.")

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
    inital_state
    state
    state_changes
    reactions
    
    function Model(initial_state :: Dict{Symbol, Int64}, reactions)
        new(initial_state,
            copy(initial_state),
            [],
            reactions)
    end
end


function Model(species :: Vector{Symbol}, reacts)
    Model(Dict{Symbol, Int64}(map(Pair,
                                  species,
                                  zeros(Int64, length(species)))),
          reacts)
end


function Model(initial :: Vector{Tuple{Symbol, Int64}}, reacts)
    Model(Dict(intial), reacts)
end


function Model(initial)
    Model(initial, [])
end


### Specie
### ======

function add_specie(model, symbol, initial = 0)
    @assert initial >= 0
    model.species[symbol] = initial
end


function set_specie(model, symbol, n)
    @assert n >= 0
    @assert symbol in keys(model.species)
    model.species[symbol] = n
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


function Reaction(model, variables, stoich, propensity_form :: Expr)
    Reaction(model, variables, (@lambda variables propensity_form), stoich)
end


"
happen(reaction :: Reaction)

Updates the overall state in the model (a slot in `reaction`)
"
function happen(reaction)
    m =  reaction.model
    vs = reaction.reactants
    s =  reaction.stoichiometry

    function reactvar_change(v, s)
        m.state[v] += s[v]
    end

    ## Note: in the following using `broadcast` would be very
    ## cute. However a reaction just involves few reactants; a simple
    ## loop over them might be more efficient.

    ## broadcast(v -> m.state[v] += s[v], vs)

    for v in vs
        m.state[v] += s[v]      # Update the model's state.
    end

    push!(m.state_changes, s)   # Record the change.

    return reaction
end
    

"
current_propensity

Computes the propensity of a reaction in the current state.
"
function current_propensity(reaction)
    m =  reaction.model
    vs = reaction.reactants
    ## s =  reaction.stoichiometry
    p =  reaction.propensity

    p(map(v -> m.state[v], vs))
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

