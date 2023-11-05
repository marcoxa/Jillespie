### -*- Mode: Julia -*-

"""
module Jillespie

An implementation of Gillespie's Stochastic/Montecarlo simulation algorithm.
"""
module Jillespie

using StatsBase
using Distributions


export greet, Model, simulate


"""
greet

The Jillespie module greeting.
"""
greet() = print("Jillespie: a Gillespie simulation algorithm  implementation.")

macro lambda(args, body)
    :($args -> $body)
end


### Model
### =====

"""
Model

A collection of species and reactions.
"""
struct Model
    reactions
    inital_state
    state
    state_changes
    t0

    ## Inner constructor.
    function Model(reactions, initial_state :: Dict{Symbol, Int64})
        new(reactions,
            initial_state,
            copy(initial_state),
            [],
            0.0
            )
    end
end


function Model(reacts, species :: Vector{Symbol})
    Model(reacts,
          Dict{Symbol, Int64}(map(Pair,
                                  species,
                                  zeros(Int64, length(species)))),
          )
end


function Model(initial :: Vector{Tuple{Symbol, Int64}}, reacts)
    Model(reacts, Dict(intial))
end


function Model(initial)
    Model([], initial)
end


### Specie
### ======

function add_specie(model, symbol, initial = 0)
    @assert initial >= 0
    if ! (v in initial_state)  # v must be in model.state
        model.initial_state[symbol] = initial
        model.state[symbol] = initial
    end
end


function set_specie(model, symbol, n)
    @assert n >= 0
    @assert symbol in keys(model.state)
    model.state[symbol] = n
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
    current_state = m.state
    p =  reaction.propensity
          

    p(map(v -> current_state[v], vs))
end


### Simulation engine
### =================

"""
    simulate(model; start, finish)

Simulate a model from `start` to `end`, intended as time.
The function returns a vector of 'times' corresponding to each
reaction occurence, the corresponding vector of `state_changes` from
one state to the next and the model `initial_state`.
The actual state trace at each time step can be reconstructed
afterward (we are not in the 80's anymore, but it is still nice to
save memory).
"""
function simulate(model; start = 0, finish = 1024)
    times = [0.0]
    rates = []
    reactions = model.reactions

    while times[end] <= finish

        rates = map(current_propensity, reactions)

        if all(rate == 0 for rate in rates)
            break
        end

        transition = sample(reactions, rates)
        happen(trasition)

        dt = rand(Exponential(sum(rates)))

        push!(times, times[end] + dt)
    end

    return times, model.state_changes, model.initial_state
end

end # module Jillespie

### end of file Jillespie.jl
