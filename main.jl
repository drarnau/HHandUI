# Load packages
using QuantEcon # For Tauchen AR(1) approximation, linear interpolation

# Parameters calibrated by Krusell et al.
const β = 0.99465       # Discount factor
const α = 0.485         # Utility cost of working
const γ_bar = 0.042     # Average search cost
const ϵ_γ = 0.030       # Standard deviation search cost

const ρ_z   = 0.996       # Persistance productivity process
const σ_ϵ = 0.096       # Standard deviation productivity process

const σ_q = 0.034       # Standard deviation match quality process

const μ = 1.0/6.0       # Average duration
const b_0 = 0.23        # Default replacement ratio
const b_bar = 0.465     # Benefits cap

const λ_e = 0.121       # Probability of finding another job for employed agents
const λ_u = 0.278       # Probability of finding a job for unemployed agents
const λ_n = 0.182       # Probability of finding a jon for OLF agents
const σ = 0.0178        # Probability of losing a job for employed agents

const θ = 0.3           # Capital share of output in the aggregate production function
const δ = 0.0067        # Capital depreciation

const τ = 0.3           # Proportional tax on labor income

# Numerical parameters
const tiny_number = 1e-10
const maxIter = 20000   # Maximum number of iterations when computing policy functions
const tolerance = 1e-4  # Tolerace for value function iteration
const showError = 100   # Frequency of value function iteration error display

const gp_a = 48         # Grid points for assets
const min_a = 0.0       # Minimum level of assets (Krusell et al. = 0.0)
const max_a = 1440.0    # Maximum level of assets (Krussel et al. = 1440.0)

const gp_z = 20         # Grid points productivity process
const cover_z = 2       # Number of standard deviations to each side the productivity process

const gp_q = 7          # Grid points quality match quality process
const cover_q = 2       # Number of standard deviations to each side the match quality process

const gp_γ = 3          # Grid points for serach effort

# Numerical funtions
"""
    loggrid(start, stop, n)

Construct a "logarithmic" grid as implemented in Krusell et al.
"""
function loggrid(start::Real, stop::Real, n::Int64)
    # Check start is positive
    if start < 0
        throw(ArgumentError("start cannot be negative"))
    end
    if stop < 0
        throw(ArgumentError("stop cannot be negative"))
    end
    if stop < start
        throw(ArgumentError("stop cannot be smaller than start"))
    end
    if n < 2
        throw(ArgumentError("n cannot be smaller than 2"))
    end
    step = (log(stop+2.0)-log(start+2.0))/(n-1)
    grid = zeros(n)
    grid[1] = start
    for element in 2:length(grid)
        grid[element] = exp(log(grid[element-1]+2.0) + step) - 2.0
    end
    grid[end] = stop
    return grid
end

# Model-specific types
"""
Stores the Value Functions N, U, W, J, and V for Krussel et al.

An instance of this type admits:
- N(gp_a,gp_z)
- U(gp_a,gp_z,gp_γ,2)
- W(gp_a,gp_z,gp_q)
And computes/creates:
- J(gp_a,gp_z,gp_γ,2)
- V(gp_a,gp_z,gp_q,gp_γ,2)
"""
mutable struct ValueFunctions
    N::Array{Float64,2}
    U::Array{Float64,4}
    W::Array{Float64,3}
    J::Array{Float64,4}
    V::Array{Float64,5}

    function ValueFunctions(N, U, W)
        J = Array{Float64}(gp_a,gp_z,gp_γ,2)
        V = Array{Float64}(gp_a,gp_z,gp_q,gp_γ,2)
        for a in 1:gp_a, z in 1:gp_z, γ in 1:gp_γ, Iᴮ in 1:2
            J[a,z,γ,Iᴮ] = max(U[a,z,γ,Iᴮ],N[a,z])
            for q in 1:gp_q
                V[a,z,q,γ,Iᴮ] = max(W[a,z,q],J[a,z,γ,Iᴮ])
            end
        end
        new(N, U, W, J, V)
    end
end

# Model-specific functions
"""
    utility(consumption[, works, searches, ind_γ])

Compute instantaneous utility as in Krusell et al.

# Arguments
- `consumption::Float64` : amount of consumption
- `works::Bool` : whether the agent works
- `searches:: Bool` : whether the agent seraches
- `ind_γ::Int64` : index of the serach shock
"""
function utility(
            consumption::Float64,
            works::Bool = false,
            searches::Bool = false,
            ind_γ::Int64 = 0
            )
    if searches
        if ind_γ > 0
            aux_search = γ_values[ind_γ]*searches
        else
            throw(ArgumentError("When searches is true,
                                realisation_γ needs to be a valid index for γ"))
        end
        if works
            throw(ArgumentError("Agent cannot work and search simultaneously"))
        end
    else
        aux_search = 0.
    end

    if consumption <= 0.
        return -1e10
    else
        return log(consumption) - (α*works) - aux_search
    end
end

"""
    w, r = prices(KLratio)

Compute aggregate prices given a capital to labor ratio
"""
function prices(KLratio::Float64)
    if KLratio <= 0.0
        throw(ArgumentError("Capital to labor ration needs to be positive!"))
    else
        wage = (1.0-θ)*(KLratio^θ)
        net_r = θ*(KLratio^(θ-1.0)) - δ
        return wage, net_r
    end
end

"""
    benefits(productivity)

Compute UI payments according to the implementation in Krusell et al.

**CAREFUL**: wage already included in the function output
"""
function benefits(productivity::Float64)
    if productivity <= 0.0
        throw(ArgumentError("Productivity needs to be positive!"))
    else
        if (productivity*b_0) < (average_z*b_bar)
            return productivity*b_0*wage
        else
            return average_z*b_bar*wage
        end
    end
end

"""
    value_N(ind_a,ind_z,VFs)

Create a *function* of assets tomorrow to interpolate value for OLF Bellman equation

# Arguments
- `ind_a::Int64`: index for assets today
- `ind_z::Int64`: index for productivity shock today
- `VFs::ValueFunctions`: current values of value functions
"""
function value_N(
            ind_a::Int64,
            ind_z::Int64,
            VFs::ValueFunctions
            )
    # Unpack some variables to improve readability
    a = a_values[ind_a]     # Value of assets today
    V = VFs.V               # Value function V
    J = VFs.J               # Value function J

    # Auxiliary vector to store the expected value for each level of assets tomorrow
    aux_exp = zeros(gp_a)
    for ind_a′ in 1:gp_a
        # sum_prob = 0.   # Auxiliary variable to check the sum of probabilities
        for ind_z′ in 1:gp_z, ind_q′ in 1:gp_q, ind_γ′ in 1:gp_γ
            prob =  z_z′[ind_z,ind_z′]*
                    q_q′[ind_q′]*
                    γ_γ′[ind_γ′]
            aux_exp[ind_a′] = aux_exp[ind_a′] + (prob*
                                                ((λ_n*V[ind_a′,ind_z′,ind_q′,ind_γ′,2])+
                                                ((1.-λ_n)*J[ind_a′,ind_z′,ind_γ′,2])))
            # sum_prob = sum_prob + prob
        end
        # if sum_prob ≉ 1.
        #     error("Probability sum in value_N not equal 1!")
        # end
    end
    # Function of assets tomorrow
    return_f =  function(a′::Float64)
                    # Check assets are in bounds

                    # Define interpolation for assets tomorrow
                    aux_inter = LinInterp(a_values, aux_exp)

                    # Compute consumption
                    c = (1.+r)*a + T - a′

                    return utility(c) + (β*aux_inter(a′))
                end
    return return_f
end

"""
    value_U(ind_a,ind_z,ind_γ,ind_Iᴮ,VFs)

Create a *function* of assets tomorrow to interpolate value for U Bellman equation

# Arguments
- `ind_a::Int64`: index for assets today
- `ind_z::Int64`: index for productivity shock today
- `ind_γ::Int64` : index for serach cost
- `ind_Iᴮ:: Int64`: index for UI indicator today
- `VFs::ValueFunctions`: current values of value functions
"""
function value_U(
            ind_a::Int64,
            ind_z::Int64,
            ind_γ::Int64,
            ind_Iᴮ::Int64,
            VFs::ValueFunctions
            )
    # Unpack some variables to improve readability
    a = a_values[ind_a]             # Value of assets today
    z = z_values[ind_z]             # Value of productivity today
    Iᴮ = Iᴮ_values[ind_Iᴮ]          # Value of indicator UI function (is a Boolean)
    V = VFs.V                       # Value function V
    J = VFs.J                       # Value function J

    # Auxiliary vector to store the expected value for each level of assets tomorrow
    aux_exp = zeros(gp_a)
    for ind_a′ in 1:gp_a
        # sum_prob = 0.   # Auxiliary variable to check the sum of probabilities
        for ind_z′ in 1:gp_z, ind_q′ in 1:gp_q, ind_γ′ in 1:gp_γ, ind_Iᴮ′ in 1:2
            prob =  z_z′[ind_z,ind_z′]*
                    q_q′[ind_q′]*
                    γ_γ′[ind_γ′]*
                    Iᴮ_Iᴮ′[ind_Iᴮ,ind_Iᴮ′]
            aux_exp[ind_a′] = aux_exp[ind_a′] + (prob*
                                                ((λ_u*V[ind_a′,ind_z′,ind_q′,ind_γ′,ind_Iᴮ′])+
                                                ((1.-λ_u)*J[ind_a′,ind_z′,ind_γ′,ind_Iᴮ′])))
            # sum_prob = sum_prob + prob
        end
        # if sum_prob ≉ 1.
        #     error("Probability sum in value_U not equal 1!")
        # end
    end
    # Function of assets tomorrow
    return_f =  function(a′::Float64)
                    # Check assets are in bounds

                    # Define interpolation for assets tomorrow
                    aux_inter = LinInterp(a_values, aux_exp)

                    # Compute consumption
                    c = (1.+r)*a + (1.-τ)*Iᴮ*benefits(z) + T - a′

                    return utility(c,false,true,ind_γ) + (β*aux_inter(a′))
                end
    return return_f
end

"""
    value_W(ind_a,ind_z,ind_q,VFs)

Create a *function* of assets tomorrow to interpolate value for W Bellman equation

# Arguments
- `ind_a::Int64`: index for assets today
- `ind_z::Int64`: index for productivity shock today
- `ind_q::Int64` : index for match qualtiy today
- `VFs::ValueFunctions`: current values of value functions
"""
function value_W(
            ind_a::Int64,
            ind_z::Int64,
            ind_q::Int64,
            VFs::ValueFunctions
            )
    # Unpack some variables to improve readability
    a = a_values[ind_a]             # Value of assets today
    z = z_values[ind_z]             # Value of productivity today
    q = q_values[ind_q]             # Value of match quality
    V = VFs.V                       # Value function V
    J = VFs.J                       # Value function J

    # Auxiliary vector to store the expected value for each level of assets tomorrow
    aux_exp = zeros(gp_a)
    for ind_a′ in 1:gp_a
        # sum_prob = 0.   # Auxiliary variable to check the sum of probabilities
        for ind_z′ in 1:gp_z, ind_q′ in 1:gp_q, ind_γ′ in 1:gp_γ
            prob =  z_z′[ind_z,ind_z′]*
                    q_q′[ind_q′]*
                    γ_γ′[ind_γ′]
            aux_exp[ind_a′] = aux_exp[ind_a′] + (prob*(
                                ((1.-σ-λ_e)*V[ind_a′,ind_z′,ind_q′,ind_γ′,2]) +
                                (λ_e*V[ind_a′,ind_z′,max(ind_q,ind_q′),ind_γ′,2]) +
                                (σ*(1.-λ_u)*J[ind_a′,ind_z′,ind_γ′,1]) +
                                (σ*λ_u*V[ind_a′,ind_z′,ind_q′,ind_γ′,1])))
            # sum_prob = sum_prob + prob
        end
        # if sum_prob ≉ 1.
        #     error("Probability sum in value_W not equal 1!")
        # end
    end
    # Function of assets tomorrow
    return_f =  function(a′::Float64)
                    # Check assets are in bounds

                    # Define interpolation for assets tomorrow
                    aux_inter = LinInterp(a_values, aux_exp)

                    # Compute consumption
                    c = (1.+r)*a + (1.-τ)*wage*q + T - a′

                    return utility(c,true) + (β*aux_inter(a′))
                end
    return return_f
end

# Grids and transition matrices
const a_values = loggrid(min_a, max_a, gp_a)

const z_process = tauchen(gp_z, ρ_z, σ_ϵ, 0.0, cover_z)
const z_values = exp.(z_process.state_values)   # Unpack vector values associated with states
const z_z′ = z_process.p                        # Unpack transition matrix

const q_process = tauchen(gp_q, 0.0, σ_q, 0.0, cover_q)
const q_values = exp.(q_process.state_values)   # Unpack vector values associated with states
const q_q′ = q_process.p[1,:]                   # Unpack probability distribution (ρ = 0.0)

const γ_values = [γ_bar - ϵ_γ, γ_bar, γ_bar + ϵ_γ]          # Vector of search costs
const γ_γ′ = fill(1.0/length(γ_values), length(γ_values))   # Unifrom probability distribution

const Iᴮ_values = [true, false]                 # Vector of UI indicator values
const Iᴮ_Iᴮ′ = [μ 1.0-μ; 0.0 1.0]               # Tranistion matrix UI eligibility


# Initial guesses
KLratio = 129.314057            # Capital to output ratio
(wage, r) = prices(KLratio)     # Wage and net interest rate
average_z = 2.5674              # Average productivity of employed agents
T = 1.40182                     # Government lump-sum transfer to households

N = rand(gp_a,gp_z)            # Value function OLF agent
Nᴺ= ones(gp_a,gp_z)

U = rand(gp_a,gp_z,gp_γ,2)     # Value function Unemployed agent
Uᴺ = ones(gp_a,gp_z,gp_γ,2)

W = rand(gp_a,gp_z,gp_q)       # Value function Employed agent
Wᴺ= ones(gp_a,gp_z,gp_q)

VFs = ValueFunctions(N, U, W)   # Instance with all value functions
VFsᴺ = ValueFunctions(Nᴺ, Uᴺ, Wᴺ)

pf_N = ones(gp_a,gp_z)          # Policy function OLF agent
pf_Nᴺ= zeros(gp_a,gp_z)

pf_U = ones(gp_a,gp_z,gp_γ,2)   # Policy function Unemployed agent
pf_Uᴺ= zeros(gp_a,gp_z,gp_γ,2)

pf_W = ones(gp_a,gp_z,gp_q)     # Policy function Employed agent
pf_Wᴺ= zeros(gp_a,gp_z,gp_q)

Compute policy functions
for iter in 1:maxIter

    for ind_a in 1:gp_a, ind_z in 1:gp_z
        # Unpack some variables to improve readability
        a = a_values[ind_a]             # Value of assets today
        z = z_values[ind_z]             # Value of productivity today

        # Compute upper-boud for assets
        aux_max_a = min(max_a-tiny_number, max(min_a,(1.+r)*a + T))

        # Solve for level of assets tomorrow that maximizes value function
        pf_Nᴺ[ind_a,ind_z], VFsᴺ.N[ind_a,ind_z] =
            golden_method(value_N(ind_a,ind_z,VFs), min_a, aux_max_a)

        for ind_q in 1:gp_q
            # Unpack some variables to improve readability
            q = q_values[ind_q]             # Value of match quality

            # Compute upper-boud for assets
            aux_max_a = min(max_a-tiny_number, max(min_a,(1.+r)*a + (1.-τ)*wage*z*q + T))

            # Solve for level of assets tomorrow that maximizes value function
            pf_Wᴺ[ind_a,ind_z], VFsᴺ.W[ind_a,ind_z] =
                golden_method(value_W(ind_a,ind_z,ind_q,VFs), min_a, aux_max_a)
        end

        for ind_γ in 1:gp_γ, ind_Iᴮ in 1:2
            # Unpack some variables to improve readability
            Iᴮ = Iᴮ_values[ind_Iᴮ]          # Value of indicator UI function (is a Boolean)

            # Compute upper-boud for assets
            aux_max_a = min(max_a-tiny_number, max(min_a,(1.+r)*a + (1.-τ)*benefits(z)*Iᴮ + T))

            # Solve for level of assets tomorrow that maximizes value function
            pf_Uᴺ[ind_a,ind_z], VFsᴺ.U[ind_a,ind_z] =
                golden_method(value_U(ind_a,ind_z,ind_γ,ind_Iᴮ,VFs), min_a, aux_max_a)
        end
    end

    # Compute errors
    error_N = copy(maximum(abs.(VFs.N-VFsᴺ.N)))
    error_U = copy(maximum(abs.(VFs.U-VFsᴺ.U)))
    error_W = copy(maximum(abs.(VFs.W-VFsᴺ.W)))
    error_J = copy(maximum(abs.(VFs.J-VFsᴺ.J)))
    error_V = copy(maximum(abs.(VFs.V-VFsᴺ.V)))

    error_pf_N = copy(maximum(abs.(pf_N-pf_Nᴺ)))
    error_pf_U = copy(maximum(abs.(pf_U-pf_Uᴺ)))
    error_pf_W = copy(maximum(abs.(pf_W-pf_Wᴺ)))

    # Check if convergence is good enoguh
    if  (error_N + error_U + error_W < tolerance) &&
        (error_pf_N + error_pf_U + error_pf_W < tolerance)
        println("Value function iteration finished at iteration $iter")
        break
    elseif rem(iter,showError) == 0
        println("Current value function errors in iteration $iter:")
        println("Error for N is $error_N")
        println("Error for U is $error_U")
        println("Error for W is $error_W")
        println("Error for N's policy function is $error_pf_N")
        println("Error for U's policy function is $error_pf_U")
        println("Error for W's policy function is $error_pf_W")
    end

    # Update value and policy functions
    VFs = ValueFunctions(copy(VFsᴺ.N), copy(VFsᴺ.U), copy(VFsᴺ.W))

    pf_N = copy(pf_Nᴺ)
    pf_U = copy(pf_Uᴺ)
    pf_W = copy(pf_Wᴺ)
end
