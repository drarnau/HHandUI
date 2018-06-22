# Load packages
using QuantEcon # For Tauchen AR(1) approximation, linear interpolation
using BenchmarkTools, Compat # Time benckmarking

# Numerical functions
"""
    loggrid(start, stop, n)

Construct a "logarithmic" grid as implemented in Krusell et al.
"""
function loggrid(start::Float64, stop::Float64, n::Int64)
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
	foo = 6.4
    return grid
end

# Model-specific types
"""
Stores and creates all unchangable fundamenals of the model:
- Calibrated parameters.
- Assigned parameters.
- Numerical parameters.
- Grids and transition matrices.

An instance of this type can be created with default values with:

`myfundamentals = Fundamentals()`

A given parameter(s) can be assigned with:

`myfundamentals = Fundamentals(β = 0.96, gp_a = 101)`
"""
immutable Fundamentals
	# Calibrated
	α::Float64				# Utility cost of working
	β::Float64       		# Discount factor
	γ_bar::Float64 			# Average search cost
	ϵ_γ::Float64			# Standard deviation search cost
	ρ_z::Float64			# Persistance productivity process
	σ_ϵ::Float64			# Standard deviation productivity process
	σ_q::Float64			# Standard deviation match quality process
	λ_e::Float64			# Probability of finding another job for employed agents
	λ_u::Float64			# Probability of finding a job for unemployed agents
	λ_n::Float64			# Probability of finding a jon for OLF agents
	σ::Float64				# Probability of losing a job for employed agents
	# Assigned
	μ::Float64				# Average duration
	b_0::Float64			# Default replacement ratio
	b_bar::Float64			# Benefits cap
	θ::Float64				# Capital share of output in the aggregate production function
	δ::Float64				# Capital depreciation
	τ::Float64				# Proportional tax on labor income
	# Numerical
	gp_a::Int64 			# Grid points for assets
	gp_z::Int64				# Grid points for match productivity process
	gp_q::Int64				# Grid points for match quality process
	gp_γ::Int64				# Grid points for match grid points for search effort
	min_a::Float64 			# Minimum level of assets
	max_a::Float64			# Maximum level of assets
	cover_z::Int64			# Number of standard deviations to each side the productivity process
	cover_q::Int64			# Number of standard deviations to each side the match quality process
	# Grids
	a_values::Vector{Float64} 	# Values for assets
	z_values::Vector{Float64} 	# Values producitivy process
	q_values::Vector{Float64} 	# Values match quality process
	γ_values::Vector{Float64}	# Values search cost process
	Iᴮ_values::Vector{Bool}		# Values UI status
	# Tranistion matrices
	z_z′::Array{Float64,2}		# Productivity process
	q_q′::Vector{Float64} 		# Match quality
	γ_γ′::Vector{Float64} 		# Search costs
	Iᴮ_Iᴮ′::Array{Float64,2}	# UI

	function Fundamentals(;
		α::Float64 = 0.485,
		β::Float64 = 0.99465,
		γ_bar::Float64 = 0.042,
		ϵ_γ::Float64 = 0.030,
		ρ_z::Float64 = 0.996,
		σ_ϵ::Float64 = 0.096,
		σ_q::Float64 = 0.034,
		λ_e::Float64 = 0.121,
		λ_u::Float64 = 0.278,
		λ_n::Float64 = 0.182,
		σ::Float64 = 0.0178,
		μ::Float64 = 1.0/6.0,
		b_0::Float64 = 0.23,
		b_bar::Float64 = 0.465,
		θ::Float64 = 0.3,
		δ::Float64 = 0.0067,
		τ::Float64 = 0.3,
		gp_a::Int64 = 48,
		gp_q::Int64 = 7,
		gp_z::Int64 = 20,
		gp_γ::Int64 = 3,
		min_a::Float64 = 0.0,
		max_a::Float64 = 1440.0,
		cover_z::Int64 = 2,
		cover_q::Int64 = 2
			)
		# Grid for assets
		a_values = loggrid(min_a, max_a, gp_a)

		# Productivity process
		z_process = tauchen(gp_z, ρ_z, σ_ϵ, 0.0, cover_z)
	 	z_values = exp.(z_process.state_values)   	# Unpack vector values associated with states
		z_z′ = z_process.p                        	# Unpack transition matrix

		# Match quality process
		q_process = tauchen(gp_q, 0.0, σ_q, 0.0, cover_q)
		q_values = exp.(q_process.state_values)   	# Unpack vector values associated with states
		q_q′ = q_process.p[1,:]                   	# Unpack probability distribution (ρ = 0.0)

		# Search cost process
		γ_values = [γ_bar - ϵ_γ, γ_bar, γ_bar + ϵ_γ]          # Vector of search costs
		γ_γ′ = fill(1.0/length(γ_values), length(γ_values))   # Unifrom probability distribution

		# UI process
		Iᴮ_values = [true, false]                 # Vector of UI indicator values
		Iᴮ_Iᴮ′ = [μ 1.0-μ; 0.0 1.0]               # Tranistion matrix UI eligibility

		new(α, β, γ_bar, ϵ_γ, ρ_z, σ_ϵ, σ_q, λ_e, λ_u, λ_n, σ,
			μ, b_0, b_bar, θ, δ, τ,
			gp_a, gp_z, gp_q, gp_γ, min_a, max_a, cover_z, cover_q,
			a_values, z_values, q_values, γ_values, Iᴮ_values,
			z_z′, q_q′, γ_γ′, Iᴮ_Iᴮ′
			)
	end
end

"""
Stores and creates model equilibirum objects:
- Prices: wage and net interest rate
- K/L ratio

A new instance requires:
- `θ`: capital share of output
- `δ`: depreciation rate

Can accept also:
- `KLratio`: capital to output ratio
- `average_z`: average productivity of employed agents
- `T`: government lump-sum transfer to households
"""
mutable struct Equilibrium
	r::Float64				# Net interest rate
	w::Float64				# Wage
	KLratio::Float64 		# Capital to output ratio
	average_z::Float64 		# Average productivity of employed agents
	T::Float64				# Government lump-sum transfer to households

	function Equilibrium(θ::Float64, δ::Float64 ;
			KLratio::Float64 = 129.314057,
			average_z::Float64 = 2.5674,
			T::Float64 = 1.40182
			)
		# Prices
		w = (1.0-θ)*(KLratio^θ)
        r = θ*(KLratio^(θ-1.0)) - δ
		new(r, w, KLratio, average_z, T)
	end
end

"""
Stores the Value Functions N, U, W, J, and V.

An instance of this type requires:
- `N(gp_a,gp_z)`: OLF value function
- `U(gp_a,gp_z,gp_γ,2)`: unemployed value function
- `W(gp_a,gp_z,gp_q)`: employed value function
- `gp_a`: number of grid points for asset values
- `gp_z`: number of grid points for productivity process
- `gp_q`: number of grid points for quality match process
- `gp_γ`: number of grid points for search cost process

Computes/creates:
- J(gp_a,gp_z,gp_γ,2)
- V(gp_a,gp_z,gp_q,gp_γ,2)
"""
mutable struct ValueFunctions
    N::Array{Float64,2}
    U::Array{Float64,4}
    W::Array{Float64,3}
    J::Array{Float64,4}
    V::Array{Float64,5}

    function ValueFunctions(
		N::Array{Float64,2},
		U::Array{Float64,4},
		W::Array{Float64,3},
		gp_a::Int64, gp_z::Int64, gp_q::Int64, gp_γ ::Int64
		)
        J = Array{Float64}(2,gp_γ,gp_z,gp_a,)
        V = Array{Float64}(2,gp_γ,gp_q,gp_z,gp_a)
        for a in 1:gp_a, z in 1:gp_z, γ in 1:gp_γ, Iᴮ in 1:2
            J[Iᴮ,γ,z,a] = max(U[Iᴮ,γ,z,a],N[z,a])
            for q in 1:gp_q
                V[Iᴮ,γ,q,z,a] = max(W[q,z,a],J[Iᴮ,γ,z,a])
            end
        end
        new(N, U, W, J, V)
    end
end

"""
Stores the decision rules pf_N, pf_U, pf_W
"""
mutable struct DecisionRules
	pf_N::Array{Float64,2}
    pf_U::Array{Float64,4}
    pf_W::Array{Float64,3}

	function DecisionRules(
		pf_N::Array{Float64,2},
		pf_U::Array{Float64,4},
		pf_W::Array{Float64,3}
		)
		new(pf_N, pf_U, pf_W)
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
utility = function(
			f::Fundamentals,
			consumption::Float64,
			works::Bool = false,
			searches::Bool = false,
			ind_γ::Int64 = 0
			)
	if searches
		if ind_γ > 0
			aux_search = f.γ_values[ind_γ]*searches
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
		return log(consumption) - (f.α*works) - aux_search
	end
end
"""
	benefits(productivity)

Compute UI payments according to the implementation in Krusell et al.

**CAREFUL**: wage already included in the function output
"""
function benefits(f::Fundamentals, e::Equilibrium, productivity::Float64)
	if productivity <= 0.0
		throw(ArgumentError("Productivity needs to be positive!"))
	else
		if (productivity*f.b_0) < (e.average_z*f.b_bar)
			return productivity*f.b_0*e.w
		else
			return e.average_z*f.b_bar*e.w
		end
	end
end

"""
Create a *function* of assets tomorrow to interpolate value for OLF Bellman equation
"""
function value_N(
            ind_a::Int64,
            ind_z::Int64,
			f::Fundamentals,
			e::Equilibrium,
            VFs::ValueFunctions
            )
    # Unpack some variables to improve readability
	a = f.a_values[ind_a]	# Value of assets
    V = VFs.V               # Value function V
    J = VFs.J               # Value function J

    # Auxiliary vector to store the expected value for each level of assets tomorrow
    aux_exp = zeros(f.gp_a)
    @inbounds for ind_a′ in 1:f.gp_a
        # sum_prob = 0.   # Auxiliary variable to check the sum of probabilities
        @inbounds for ind_z′ in 1:f.gp_z, ind_q′ in 1:f.gp_q, ind_γ′ in 1:f.gp_γ
            prob =  f.z_z′[ind_z,ind_z′]*
                    f.q_q′[ind_q′]*
                    f.γ_γ′[ind_γ′]
            aux_exp[ind_a′] = aux_exp[ind_a′] + (prob*
                                                ((f.λ_n*V[2,ind_γ′,ind_q′,ind_z′,ind_a′])+
                                                ((1.-f.λ_n)*J[2,ind_γ′,ind_z′,ind_a′])))
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
					aux_inter = LinInterp(f.a_values, aux_exp)

                    # Compute consumption
                    c = (1.+e.r)*a + e.T - a′

                    return utility(f,c) + (f.β*aux_inter(a′))
                end
    return return_f
end

# function solve_N(f::Fundamentals, e::Equilibrium, VFs::ValueFunctions)
# 	pf_N = zeros(f.gp_a, f.gp_z)
# 	N = zeros(f.gp_a, f.gp_z)
# 	for ind_a in 1:f.gp_a, ind_z in 1:f.gp_z
# 		# Unpack some variables to improve readability
# 		a = f.a_values[ind_a]             # Value of assets today
# 		z = f.z_values[ind_z]             # Value of productivity today
#
# 		# Compute upper-boud for assets
# 		aux_max_a = min(f.max_a-tiny, max(f.min_a,(1.+e.r)*a + e.T))
#
# 		# Get function to optimize
# 		aux_f = value_N(ind_a,ind_z,f,e,VFs)
#
# 		# Solve for level of assets tomorrow that maximizes value function
# 		pf_N[ind_a,ind_z], N[ind_a,ind_z] = golden_method(aux_f, f.min_a, aux_max_a)
# 	end
# 	return pf_N, N
# end

"""
Create a *function* of assets tomorrow to interpolate value for U Bellman equation
"""
function value_U(
            ind_a::Int64,
            ind_z::Int64,
            ind_γ::Int64,
            ind_Iᴮ::Int64,
			f::Fundamentals,
			e::Equilibrium,
            VFs::ValueFunctions
            )
    # Unpack some variables to improve readability
	a = f.a_values[ind_a]			# Value of assets
    z = f.z_values[ind_z]           # Value of productivity today
    Iᴮ = f.Iᴮ_values[ind_Iᴮ]        # Value of indicator UI function (is a Boolean)
    V = VFs.V                       # Value function V
    J = VFs.J                       # Value function J

    # Auxiliary vector to store the expected value for each level of assets tomorrow
    aux_exp = zeros(f.gp_a)
    @inbounds for ind_a′ in 1:f.gp_a
        # sum_prob = 0.   # Auxiliary variable to check the sum of probabilities
        @inbounds for  ind_z′ in 1:f.gp_z, ind_q′ in 1:f.gp_q, ind_γ′ in 1:f.gp_γ, ind_Iᴮ′ in 1:2
            prob =  f.z_z′[ind_z,ind_z′]*
                    f.q_q′[ind_q′]*
                    f.γ_γ′[ind_γ′]*
                    f.Iᴮ_Iᴮ′[ind_Iᴮ,ind_Iᴮ′]
            aux_exp[ind_a′] = aux_exp[ind_a′] + (prob*
                                                ((f.λ_u*V[ind_Iᴮ′,ind_γ′,ind_q′,ind_z′,ind_a′])+
                                                ((1.-f.λ_u)*J[ind_Iᴮ′,ind_γ′,ind_z′,ind_a′])))
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
					aux_inter = LinInterp(f.a_values, aux_exp)
                    # Compute consumption
                    c = (1.+e.r)*a + (1.-f.τ)*Iᴮ*benefits(f,e,z) + e.T - a′

                    return utility(f,c,false,true,ind_γ) + (f.β*aux_inter(a′))
                end
    return return_f
end

# function solve_U(f::Fundamentals, e::Equilibrium, VFs::ValueFunctions)
# 	pf_U = zeros(f.gp_a,f.gp_z,f.gp_γ,2)
# 	U = zeros(f.gp_a,f.gp_z,f.gp_γ,2)
#     for ind_a in 1:f.gp_a, ind_z in 1:f.gp_z, ind_q in 1:f.gp_q, ind_γ in 1:f.gp_γ, ind_Iᴮ in 1:2
#         # Unpack some variables to improve readability
#         a = f.a_values[ind_a]             # Value of assets today
#         z = f.z_values[ind_z]             # Value of productivity today
# 		q = f.q_values[ind_q]             # Value of match quality
#      	Iᴮ = f.Iᴮ_values[ind_Iᴮ]          # Value of indicator UI function (is a Boolean)
#
# 		# Compute upper-boud for assets
#         aux_max_a = min(f.max_a-tiny, max(f.min_a,(1.+e.r)*a + (1.-f.τ)*benefits(f,e,z)*Iᴮ + e.T))
#
# 		# Get function to optimize
# 		aux_f = value_U(ind_a,ind_z,ind_γ,ind_Iᴮ,f,e,VFs)
#
#         # Solve for level of assets tomorrow that maximizes value function
#         pf_U[ind_a,ind_z,ind_γ,ind_Iᴮ], U[ind_a,ind_z,ind_γ,ind_Iᴮ] =
# 			golden_method(aux_f, f.min_a, aux_max_a)
#     end
# 	return pf_U, U
# end

"""
Create a *function* of assets tomorrow to interpolate value for W Bellman equation
"""
function value_W(
            ind_a::Int64,
            ind_z::Int64,
			ind_q::Int64,
			f::Fundamentals,
			e::Equilibrium,
            VFs::ValueFunctions
            )
    # Unpack some variables to improve readability
	a = f.a_values[ind_a]			# Value of assets
    z = f.z_values[ind_z]           # Value of productivity today
	q = f.q_values[ind_q]			# Value of match quality
    V = VFs.V                       # Value function V
    J = VFs.J                       # Value function J

    # Auxiliary vector to store the expected value for each level of assets tomorrow
    aux_exp = zeros(f.gp_a)
    @inbounds for ind_a′ in 1:f.gp_a
        # sum_prob = 0.   # Auxiliary variable to check the sum of probabilities
        @inbounds for ind_z′ in 1:f.gp_z, ind_q′ in 1:f.gp_q, ind_γ′ in 1:f.gp_γ
            prob =  f.z_z′[ind_z,ind_z′]*
                    f.q_q′[ind_q′]*
                    f.γ_γ′[ind_γ′]
            aux_exp[ind_a′] = aux_exp[ind_a′] + (prob*(
                                ((1.-f.σ-f.λ_e)*V[2,ind_γ′,ind_q′,ind_z′,ind_a′]) +
                                (f.λ_e*V[2,ind_γ′,max(ind_q,ind_q′),ind_z′,ind_a′,]) +
                                (f.σ*(1.-f.λ_u)*J[1,ind_γ′,ind_z′,ind_a′]) +
                                (f.σ*f.λ_u*V[1,ind_γ′,ind_q′,ind_z′,ind_a′])))

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
					aux_inter = LinInterp(f.a_values, aux_exp)

                    # Compute consumption
                    c = (1.+e.r)*a + (1.-f.τ)*e.w*q + e.T - a′

                    return utility(f,c,true) + (f.β*aux_inter(a′))
                end
    return return_f
end

# function solve_W(f::Fundamentals, e::Equilibrium, VFs::ValueFunctions)
# 	pf_W = zeros(f.gp_a, f.gp_z, f.gp_q)
# 	W = zeros(f.gp_a, f.gp_z, f.gp_q)
# 	for ind_a in 1:f.gp_a, ind_z in 1:f.gp_z, ind_q in 1:f.gp_q
# 		# Unpack some variables to improve readability
# 		a = f.a_values[ind_a]             # Value of assets today
# 		z = f.z_values[ind_z]             # Value of productivity today
# 		q = f.q_values[ind_q]             # Value of match quality
#
# 		# Compute upper-boud for assets
# 		aux_max_a = min(f.max_a-tiny, max(f.min_a,(1.+e.r)*a + (1.-f.τ)*e.w*z*q + e.T))
#
# 		# Get function to optimize
# 		aux_f = value_W(ind_a,ind_z,ind_q,f,e,VFs)
#
# 		# Solve for level of assets tomorrow that maximizes value function
# 		pf_W[ind_a,ind_z,ind_q], W[ind_a,ind_z,ind_q] = golden_method(aux_f, f.min_a, aux_max_a)
# 	end
# 	return pf_W, W
# end

"""
Computes decision rules
"""
function find_DRs(f::Fundamentals, e::Equilibrium,
		maxIter::Int64, showError::Int64, tiny::Float64, tolerance::Float64)

	# Create a guess of value functions
	VFs = ValueFunctions(
		rand(f.gp_z, f.gp_a),
		rand(2, f.gp_γ, f.gp_z, f.gp_a),
		rand(f.gp_q, f.gp_z, f.gp_a),
		f.gp_a, f.gp_z, f.gp_q, f.gp_γ
		)

	# Create "empty" value functions
	VFsᴺ = ValueFunctions(
		zeros(f.gp_z, f.gp_a),
		zeros(2, f.gp_γ, f.gp_z, f.gp_a),
		zeros(f.gp_q, f.gp_z, f.gp_a),
		f.gp_a, f.gp_z, f.gp_q, f.gp_γ
		)

	# Create a guess for decision rules
	DRs = DecisionRules(
		rand(f.gp_z, f.gp_a),
		rand(2, f.gp_γ, f.gp_z, f.gp_a),
		rand(f.gp_q, f.gp_z, f.gp_a)
		)

	# Create "empty" decision rules
	DRsᴺ = DecisionRules(
		zeros(f.gp_z, f.gp_a),
		zeros(2, f.gp_γ, f.gp_z, f.gp_a),
		zeros(f.gp_q, f.gp_z, f.gp_a)
		)

	# Value function iteration
	@inbounds for iter in 1:maxIter
	    @inbounds for ind_a in 1:f.gp_a, ind_z in 1:f.gp_z
	        # Unpack some variables to improve readability
	        a = f.a_values[ind_a]             # Value of assets today
	        z = f.z_values[ind_z]             # Value of productivity today

	        # Compute upper-boud for assets
	        aux_max_a = min(f.max_a-tiny, max(f.min_a,(1.+e.r)*a + e.T))

	        # Solve for level of assets tomorrow that maximizes value function
	        DRsᴺ.pf_N[ind_z,ind_a], VFsᴺ.N[ind_z,ind_a] =
	            golden_method(value_N(ind_a,ind_z,f,e,VFs), f.min_a, aux_max_a)

	        @inbounds for ind_q in 1:f.gp_q
	            # Unpack some variables to improve readability
	            q = f.q_values[ind_q]             # Value of match quality

	            # Compute upper-boud for assets
	            aux_max_a = min(f.max_a-tiny, max(f.min_a,(1.+e.r)*a + (1.-f.τ)*e.w*z*q + e.T))

	            # Solve for level of assets tomorrow that maximizes value function
	            DRsᴺ.pf_W[ind_q,ind_z,ind_a], VFsᴺ.W[ind_q,ind_z,ind_a] =
	                golden_method(value_W(ind_a,ind_z,ind_q,f,e,VFs), f.min_a, aux_max_a)
	        end

	        @inbounds for ind_γ in 1:f.gp_γ, ind_Iᴮ in 1:2
	            # Unpack some variables to improve readability
	            Iᴮ = f.Iᴮ_values[ind_Iᴮ]          # Value of indicator UI function (is a Boolean)

	            # Compute upper-boud for assets
	            aux_max_a = min(f.max_a-tiny, max(f.min_a,(1.+e.r)*a + (1.-f.τ)*benefits(f,e,z)*Iᴮ + e.T))

	            # Solve for level of assets tomorrow that maximizes value function
	            DRsᴺ.pf_U[ind_Iᴮ,ind_γ,ind_z,ind_a], VFsᴺ.U[ind_Iᴮ,ind_γ,ind_z,ind_a] =
	                golden_method(value_U(ind_a,ind_z,ind_γ,ind_Iᴮ,f,e,VFs), f.min_a, aux_max_a)
	        end
	    end

	    # Compute errors
	    error_N = copy(maximum(abs.(VFs.N-VFsᴺ.N)))
	    error_U = copy(maximum(abs.(VFs.U-VFsᴺ.U)))
	    error_W = copy(maximum(abs.(VFs.W-VFsᴺ.W)))
	    error_J = copy(maximum(abs.(VFs.J-VFsᴺ.J)))
	    error_V = copy(maximum(abs.(VFs.V-VFsᴺ.V)))

	    error_pf_N = copy(maximum(abs.(DRs.pf_N-DRsᴺ.pf_N)))
	    error_pf_U = copy(maximum(abs.(DRs.pf_U-DRsᴺ.pf_U)))
	    error_pf_W = copy(maximum(abs.(DRs.pf_W-DRsᴺ.pf_W)))

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
	    VFs = ValueFunctions(copy(VFsᴺ.N),copy(VFsᴺ.U), copy(VFsᴺ.W),
							f.gp_a, f.gp_z, f.gp_q, f.gp_γ)
		DRs = DecisionRules(copy(DRsᴺ.pf_N), copy(DRsᴺ.pf_U), copy(DRsᴺ.pf_W))
	end

	return VFs, DRs
end

# Define precision parameters
const maxIter = 20000			# Maximum number of iterations value function iteration
const showError = 100			# Frequence to display error in value function iteration
const tiny = 1e-10				# A tiny positive number
const tolerance = 1e-4			# Tolerance value function iteration

# Create an instance of fundamentals with default values
myFund = Fundamentals()

# Create an instance of Equilibrium with guessed/defult values
myEq = Equilibrium(myFund.θ, myFund.δ)

# Create a guess of value functions
VFs = ValueFunctions(
	rand(myFund.gp_z, myFund.gp_a),
	rand(2, myFund.gp_γ, myFund.gp_z, myFund.gp_a),
	rand(myFund.gp_q, myFund.gp_z, myFund.gp_a),
	myFund.gp_a, myFund.gp_z, myFund.gp_q, myFund.gp_γ
	)

# Compute value functions and decision rules
@time VFs, DRs = find_DRs(myFund, myEq, maxIter, showError, tiny, tolerance)
