
#Copyright (C) 2020 Danielle Rager

using Distributed, Hwloc

addprocs(Hwloc.num_physical_cores()) #adds worker per physical computer core

using Dates, Random, UUIDs, PyPlot, JLD2, HDF5, Statistics, Printf, DataFrames, CSV, Parameters, BenchmarkTools

@everywhere using Distributed, Pkg
@everywhere Pkg.activate(".")
@everywhere using PyPlot, JLD2, HDF5, Statistics, Printf, DataFrames, CSV, Parameters

@everywhere @with_kw struct SimParams
	T::Int64 #total sim length (ms)
end

@everywhere @with_kw struct NCount
	Ne::Int64 #Num E neurons, Rec layer
	Ni::Int64 #Num I neurons, Rec layer
	N0::Int64 #Num FF neurons (all E)
end

@everywhere @with_kw struct ConnProbs
	#connection probabilities. pjk is prob of connection from cell type k to cell type j
	#FF connections
	pe0::Float64
	pi0::Float64
	#Rec connections
	pee::Float64
	pei::Float64
	pie::Float64
	pii::Float64
end


@everywhere @with_kw struct ConnStrength
	strong::Bool
	JR::Float64
	q::Float64
	asymmCoeff::Float64
end

@everywhere @with_kw struct TimeConstants
	# all in ms
	tauerise::Int64
	tauedecay::Int64
	tauirise::Int64
	tauidecay::Int64
	taue::Int64
	taui::Int64
	refrac:: Int64
end

@everywhere @with_kw struct OU
	mulambda0::Float64
	sigma0::Float64
	tau0::Float64
end

@everywhere @with_kw struct Bias
	muemin::Float64
	muemax::Float64
	muimin::Float64
	muimax::Float64
	sigPriv::Float64
end


#simParams = SimParams(12000);

simParams = SimParams(2000);


sysSize = NCount(4000,1000,4000);
connProbs = ConnProbs(0.2,0.5,0.2,0.5,0.5,0.5);
taus = TimeConstants(1,3,1,2,15,10,5);
v4OU = OU(4.0,0.71,60.)


Ncells = sysSize.Ne + sysSize.Ni + sysSize.N0


include("simTwoPopHemiInputWeakRecurrentCoupling.jl")
include("simTwoPopHemiInputNoCoupling.jl")
include("simTwoPopHemiInputUnpack_StrongRec.jl")
include("simTwoPopHemiInput.jl")
include("simTwoPopHemiInput_FreezeConnections_synInput.jl")
include("runSimsFrozenGraph.jl")
include("runSimsFrozenGraphLinResp.jl")
include("currentMethods.jl")

times,ns,times0,ns0,weights,bias,connStrength = simTwoPopHemiInputUnpack_NoCoupleInit(simParams,sysSize,connProbs,taus,v4OU)

times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength, eFFOverTime, eRecOverTime, iRecOverTime = simTwoPopHemiInputUnpack_WeakCoupleInit(simParams,sysSize,connProbs,taus,v4OU)

times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength, numInputs = simTwoPopHemiInputUnpack_StrongRecSymmClusters(simParams,sysSize,connProbs,taus,v4OU)

times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength,v4OU = simTwoPopHemiInputUnpack_StrongRecAsymmClusters(simParams,sysSize,connProbs,taus,v4OU)

println("mean excitatory firing rate: ",mean(1000*ns[1:sysSize.Ne]/simParams.T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(sysSize.Ne+1):(Ncells-sysSize.N0)]/simParams.T)," Hz")

numSims = 50

@time runSimsFrozenGraphLinResp(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)

@time runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)
