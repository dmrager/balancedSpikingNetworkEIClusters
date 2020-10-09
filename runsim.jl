
#this file is part of litwin-kumar_doiron_cluster_2012
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


simParams = SimParams(12000);
sysSize = NCount(4000,1000,4000);
connProbs = ConnProbs(0.2,0.5,0.2,0.5,0.5,0.5);
taus = TimeConstants(1,3,1,2,15,10,5);
v4OU = OU(4.0,0.71,60.)


Ncells = sysSize.Ne + sysSize.Ni + sysSize.N0




#inputTrains = readdlm("V4BothHemiInput_10sChunk_sim1.txt", ',')

#nsim = 0
#sigma0Array = collect(0:0.4:2.4)
#sigma0Array = collect(0.8:0.4:2.4)
#JRatioArray = collect(1:0.2:2.4)

#for sigma0i = 1:length(sigma0Array)
	#sigma0 = sigma0Array[sigma0i]
	#for JRatioi = 1:length(JRatioArray)
		#JR = JRatioArray[JRatioi]
		#for iRep = 1:3
			#nsim+=1
#for iWiring = 1:1
# 	iWiring = 1
# 	sigma0 = 0.71# Sigma of OU process. Controls degree of correlation in lateralized inputs.
# 	tau0 = 60.
# 	JR = 1.0
#
# mue = -1.1
# mui = -1.
# q = 0.75

include("simTwoPopHemiInputWeakRecurrentCoupling.jl")
include("simTwoPopHemiInputNoCoupling.jl")
include("simTwoPopHemiInputUnpack_StrongRec.jl")
include("simTwoPopHemiInput.jl")
include("runSimsFrozenGraph.jl")
include("runSimsFrozenGraphLinResp.jl")

times,ns,times0,ns0,weights,bias,connStrength = simTwoPopHemiInputUnpack_NoCoupleInit(simParams,sysSize,connProbs,taus,v4OU)

times,ns,times0,ns0,weights,voltageOverTime,bias,connStrength = simTwoPopHemiInputUnpack_WeakCoupleInit(simParams,sysSize,connProbs,taus,v4OU)

times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength = simTwoPopHemiInputUnpack_StrongRecSymmClusters(simParams,sysSize,connProbs,taus,v4OU)

times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength,v4OU = simTwoPopHemiInputUnpack_StrongRecAsymmClusters(simParams,sysSize,connProbs,taus,v4OU)

println("mean excitatory firing rate: ",mean(1000*ns[1:sysSize.Ne]/simParams.T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(sysSize.Ne+1):(Ncells-sysSize.N0)]/simParams.T)," Hz")

numSims = 200

@time runSimsFrozenGraphLinResp(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)

@time runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)








			times,ns,times0,ns0,weights,voltageOverTime = simTwoPopHemiInputWeakRecurrentCoupling(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR, q)

			tic()
			times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime,voltageOverTime,privNoiseOverTime = simTwoPopHemiInput_privateNoise(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR)
			toc()
			writecsv("spHemi_JR_1_6_JEI_1_1_sig01_4_tau0_60_init$iWiring.csv",times)
			writecsv("nsHemi_JR_1_6_JEI_1_1_sig01_4_tau0_60_init$iWiring.csv",ns)
			writecsv("weights_JR_1_6_JEI_1_1_sig01_4_tau0_60_init$iWiring.csv",weights)
			writecsv("sp0_JR_1_6_JEI_1_1_sig01_4_tau0_60_init$iWiring.csv",times0)
			writecsv("ns0_JR_1_6_JEI_1_1_sig01_4_tau0_60_init$iWiring.csv",ns0)
			writecsv("meanInput_JR_1_6_JEI_1_1_sig01_4_tau0_60_init$iWiring.csv",synInputPerNeuronOverTime)

			#writecsv("ns0_sim$nsim.csv",ns0)
			figure(figsize=(4,4))
			for ci = 1:Ne
				vals = times[ci,1:ns[ci]]
				y = ci*ones(length(vals))
				scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
			end
			xlim(T/4,T/2)
			ylim(0,Ne)
			ylabel("Neuron")
			xlabel("Time")
			tight_layout()
			savefig("CheckCheck2.png",dpi=150)
			#"PFCRaster_5_28_20_JR_1_4_JIE_1_1_sig02_0_tau0_60_freeze1$iWiring.png"

			figure(figsize=(4,4))
			for ci0 = 1:N0
				vals = times0[ci0,1:ns0[ci0]]
				y = ci0*ones(length(vals))
				scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
			end
			xlim(0,T/4)
			ylim(0,Ne)
			ylabel("Neuron")
			xlabel("Time")
			tight_layout()
			savefig("InputCheck.png",dpi=150)




#end
CSV.write("privNoiseOverTime_p1_p8K_sig_1_2.csv",DataFrame(privNoiseOverTime),writeheader=false)





#writecsv("spHemi_R2_$nsim.csv",times)
writecsv("spHemi_RJ2point0_10s_OUsig0_sim38.csv",times)
writecsv("nsHemi_RJ2point0_10s_OUsig0_sim38.csv",ns)
writecsv("weights_sim38_OUsig0.csv",weights)
writecsv("sp0_RJ2point0_10s_OUsig0_sim38.csv",times0)
writecsv("ns0_RJ2point0_10s_0Usig0_sim38.csv",ns0)

if doplot
	println("creating raster plot")
	figure(figsize=(4,4))
	for ci = 1:Ne
		vals = times[ci,1:ns[ci]]
		y = ci*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	end
	xlim(0,T/4)
	ylim(0,Ne)
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()
	savefig("uncoupledTest_mue_minus1_mui_minus08_7_20_2020.png",dpi=150)

	figure(figsize=(4,4))
	for ci0 = 1:N0
		vals = times0[ci0,1:ns0[ci0]]
		y = ci0*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	end
	xlim(0,T/4)
	ylim(0,Ne)
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()
	savefig("UncouopledTest_7_2020.png",dpi=150)




	println("creating weight matrix plot")
	figure()
	pcolormesh(weights)
	colorbar()
	tight_layout()
	savefig("NoCoupleLinResp.png",dpi=150)

	#savefig("weightStatic_R1_$nsim.png",dpi=150)
end

end


@save "bothHemi_R1_SimSet0Balance.jld2" K KI N0 N0pop Ncells Ne Nepop Ni Nipop T taue tauedecay tauerise taui tauidecay tauirise weights


check = @load "E:\\Doiron Lab\\Sim Data\\08_11_2020_244b2a50-db94-11ea-15b9-c9962d3d462d\\08_11_2020_freeze244b2a50-db94-11ea-15b9-c9962d3d462d.jld2"
