
#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

#uncomment the line below and set doplot=true to plot a raster
using PyPlot
#using Plots

using JLD2

using HDF5

using Printf
using DataFrames
using CSV

#inputTrain = readdlm("V4Hemi1Input.txt", ',')

#inputTrain = readdlm("Users/dmrhome/Documents/ClusterSimsNew/BothHemiInput_600s_10sChunks_Sim0/V4BothHemiInput_10sChunk_sim1.txt", ',')


doplot = true

include("simTwoPopHemiInput.jl")
include("simTwoPopHemiInput_FreezeConnections_synInput.jl")

#for nsim = 2:2

T = 12000#simulation time (ms)

Ne = 4000
#Ne = 800
Ni = 1000
#Ni = 200

#N0 = 800

N0 = 4000

K = 800 #average number of E->E connections per neuron pee * NE
#K = 160
KI = 500
#KI = 40

Nepop = 2000 #fix to be fxn of Ne
#Nepop = 400
Nipop = 500
#Nipop = 100
N0pop = 2000

#N0pop = 400

tauerise = 1
tauedecay = 3
tauirise = 1
tauidecay = 2

taue = 15 #membrane time constant for exc. neurons (ms)
taui = 10

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
	iWiring = 1
	sigma0 = 0.71 # Sigma of OU process. Controls degree of correlation in lateralized inputs.
		JR = 1.1

			tic()
			times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime = simTwoPopHemiInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR)
			toc()
			writecsv("spHemi_JR_1_4_sig02_0_tau0_24_init$iWiring.csv",times)
			writecsv("nsHemi_JR_1_4_sig02_0_tau0_24_init$iWiring.csv",ns)
			writecsv("weights_JR_1_4_sig02_0_tau0_24_init$iWiring.csv",weights)
			writecsv("sp0_JR_1_4_sig02_0_tau0_24_init$iWiring.csv",times0)
			writecsv("ns0_JR_1_4_sig02_0_tau0_24_init$iWiring.csv",ns0)
			writecsv("meanInput_JR_3_4_sig02_0_tau0_24_init$iWiring.csv",synInputPerNeuronOverTime)

			#writecsv("ns0_sim$nsim.csv",ns0)
			@save "1_25_20_JR_2_3_sig00_71_tau0_60_init2.jld2" K KI N0 N0pop Ncells Ne Nepop Ni Nipop T taue tauedecay tauerise taui tauidecay tauirise weights sigma0 JR
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
			savefig("PFCRaster_11_27_19_JR_1_4_sig02_0_tau0_60_init5_$iWiring.png",dpi=150)

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
			savefig("PFCRaster_7_22_19_JR_3_2_sig00_7_NoCouple_tau0_60_init$iWiring.png",dpi=150)




#end

iSimRepeat = 1

iWiring = 1

start = time()
for iSimRepeat = 101:300


	times,ns,times0,ns0,Ne_loc,Ncells,T_loc,weights_loc,synInputPerNeuronOverTime = simTwoPopHemiInputFreeze_synInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui,weights,sigma0,JR)
	CSV.write("spHemi_JR_1_0_sig01_0_tau0_60_freeze$iWiring($iSimRepeat)_1_29_20.csv",DataFrame(times),writeheader=false)
	CSV.write("nsHemi_JR_1_0_sig01_0_tau0_60_freeze$iWiring($iSimRepeat)_1_29_20.csv",DataFrame(ns'),writeheader=false)
	#writecsv("weights_JR_3_2_sig00_7_tau0_60_freeze$iWiring($iSimRepeat)_1_20_20.csv",weights)
	CSV.write("sp0_JR_1_0_sig01_0_tau0_60_freeze$iWiring($iSimRepeat)_1_29_20.csv",DataFrame(times0),writeheader=false)
	CSV.write("ns0_JR_1_0_sig01_0_tau0_60_freeze$iWiring($iSimRepeat)_1_29_20.csv",DataFrame(ns0'),writeheader=false)
	CSV.write("synInput_JR_1_0_sig01_0_tau0_60_freeze$iWiring($iSimRepeat)_1_29_20.csv",DataFrame(synInputPerNeuronOverTime),writeheader=false)

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
	savefig("PFCRaster_1_29_20_JR_1_1_sig00_71_tau0_60_freeze$iWiring($iSimRepeat).png",dpi=150)
	PyPlot.close()




end
elapsed=time()-start
println(elapsed)




println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):(Ncells-N0)]/T)," Hz")

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
	savefig("RJ2point0_OUInput_lambda0_4_sigma0_0_sim38.png",dpi=150)

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
	savefig("OUInput_lambda0_4_sig_1point0.png",dpi=150)




	println("creating weight matrix plot")
	figure()
	pcolormesh(weightsNoCoupleR)
	colorbar()
	tight_layout()
	savefig("NoCoupling_1_20_20.png",dpi=150)

	#savefig("weightStatic_R1_$nsim.png",dpi=150)
end

end


@save "bothHemi_R1_SimSet0Balance.jld2" K KI N0 N0pop Ncells Ne Nepop Ni Nipop T taue tauedecay tauerise taui tauidecay tauirise weights

h5open("synInputH5Test.h5","w") do file
	write(file,"synInputMat",synInputPerNeuronOverTime)
end

@load "1_24_20_JR_1_0_sig00_71_tau0_60_init1.jld2"
