
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


check = @load "E:\\Doiron Lab\\Sim Data\\Symm\\sig0_p71\\Lin Resp State Division\\10_14_2020_93e550fe-0e94-11eb-0cdc-570974db1d4f\\10_14_2020_freeze93e550fe-0e94-11eb-0cdc-570974db1d4f.jld2"


@load "..\\09_15_2020_freezed5a16b30-f70e-11ea-22f1-5b042a6a80b1.jld2"

@load "..\\09_14_2020_freezeaa08a862-f706-11ea-0267-f177c298400e.jld2"
