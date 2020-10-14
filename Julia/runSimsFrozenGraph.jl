function runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)

	iWiring = string(uuid1(MersenneTwister()))

	dateStr = Dates.format(today(),"mm_dd_yyyy")

	simSetDir = string(dateStr,"_",iWiring)
	stringMkDir = string("..\\",simSetDir)
	mkdir(stringMkDir)
	pathString = string("..\\",simSetDir,"\\")

	include("simTwoPopHemiInput_FreezeConnections_synInput.jl")
	@everywhere include("simTwoPopHemiInput_FreezeConnections_synInput.jl")


	#start = time()
	@sync @distributed for iSimRepeat = 1:numSims
		times,ns,times0,ns0,weights_loc,synInputPerNeuronOverTime = simTwoPopHemiInputUnpack_FreezeConnInit(simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)

		# spStr = string(pathString,"spHemi_",dateStr,"_freeze",iWiring,"($iSimRepeat).h5")
		# h5open(spStr,"w") do file
		# 	write(file,"spHemi",times)
		# end
		#
		# nsStr = string(pathString,"nsHemi_",dateStr,"_freeze",iWiring,"($iSimRepeat).h5")
		# h5open(nsStr,"w") do file
		# 	write(file,"nsHemi",ns)
		# end
		#
		# sp0Str = string(pathString,"sp0_",dateStr,"_freeze",iWiring,"($iSimRepeat).h5")
		# h5open(sp0Str,"w") do file
		# 	write(file,"sp0",times0)
		# end
		#
		# ns0Str = string(pathString,"ns0_",dateStr,"_freeze",iWiring,"($iSimRepeat).h5")
		# h5open(ns0Str,"w") do file
		# 	write(file,"ns0",ns0)
		# end

		meanIStr = string(pathString,"meanInput_",dateStr,"_freeze",iWiring,"($iSimRepeat).h5")
		h5open(meanIStr,"w") do file
			write(file,"meanInput",synInputPerNeuronOverTime)
		end

		CSV.write(string(pathString,"spHemi_",dateStr,"_freeze",iWiring,"($iSimRepeat).csv"),DataFrame(times),writeheader=false)
		CSV.write(string(pathString,"nsHemi_",dateStr,"_freeze",iWiring,"($iSimRepeat).csv"),DataFrame(ns'),writeheader=false)
		#writecsv("weights_JR_3_2_sig00_7_tau0_60_freeze$iWiring($iSimRepeat)_1_20_20.csv",weights)
		CSV.write(string(pathString,"sp0_",dateStr,"_freeze",iWiring,"($iSimRepeat).csv"),DataFrame(times0),writeheader=false)
		CSV.write(string(pathString,"ns0_",dateStr,"_freeze",iWiring,"($iSimRepeat).csv"),DataFrame(ns0'),writeheader=false)
		#CSV.write(string(pathString,"meanInput_",dateStr,"_freeze",iWiring,"($iSimRepeat).csv"),DataFrame(synInputPerNeuronOverTime),writeheader=false)

		#weightsAbs = broadcast(abs,weights);
		#maximum(maximum(weightsAbs,dims=1))
		#weights_1=weights;
		#CSV.write("gephiWeightsAbs_R_1_0.csv",DataFrame(weightsAbs),writeheader=false)
		#CSV.write("voltage_weakRecCouple_Actualsig00_71_tau0_60_7_27_2020_freeze$iWiring($iSimRepeat).csv",DataFrame(voltageOverTime),writeheader=false)

		figure(figsize=(4,4))
		for ci = 1:sysSize.Ne
			vals = times[ci,1:ns[ci]]
			y = ci*ones(length(vals))
			scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		end
		xlim(simParams.T/4,simParams.T/2)
		ylim(0,sysSize.Ne)
		ylabel("Neuron")
		xlabel("Time")
		tight_layout()
		savefig(string(pathString,dateStr,"_PFCRaster_freeze",iWiring,"($iSimRepeat).png"),dpi=150)
		PyPlot.close()
	end
	#elapsed=time()-start
	#println(elapsed)

	figure()
	pcolormesh(weights)
	colorbar()
	tight_layout()
	savefig(string(pathString,dateStr,"weights_",iWiring,".png"),dpi=150)

	weightsStr = string(pathString,"weights_",dateStr,"_freeze",iWiring,".h5")
	h5open(weightsStr,"w") do file
		write(file,"weights",weights)
	end

	paramSaveStr = string(pathString,dateStr,"_freeze",iWiring,".jld2")
	@save paramSaveStr simParams sysSize connProbs connStrength taus v4OU bias weights dateStr iWiring

end
