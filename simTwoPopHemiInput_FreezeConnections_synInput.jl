function simTwoPopHemiInputFreeze_synInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui,weights,sigma0, JR)
	println("setting up parameters")
	#Ne = 4000
	#Ni = 1000
	#N0 = 4000 #was 2500
	Ncells = Ne + Ni + N0
	NRec = Ne + Ni

	#T = 60000 #simulation time (ms)
	#inputTrainsH1 = readdlm("V4Hemi1Input_Test2.txt", ',')
	#inputTrains = readdlm(inputFile, ',')
	#inputTrains = hcat(inputTrainsH1,inputTrainsH2)


	#taue = 15 #membrane time constant for exc. neurons (ms)
	#taui = 10

	#connection probabilities
	#pee = .2
	#pei = .5
	#pie = .5
	#pii = .5
	#pe0 = 0.1
	#pi0 = 0.05

	#K = 800 #average number of E->E connections per neuron pee * NE
	#KI = 500 #pii * NI
	sqrtK = sqrt(K)

	#Nepop = 2000 #fix to be fxn of Ne
	#Nipop = 500
	#N0pop = 2000

	#constant bias to each neuron type
	muemin = 1.1
	muemax = 1.2 # was 1.2
	#muemin = -1.1
	#muemax = -1.1
	#muimin = -1.0 # was 1
	muimin = 1
	#muimax = -1.0 # was 1.05
	muimax = 1.05 #.05

	vre = 0. #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = .1 #simulation timestep (ms)
	refrac = 5 #refractory period (ms)

	#synaptic time constants (ms)
	#tauerise = 1
	#tauedecay = 3
	#tauirise = 1
	#tauidecay = 2

	maxrate = 100 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

	mu = zeros(NRec)
	thresh = zeros(NRec)
	tau = zeros(NRec)

	mu[1:Ne] = (muemax-muemin)*rand(Ne) + muemin
	mu[(Ne+1):(Ncells-N0)] = (muimax-muimin)*rand(Ni) + muimin

	thresh[1:Ne] = threshe
	thresh[(1+Ne):(Ncells-N0)] = threshi

	tau[1:Ne] = taue
	tau[(1+Ne):(Ncells-N0)] = taui

	tau0 = 60.0
	mulambda0 = 4.0 / 1000.0;
	sig0 = sigma0 / 1000.0;


	Npop = round(Int,Ne/Nepop)
	Npopin = round(Int,Ni/Nipop)

	#weights = zeros(Ncells,Ncells)



	maxTimes = round(Int,maxrate*3*T/1000)
	times = zeros(NRec,maxTimes)
	times0 = zeros(N0,maxTimes)
	ns = zeros(Int,NRec)
	ns0 = zeros(Int,N0)

	forwardInputsE = zeros(NRec) #summed weight of incoming E spikes
	forwardInputsI = zeros(NRec)
	forwardInputsEPrev = zeros(NRec) #as above, for previous timestep
	forwardInputsIPrev = zeros(NRec)

	xerise = zeros(NRec) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(NRec)
	xirise = zeros(NRec)
	xidecay = zeros(NRec)
	lambda0 = zeros(N0)

	Nstim_1 = 2000
	Nstim_2 = 4500
	stimstr = 0/taue
	stimstart = 500
	stimend = T

	v = rand(NRec) #membrane voltage

	lastSpike = -100*ones(NRec) #time of last spike

	Nsteps = round(Int,T/dt)
	dSampAmount = 100
	NSynapticDownsamp = round(Int,Nsteps/dSampAmount)

	synInputPerNeuronOverTime = zeros(NRec,NSynapticDownsamp)

	println("starting simulation")

	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] = 0
		forwardInputsI[:] = 0
		#spikingInputs = inputTrains[2,find(inputTrains[1,:] .== t)]

		OUHemi1 = randn();
		OUHemi2 = randn();
		N0Hemi1End = Int(N0/2);
		N0Hemi2Start = Int(N0Hemi1End + 1);

		for c0i = 1:N0Hemi1End
			cat = 3;
			lambda0[c0i] += dt * (1/tau0) * (mulambda0 - lambda0[c0i]) + sig0 * sqrt(dt/tau0) * OUHemi1;
			if rand() < lambda0[c0i]*dt
				#build structure for input spike trains here
				ns0[c0i] = ns0[c0i]+1
				times0[c0i,ns0[c0i]] = t
				for cinputTarget = 1:NRec
					forwardInputsE[cinputTarget] += weights[cinputTarget,(Ne+Ni+Int(c0i))]
				end
			end
		end
		for c0i = N0Hemi2Start:N0
			lambda0[c0i] += dt * (1/tau0) * (mulambda0 - lambda0[c0i]) + sig0 * sqrt(dt/tau0) * OUHemi2;
			if rand() < lambda0[c0i]*dt
				#build structure for input spike trains here
				ns0[c0i] = ns0[c0i]+1
				times0[c0i,ns0[c0i]] = t
				for cinputTarget = 1:NRec
					forwardInputsE[cinputTarget] += weights[cinputTarget,(Ne+Ni+Int(c0i))]
				end
			end
		end

		###if ~isempty(spikingInputs)
			###for c0i in spikingInputs
				###for cinputTarget = 1:NRec
					###forwardInputsE[cinputTarget] += weights[cinputTarget,(Ne+Ni+Int(c0i))]
				###end
			###end
		###end



			for ci = 1:NRec
				#forwardInputsE[ci] += weights[ci,(Ne+Ni+c0i)]
				xerise[ci] += -dt*xerise[ci]/tauerise + forwardInputsEPrev[ci]
				xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]
				xirise[ci] += -dt*xirise[ci]/tauirise + forwardInputsIPrev[ci]
				xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]

				synInput = (xedecay[ci] - xerise[ci])/(tauedecay - tauerise) + (xidecay[ci] - xirise[ci])/(tauidecay - tauirise)

				if mod(ti,dSampAmount) == 0
					synInputPerNeuronOverTime[ci,Int(ti/dSampAmount)] = synInput*tau[ci]
				end

				if (ci < Nstim_1 || ci > Nstim_2) && (t > stimstart) && (t < stimend)
					synInput += stimstr;
				end

				if t > (lastSpike[ci] + refrac)  #not in refractory period
					v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput)

					if v[ci] > thresh[ci]  #spike occurred
						v[ci] = vre
						lastSpike[ci] = t
						ns[ci] = ns[ci]+1
						if ns[ci] <= maxTimes
							times[ci,ns[ci]] = t
						end

						for j = 1:NRec
							if weights[j,ci] > 0  #E synapse
								forwardInputsE[j] += weights[j,ci]
							elseif weights[j,ci] < 0  #I synapse
								forwardInputsI[j] += weights[j,ci]
							end
						end #end loop over synaptic projections
					end #end if(spike occurred)
				end #end if(not refractory)
			end #end loop over neurons

		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
	end #end loop over time
	@printf("\r")

	times = times[:,1:maximum(ns)]
	times0 = times0[:,1:maximum(ns0)]

	cat = 3

	return times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime
end
