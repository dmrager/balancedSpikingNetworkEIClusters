
function simTwoPopHemiInputFreeze_synInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui,weights,v4OU, JR, bias, progressReport)

	@unpack_Bias bias;
	@unpack_OU v4OU;


	Ncells = Ne + Ni + N0
	NRec = Ne + Ni

	vre = 0. #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = .1 #simulation timestep (ms)
	refrac = 5 #refractory period (ms)

	maxrate = 100 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

	mu = zeros(NRec)
	thresh = zeros(NRec)
	tau = zeros(NRec)

	mu[1:Ne] = (muemax-muemin)*rand(Ne) .+ muemin
	mu[(Ne+1):(Ncells-N0)] = (muimax-muimin)*rand(Ni) .+ muimin

	thresh[1:Ne] .= threshe
	thresh[(1+Ne):(Ncells-N0)] .= threshi

	tau[1:Ne] .= taue
	tau[(1+Ne):(Ncells-N0)] .= taui

	#tau0 is unpacked now
	mulambda0 = mulambda0 / 1000.0;
	sig0 = sigma0 / 1000.0;


	Npop = round(Int,Ne/Nepop)
	Npopin = round(Int,Ni/Nipop)


	maxTimes = round(Int,maxrate*T/1000)
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



	v = rand(NRec) #membrane voltage

	lastSpike = -100*ones(NRec) #time of last spike

	Nsteps = round(Int,T/dt)
	dSampAmount = 100
	NSynapticDownsamp = round(Int,Nsteps/dSampAmount)

	synInputPerNeuronOverTime = zeros(NRec,NSynapticDownsamp)

	progressReport == true && println("starting simulation")

	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			progressReport == true && @printf("\r%d%%",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] .= 0
		forwardInputsI[:] .= 0

		OUHemi1 = randn();
		OUHemi2 = randn();
		N0Hemi1End = Int(N0/2);
		N0Hemi2Start = Int(N0Hemi1End + 1);

		for c0i = 1:N0Hemi1End
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


			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput) + (sigPriv * sqrt(dt/tau[ci]) * randn())

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
	progressReport == true && @printf("\r")

	times = times[:,1:maximum(ns)]
	times0 = times0[:,1:maximum(ns0)]

	return times,ns,times0,ns0,weights,synInputPerNeuronOverTime
end


function simTwoPopHemiInputUnpack_FreezeConnInit(simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)

	@unpack_SimParams simParams;
	@unpack_NCount sysSize;
	@unpack_ConnProbs connProbs;
	@unpack_ConnStrength connStrength;
	@unpack_TimeConstants taus;
	@unpack_OU v4OU;

	K = Int(Ne * pee)
	KI = Int(Ni * pii)

	if connStrength.strong == true
		ffScaling = sqrt(K)
		recScaling = sqrt(K)
	else
		ffScaling = sqrt(K)
		recScaling = K * q
	end

	Nepop = Int(Ne/2)
	Nipop = Int(Ni/2)
	N0pop = Int(N0/2)

	progressReport = false

	#Juno.@enter simTwoPopHemiInputWeakRecurrentCoupling(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR, q, recScaling, ffScaling, connProbs)
	times,ns,times0,ns0,weights,synInputPerNeuronOverTime = simTwoPopHemiInputFreeze_synInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui,weights,v4OU, JR, bias, progressReport)

	return times,ns,times0,ns0,weights,synInputPerNeuronOverTime


end
