#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information


function simTwoPopHemiInputWeakRecurrentCoupling(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR, q, recScaling, ffScaling, connProbs, bias)

	@unpack_Bias bias;

	println("setting up parameters")

	Ncells = Ne + Ni + N0
	NRec = Ne + Ni

	sqrtK = sqrt(K)

	#set up connection probabilities within and without blocks
	ratioejee = JR
	jeeout = 10. /(taue*recScaling)
	jeein = ratioejee*jeeout

	#je0 = jeeout * (14./8.) #ratios taken from Chung Chung's paper


	ratioijii = JR
	jii_out = -16. /(taui*recScaling)
	jii_in= ratioijii*jii_out

	ratiojei = JR
	jei_out = -16. *1.2 /(taue*recScaling)
	jei_in = ratiojei*jei_out

	ratiojie = JR
	jie_out = 4. /(taui*recScaling)
	jie_in= ratiojie*jie_out

	#ji0 = jie_out * (12./4.)
	je0 = 10. *3 /(taue*ffScaling)
	ji0 = 10. /(taui*ffScaling)



	ratiopee = 1
	peeout = connProbs.pee
	peein = ratiopee*peeout

	#pe0 = peeout * 1
	pe0 = connProbs.pe0

	ratiopii = 1
	pii_out = connProbs.pii
	pii_in = ratiopii * pii_out

	ratiopei = 1
	pei_out = connProbs.pei
	pei_in = ratiopei * pei_out

	ratiopie = 1
	pie_out = connProbs.pie
	pie_in = ratiopie * pie_out

	pi0 = connProbs.pi0

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

	mu[1:Ne] = (muemax-muemin)*rand(Ne) .+ muemin
	mu[(Ne+1):(Ncells-N0)] = (muimax-muimin)*rand(Ni) .+ muimin

	thresh[1:Ne] .= threshe
	thresh[(1+Ne):(Ncells-N0)] .= threshi

	tau[1:Ne] .= taue
	tau[(1+Ne):(Ncells-N0)] .= taui

	tau0 = 60.0
	mulambda0 = 4.0 / 1000.0;
	sig0 = sigma0 / 1000.0;

	Npop = round(Int,Ne/Nepop)
	Npopin = round(Int,Ni/Nipop)

	weights = zeros(Ncells,Ncells)

	#random connections
	weights[1:Ne,1:Ne] = jeeout*(rand(Ne,Ne) .< peeout)
	weights[1:Ne,(1+Ne):(Ncells-N0)] = jei_out*(rand(Ne,Ni) .< pei_out)
	weights[(1+Ne):(Ncells-N0),1:Ne] = jie_out*(rand(Ni,Ne) .< pie_out)
	weights[(1+Ne):(Ncells-N0),(1+Ne):(Ncells-N0)] = jii_out*(rand(Ni,Ni) .< pii_out)

	weights[1:Nepop,(Ncells-N0+1):(Ncells-N0pop)] = je0*(rand(Nepop,N0pop) .< pe0)
	weights[(1+Ne):(Ncells-N0-Nipop),(Ncells-N0+1):(Ncells-N0pop)] = ji0*(rand(Nipop,N0pop) .< pi0)

	weights[(1+Nepop):Ne,(Ncells-N0pop+1):Ncells] = je0*(rand(Nepop,N0pop) .< pe0)
	weights[(1+Ne+Nipop):(Ne+Ni),(Ncells-N0pop+1):Ncells] = ji0*(rand(Nipop,N0pop) .< pi0)

	#connections within cluster
	for pjk = 1:Npop
		ipopstart = 1 + Nepop*(pjk-1)
		ipopend = pjk*Nepop
		eClustMemb = rand(Nepop,Nepop) .< peein
		weights[ipopstart:ipopend,ipopstart:ipopend] = jeein*(eClustMemb)
	end

	for pjk_in = 1:Npopin
		ipopstart_in = (1+Ne) + Nipop*(pjk_in-1)
		ipopend_in = Ne + (pjk_in*Nipop)
		iClustMemb = rand(Nipop,Nipop) .< pii_in;
		weights[ipopstart_in:ipopend_in,ipopstart_in:ipopend_in] = jii_in*iClustMemb
	end

	for pjk_mixed = 1:Npopin
		ipopstart_e = 1 + Nepop*(pjk_mixed-1)
		ipopend_e = pjk_mixed*Nepop
		ipopstart_i = (1+Ne) + Nipop*(pjk_mixed-1)
		ipopend_i = Ne + (pjk_mixed*Nipop)
		weights[ipopstart_e:ipopend_e,ipopstart_i:ipopend_i] = jei_in*(rand(Nepop,Nipop) .< pei_in);
		weights[ipopstart_i:ipopend_i,ipopstart_e:ipopend_e] = jie_in*(rand(Nipop,Nepop) .< pie_in);


	end

	for ci = 1:Ncells
		weights[ci,ci] = 0
	end

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

	v = rand(NRec) #membrane voltage

	lastSpike = -100*ones(NRec) #time of last spike

	Nsteps = round(Int,T/dt)
	dSampAmount = 100
	NSynapticDownsamp = round(Int,Nsteps/dSampAmount)

	synInputPerNeuronOverTime = zeros(NRec,NSynapticDownsamp)
	voltageOverTime = zeros(5,NSynapticDownsamp)
	privNoiseOverTime = zeros(5,NSynapticDownsamp)

	println("starting simulation")


	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] .= 0
		forwardInputsI[:] .= 0
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

				if t > (lastSpike[ci] + refrac)  #not in refractory period
					v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput) + (sigPriv * sqrt(dt/tau[ci]) * randn());

					if mod(ti,dSampAmount) == 0
						if ci < 6
							voltageOverTime[ci,Int(ti/dSampAmount)] = v[ci]
							privNoiseOverTime[ci,Int(ti/dSampAmount)] = (sigPriv * sqrt(dt/tau[ci]) * randn())
						end
					end

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

	return times,ns,times0,ns0,weights,voltageOverTime
end

function simTwoPopHemiInputUnpack_WeakCoupleInit(simParams,sysSize,connProbs,taus,v4OU)

	connStrength = ConnStrength(false,1,0.75,0.0);

	bias = Bias(-1.13,-1.13,-0.8,-0.8,0.1);

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

	#Juno.@enter simTwoPopHemiInputWeakRecurrentCoupling(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR, q, recScaling, ffScaling, connProbs)
	times,ns,times0,ns0,weights,voltageOverTime = simTwoPopHemiInputWeakRecurrentCoupling(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, sigma0, JR, q, recScaling, ffScaling, connProbs, bias)
	return times,ns,times0,ns0,weights,voltageOverTime,bias,connStrength


end
