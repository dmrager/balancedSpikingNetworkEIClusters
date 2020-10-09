function simTwoPopHemiInputUnpack_StrongRecSymmClusters(simParams,sysSize,connProbs,taus,v4OU)

	connStrength = ConnStrength(true,1.25,0.0);

	bias = Bias(1.1,1.2,1,1.05,0.0);

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
	times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime = simTwoPopHemiInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, v4OU, JR, bias)
	return times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength


end

function simTwoPopHemiInputUnpack_StrongRecAsymmClusters(simParams,sysSize,connProbs,taus,v4OU)

	connStrength = ConnStrength(true,1.25,0.0);

	bias = Bias(0.99,0.98,0.984,0.984,0.0);

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
	times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime = simTwoPopHemiInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, v4OU, JR, bias)
	return times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength


end
