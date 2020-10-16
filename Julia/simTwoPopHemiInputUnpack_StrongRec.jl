function simTwoPopHemiInputUnpack_StrongRecSymmClusters(simParams,sysSize,connProbs,taus,v4OU)

	#connStrength = ConnStrength(true,2.3,0.0,1.0);

	connStrength = ConnStrength(true,1.5,0.0,1.0);


	#bias = Bias(1.0,1.1,1,1.05,0.05); #linn response segment R=1.25 **w/ J0 connections 3.2 instead of 3!!!

	bias = Bias(1.04,1.14,1,1.05,0.05);


	#bias = Bias(1.1,1.1,0.9,0.95,0.0); #works for FA


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
	times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime = simTwoPopHemiInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, v4OU, JR, bias,asymmCoeff)
	return times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength


end

function simTwoPopHemiInputUnpack_StrongRecAsymmClusters(simParams,sysSize,connProbs,taus,v4OU)

	#connStrength = ConnStrength(true,1.4,0.0,1.1);
	connStrength = ConnStrength(true,1.7,0.0,1.0);

	#v4OU = OU(4.0,2.0,60.);
	v4OU = OU(4.0,0.71,60.);

	#bias = Bias(1.2,1.3,0.984,0.984,0.0);

	bias = Bias(1.1,1.2,0.984,0.984,0.0);


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
	times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime = simTwoPopHemiInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, v4OU, JR, bias, asymmCoeff)
	return times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength, v4OU
end
