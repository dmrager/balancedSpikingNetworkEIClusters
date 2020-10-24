
function simTwoPopHemiInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui, v4OU, JR, bias, asymmCoeff)

	@unpack_Bias bias;
	@unpack_OU v4OU;

	println("setting up parameters")

	Ncells = Ne + Ni + N0
	NRec = Ne + Ni
	Npop = round(Int,Ne/Nepop)
	Npopin = round(Int,Ni/Nipop)

	#connection probabilities
	pee = .2
	pei = .5
	pie = .5
	pii = .5

	sqrtK = sqrt(K)

	jie = 4. /(taui*sqrtK)
	jei = -16. *1.2/(taue*sqrtK)
	jii = -16. /(taui*sqrtK)

	#set up connection probabilities within and without blocks
	ratioejee = JR
	jeeout = 10. /(taue*sqrtK)
	jeein = ratioejee*10. /(taue*sqrtK)

	ratioijii = JR
	jii_out = -16. /(taui*sqrtK)
	jii_in= ratioijii*-16. /(taui*sqrtK)

	ratiojei = JR * asymmCoeff;
	jei_out = -16. *1.2/(taue*sqrtK)
	jei_in = ratiojei*-16. *1.2/(taue*sqrtK)

	ratiojie = JR #* 1.1
	jie_out = 4. /(taui*sqrtK)
	jie_in= ratiojie*4. /(taui*sqrtK)

	je0 = 10. *3.2 /(taue*sqrtK) #WAS 3!!! #Change to 3.2 for lin resp!!
	ji0 = 10. /(taui*sqrtK)



	ratiopee = 1
	peeout = K/(Nepop*(ratiopee-1) + Ne)
	peein = ratiopee*peeout

	pe0 = 0.2

	ratiopii = 1
	pii_out = KI/(Nipop*(ratiopii-1)+Ni)
	pii_in = ratiopii * pii_out

	ratiopei = 1
	pei_out = KI/(Nipop*(ratiopei-1)+Ni)
	pei_in = ratiopei * pei_out

	ratiopie = 1
	pie_out = KI/(Nipop*(ratiopie-1)+Ni)
	pie_in = ratiopie * pie_out

	pi0 = 0.5

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

	progressReport = true

	times,ns,times0,ns0,weights,synInputPerNeuronOverTime = simTwoPopHemiInputFreeze_synInput(T,Ne,Ni,N0,K,KI,Nepop,Nipop,N0pop,tauerise,tauedecay,tauirise,tauidecay,taue,taui,weights,v4OU, JR, bias, progressReport)

	return times,ns,times0,ns0,Ne,Ncells,T,weights,synInputPerNeuronOverTime
end
