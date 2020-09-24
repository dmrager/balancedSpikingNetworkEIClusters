function runSimsFrozenGraphLinResp(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)

  runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)
  weights[1:5000,1:5000].=0
  bias = Bias(-1.13,-1.13,-0.8,-0.8,0.0);
  numSims = 20;
  runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)


end
