function getPSP(pspArray,plotToggle)
  psp = []
  for neuronPSP in eachrow(pspArray)
    nonZeroIdx = findall(x->x>0,first.(neuronPSP))
    if first.(neuronPSP)[nonZeroIdx[1]] == 1
      psp = last.(neuronPSP[(nonZeroIdx[1]-3):(nonZeroIdx[1]+250)])
      plotToggle == true && plot(psp)
      break
    end
  end
  return psp
end


function getCurr(pspArray,neuronIdx,plotToggle)
  neuronInputCurr = last.(pspArray[neuronIdx,:])
  plotToggle == true && plot(neuronInputCurr)
  return neuronInputCurr
end

function getMeanCurr(pspArray,neuronIdxArray,plotToggle)
  neuronInputCurrs = last.(pspArray[neuronIdxArray,:])
  popMeanInputCurr = mean(neuronInputCurrs,dims=1)
  plotToggle == true && plot(vec(popMeanInputCurr))
  return vec(popMeanInputCurr)
end
