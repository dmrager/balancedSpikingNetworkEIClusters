clear all

run('analyzeJuliaModelSpikeTrain.m')
clear all

run('analyzeJuliaModelSpikeTrainFASimState.m')
save('preFApreGaussianMixture.mat')
run('gaussianMixtureClusterForStateDetection.m')

run('FAManyRamdomSamples_Oct2020.m')



