topLevel = pwd;

[subDirsNames] = GetSubDirsFirstLevelOnly(pwd);
cd(subDirsNames{2});

cat = 3;


run('analyzeJuliaModelSpikeTrain_LinResponse_CollapseBins_NoStates.m')

save('preLinResp_OctReformat.mat')
run('DanielleFILinResponse.m')

diagCovY0 = diag(covY);
clearvars -except diagCovY0 subDirsNames
save('diagCovY0_OctReformat.mat','diagCovY0')

cd('..')
cd(subDirsNames{1});

cat = 3;


run('analyzeJuliaModelSpikeTrain_LinResponse_CollapseBins_NoStates.m')

save('linRespInitialPrep_OctReformat.mat')
run('analyzeJuliaModelSpikeTrainFASimState.m');
run('gaussianMixtureClusterForStateDetection.m');
close all;
run('gaussianMixtureDivisionLinResponse.m');
save('preLinRespandGauss_OctReformat.mat')
run('LinResponseState_SLMFit.m');

cat = 3;

function [subDirsNames] = GetSubDirsFirstLevelOnly(parentDir)
% Get a list of all files and folders in this folder.
files    = dir(parentDir);
names    = {files.name};
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir] & ~strcmp(names, '.') & ~strcmp(names, '..');
% Extract only those that are directories.
subDirsNames = names(dirFlags);
end