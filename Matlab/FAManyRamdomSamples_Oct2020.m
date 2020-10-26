numTrials = 25;

maxSampleLengthTest = max(length(cluster1Good),length(cluster2Good));
meanSampleLengthTest = mean(length(cluster1Good),length(cluster2Good));
FAPartitionResults(numTrials).cluster1Good = cluster1Good;

parfor iTrial = 1:numTrials
    
    s = RandStream('mt19937ar','Seed','shuffle');
    
    randHemiL = randsample(s,numNeuronsHemiL,Nc_FA/2);
    randHemiR = randsample(s,(numNeuronsHemiL+1):length(I1),Nc_FA/2);
    
    
    dim = faDMR(countsAllFullPop([randHemiL; randHemiR'],1:maxSampleLengthTest),1,'blah');
    dimHalf = faDMR(countsAllFullPop([randHemiL; randHemiR'],1:meanSampleLengthTest),1,'blah');
    dimC1 = faDMR(countsAllFullPop([randHemiL; randHemiR'],cluster1Good),1,'blah');
    dimC2 = faDMR(countsAllFullPop([randHemiL; randHemiR'],cluster2Good),1,'blah');
    
    FAPartitionResults(iTrial).dim = dim;
    FAPartitionResults(iTrial).dimMeanLength = dimHalf;
    FAPartitionResults(iTrial).dimC1 = dimC1;
    FAPartitionResults(iTrial).dimC2 = dimC2;
    FAPartitionResults(iTrial).cluster1Good = cluster1Good;
    FAPartitionResults(iTrial).cluster2Good = cluster2Good;
    
    close all;
    
    

    
end

save('FAPartition_225Trials.mat','FAPartitionResults','countsAllFullPop','I1');