Nc_FA = 100;

numNeuronsHemiL = sum(I1<2001);


countsAllMeanL=mean(countsAllFullPop(1:numNeuronsHemiL,:));
countsAllMeanR=mean(countsAllFullPop((numNeuronsHemiL+1):end,:));
countsAllStdR=var(countsAllFullPop((numNeuronsHemiL+1):end,:));
countsAllStdL=var(countsAllFullPop(1:numNeuronsHemiL,:));
countsAllSumStats = [countsAllMeanL;countsAllStdL;countsAllMeanR;countsAllStdR];

idxLowActivity = find(countsAllSumStats(1,:)<0.02&countsAllSumStats(3,:)<0.02);
countsAllSumStats(:,idxLowActivity)=[];
countsAllFullPop(:,idxLowActivity)=[];


rng('default') % set seed for reproducibility
options = statset('Display','final','MaxIter',1000); gm = fitgmdist(countsAllSumStats',2,'CovarianceType','diagonal','SharedCovariance',true,'Options',options);
idx = cluster(gm,countsAllSumStats');
%idx = cluster(gm,Z.mean');
cluster1 = (idx == 1);
cluster2 = (idx == 2);

X = countsAllSumStats';
%X = Z.mean';
P = posterior(gm,X);

figure;
scatter(X(cluster2,1),X(cluster2,3),20,P(cluster2,1),'x')
hold on
scatter(X(cluster1,1),X(cluster1,3),20,P(cluster1,1),'o')
xlabel('Mean Spike Count Hemi L');
ylabel('Mean Spike Count Hemi R')

figure;
scatter(X(cluster2,2),X(cluster2,4),20,P(cluster2,1),'x')
hold on
scatter(X(cluster1,2),X(cluster1,4),20,P(cluster1,1),'o')
xlabel('Var Across Neuron SCs Hemi L');
ylabel('Var Across Neuron SCs Hemi R')


cluster1Idx = find(idx==1);
cluster2Idx = find(idx==2);

cluster1Good = cluster1Idx(P(cluster1,1)>0.93);
cluster2Good = cluster2Idx(P(cluster2,2)>0.93);

countsAllMeanLFilt = countsAllMeanL;
countsAllMeanLFilt(:,idxLowActivity) = [];
countsAllMeanRFilt = countsAllMeanR;
countsAllMeanRFilt(:,idxLowActivity) = [];

%[topNeuronsDiffFR,topRandNeuronsDiffFR,topNeuronsFF,topRandNeuronsFF,topNeuronsChorus,topRandNeuronsChorus] = rankNeuronsStateDivision(countsAllFullPop,I1,cluster1Idx,cluster2Idx,countsAllMeanLFilt,countsAllMeanRFilt);


transitions = setdiff(1:length(X),[cluster1Good; cluster2Good]);

minSampleLengthTest = min(length(cluster1Good),length(cluster2Good));



uniformRandTimeSamps = zeros(10,length(cluster1Good));
fracNoTransition = (length(cluster1Good)+length(cluster2Good))./length(X);
fracTransition = 1 - fracNoTransition;
nTimeBinsTransSample = floor(fracTransition * length(cluster1Good));
nTimeBinsState = round((fracNoTransition * length(cluster1Good))/2);


for i = 1:size(uniformRandTimeSamps,1)

randTimeStateL = randsample(cluster1Good,nTimeBinsState);
randTimeStateR = randsample(cluster2Good,nTimeBinsState);
randTimeTrans = randsample(transitions,nTimeBinsTransSample);

if length([randTimeStateL; randTimeStateR; randTimeTrans']) > length(cluster1Good)
    randTimeStateL = randTimeStateL(2:end);
elseif length([randTimeStateL; randTimeStateR; randTimeTrans']) < length(cluster1Good) 
    randTimeStateL = [randTimeStateL; randsample(cluster1Good,1)];
end

uniformRandTimeSamps(i,:) = [randTimeStateL; randTimeStateR; randTimeTrans'];


end




numNeuronsHemiR = sum(I1>2000);

randHemiL1 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL2 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL3 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL4 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL5 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL6 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL7 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL8 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL9 = randsample(numNeuronsHemiL,Nc_FA/2);
randHemiL10 = randsample(numNeuronsHemiL,Nc_FA/2);

randHemiR1 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR2 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR3 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR4 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR5 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR6 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR7 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR8 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR9 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);
randHemiR10 = randsample((numNeuronsHemiL+1):length(I1),Nc_FA/2);

greyV = [0.7 0.7 0.7];

figure;
subplot(2,2,1)
plot(mean(countsAllFullPop(randHemiL1,cluster1Good(300:400))),'k','linewidth',2)
hold on
plot(mean(countsAllFullPop(randHemiR1,cluster1Good(300:400))),'color',greyV,'linewidth',2)
subplot(2,2,2)
plot(mean(countsAllFullPop(randHemiL1,cluster2Good(300:400))),'k','linewidth',2)
hold on
plot(mean(countsAllFullPop(randHemiR1,cluster2Good(500:600))),'color',greyV,'linewidth',2)
subplot(2,2,3)
plot(mean(countsAllFullPop(randHemiL2,cluster1Good(500:600))),'k','linewidth',2)
hold on
plot(mean(countsAllFullPop(randHemiR2,cluster1Good(500:600))),'color',greyV,'linewidth',2)
subplot(2,2,4)
plot(mean(countsAllFullPop(randHemiL2,cluster2Good(500:600))),'k','linewidth',2)
hold on
plot(mean(countsAllFullPop(randHemiR2,cluster2Good(500:600))),'color',greyV,'linewidth',2)





function [topNeuronsDiffFR,topRandNeuronsDiffFR,topNeuronsFF,topRandNeuronsFF,topNeuronsChorus_Final,topRandNeuronsChorus_Final] = rankNeuronsStateDivision(countsAllFullPop,I1,cluster1Idx,cluster2Idx,countsAllMeanL,countsAllMeanR)

meanFRState1 = mean(countsAllFullPop(:,cluster1Idx)');
meanFRState2 = mean(countsAllFullPop(:,cluster2Idx)');

stateFRDiff = abs(meanFRState1 - meanFRState2) ./ sqrt(meanFRState1 .* meanFRState2);
FAidx = 1:length(I1);

[topNeuronsDiffFR, topRandNeuronsDiffFR] = rankSort(stateFRDiff,I1,FAidx);

meanFRAllStates = mean(countsAllFullPop');

FFAllStates = var(countsAllFullPop') ./ meanFRAllStates;

[topNeuronsFF, topRandNeuronsFF] = rankSort(FFAllStates,I1,FAidx);

topCommon_diffFR_FF = intersect(topNeuronsDiffFR,topNeuronsFF);

topNeuronsChorus = NaN(100,2);
topRandNeuronsChorus = NaN(100,2);

for hemiI = 1:2
    
if hemiI == 1
    hemiPSTH = countsAllMeanL;
else
    hemiPSTH = countsAllMeanR;
end

crossCorrPSTH = NaN(size(meanFRAllStates));
crossCorrLag = NaN(size(meanFRAllStates));

for i = 1:length(crossCorrPSTH)
    y = sgolayfilt(countsAllFullPop(i,:),3,5);
    [r,lags] = xcorr(y,hemiPSTH,2,'normalized');        
    [rNoNorm,lagsNoNorm] = xcorr(y,hemiPSTH,2);
    [maxrNoNorm,maxidxNoNorm] = max(rNoNorm);
    [maxr,maxidx] = max(r);
    crossCorrPSTH(i) = maxr; 
    crossCorrLag(i) = lags(maxidx);
    crossCorrPSTHNoNorm(i) = maxrNoNorm; 
    crossCorrLagNoNorm(i) = lags(maxidxNoNorm);
end

[topNeuronsChorus(:,hemiI), topRandNeuronsChorus(:,hemiI)] = rankSort(crossCorrPSTH,I1,FAidx);
%topNeuronsChorusNoNorm(:,hemiI) = rankSort(crossCorrPSTHNoNorm,I1,FAidx);



end


topNeuronsChorus_Final = [topNeuronsChorus(1:50,1) ; topNeuronsChorus(51:end,2)];

topRandNeuronsChorus_Final = [topRandNeuronsChorus(1:50,1) ; topRandNeuronsChorus(51:end,2)];

%topNeuronsChorus_FinalNoNorm = [topNeuronsChorusNoNorm(1:50,1); topNeuronsChorusNoNorm(51:end,2)];

%topCommon_diffFR_Chorus = intersect(topNeuronsChorus_FinalNoNorm,topNeuronsChorus_Final);


topCommon_diffFR_Chorus = intersect(topNeuronsDiffFR,topNeuronsChorus_Final);

topCommon_FF_Chorus = intersect(topNeuronsFF,topNeuronsChorus_Final);


topNeuronsPlot(topNeuronsDiffFR,countsAllMeanL,countsAllMeanR,countsAllFullPop);

topNeuronsPlot(topRandNeuronsDiffFR,countsAllMeanL,countsAllMeanR,countsAllFullPop);


%topNeuronsPlot(topNeuronsFF,countsAllMeanL,countsAllMeanR,countsAllFullPop);

topNeuronsPlot(topNeuronsChorus_Final,countsAllMeanL,countsAllMeanR,countsAllFullPop);

topNeuronsPlot(topRandNeuronsChorus_Final,countsAllMeanL,countsAllMeanR,countsAllFullPop);


cat = 3;

















end

function [topNeurons,topRandom] = rankSort(stat,I1,FAidx)

topRandCount = 1000;

cat = 4;

fullStatMat = [stat' I1 FAidx'];
fullStatSort = sortrows(fullStatMat,1,'descend');
sortL = fullStatSort((fullStatSort(:,2)<2001),:);
sortR = fullStatSort((fullStatSort(:,2)> 2000),:);


%s = RandStream('mt19937ar','Seed','shuffle');

randTopQuarterL = randsample(topRandCount,50);
randTopQuarterR = randsample(topRandCount,50);


topNeurons = [sortL(1:50,3); sortR(1:50,3)];

topRandom = [sortL(randTopQuarterL,3); sortR(randTopQuarterR,3)];



end

function topNeuronsPlot(topNeurons,countsAllMeanL,countsAllMeanR,countsAllFullPop)

timeSamp = randsample((length(countsAllMeanR) - 201),1);

figure;
plot(countsAllMeanL(timeSamp:(timeSamp+100)));
hold on;
plot(countsAllMeanR(timeSamp:(timeSamp+100)));
plot(mean(countsAllFullPop(topNeurons(1:50),timeSamp:(timeSamp+100))));
plot(mean(countsAllFullPop(topNeurons(51:end),timeSamp:(timeSamp+100))));


end


