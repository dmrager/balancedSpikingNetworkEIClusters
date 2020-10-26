




downsampRate = 100;
dt = 0.1; %in ms
downsampTime = downsampRate * dt;

Ne = 2000;
Ni = 500;
Nc = 5000;

Ninput = 4000;

T = 12000;
time=0:1:(T-1);
Tw=150;
Tw_bin=50;
Tburn = 500;
binedges=Tburn+1:Tw_bin:T;


%figure;

searchStr = 'spHemi_.*.csv';
searchStrNs = 'nsHemi_.*.csv';
synInStr = 'meanInput_.*.csv';
%synInStr = 'meanInput_.*.h5';

inputStrSp = 'sp0_.*.csv';
inputStrNs = 'ns0_.*.csv';


currentDir = pwd;
fnames = dir(fullfile(currentDir,'*.csv'));
fnames = {fnames.name}.';

%fnamesH5 = dir(fullfile(currentDir,'*.h5'));
%fnamesH5 = {fnamesH5.name}.';

FIND = @(str)cellfun(@(c)~isempty(c),regexp(fnames,str,'once'));
%FIND_H5 = @(str)cellfun(@(c)~isempty(c),regexp(fnamesH5,str,'once'));


fnames_subset = fnames(FIND(searchStr));
fnames_subsetNs =  fnames(FIND(searchStrNs));
fnames_subsetSyn = fnames(FIND(synInStr));
%fnames_subsetSyn = fnamesH5(FIND_H5(synInStr));
fnames_inputNs = fnames(FIND(inputStrNs));
fnames_inputSp = fnames(FIND(inputStrSp));



re_all = [];
re_allInhib = [];
R_perPair_all = [];
R_perPairBoth_all = [];
R_perPair_allA = [];
R_perPair_allB = [];
countsAll = [];

reMeanHemi = [];
reMeanHemiInhib = [];


R_perPairBetweenHistAll = [];
R_perPairWithinHistAll = [];

currentEL_StateA_All = [];
currentIL_StateA_All = [];
currentER_StateA_All = [];
currentIR_StateA_All = [];

currentEL_StateB_All = [];
currentIL_StateB_All = [];
currentER_StateB_All = [];
currentIR_StateB_All = [];

FR_EL_StateA_All = [];
FR_IL_StateA_All = [];
FR_ER_StateA_All = [];
FR_IR_StateA_All = [];

FR_EL_StateB_All = [];
FR_IL_StateB_All = [];
FR_ER_StateB_All = [];
FR_IR_StateB_All = [];

inputAll = [];

lambdaIn = 4;
sigmaIn = 0.71;

covEstIn = ((sigmaIn * (Tw_bin*0.001)).^2)/2; 
varEstIn = (lambdaIn * (Tw_bin*0.001)) + covEstIn;
 
tic

simStruct(length(fnames_subset)).times = [];

numSims = length(fnames_subset);


parfor iSim = 1:numSims
    
    %FRFig = figure;
    
    
    fileStr = fnames_subset{iSim};
    fileStrNs = fnames_subsetNs{iSim};
    fileStrSyn = fnames_subsetSyn{iSim};
    %fileStrNs0 = fnames_inputNs{iSim};
    fileStrSp0 = fnames_inputSp{iSim};

    simStruct(iSim).times = csvread(fileStr);
    simStruct(iSim).ns = csvread(fileStrNs);
    simStruct(iSim).synInputMat = csvread(fileStrSyn);
    %simStruct(iSim).synInputMat = h5read(fileStrSyn,'/meanInput');

    simStruct(iSim).timesInput = csvread(fileStrSp0);
    

    
    
end
for iSim = 1:numSims
    %nsInput = csvread(fileStrNs0);
    %subplot(5,1,iSim);
    
    reBoth = [];
    reBothInhib = [];
    scBoth = [];
    
    times = simStruct(iSim).times;
    ns = simStruct(iSim).ns;
    synInputMat = simStruct(iSim).synInputMat;
    timesInput = simStruct(iSim).timesInput;
    
    timesAll = [];
    times0All = [];
    nNumAll = [];
    nNum0All = [];
    meanR = [];
    re1 = [];
    Y = [];
    counts = [];
    reMeanHemi = [];
    reMeanHemiInhib = [];

    
    for ci = 1:(Ne*2+Ni*2)
        timevals = times(ci,1:ns(ci));
        nNum = ci * ones(1,length(timevals));
        timesAll = [timesAll, timevals];
        nNumAll = [nNumAll,nNum];
%         if  ci < Ninput
%             timevalsInput = timesInput(ci,1:nsInput(ci));
%             nNumInput = ci * ones(1,length(timevalsInput));
%             times0All = [times0All,timevalsInput];
%             nNum0All = [nNum0All,nNumInput];
%         end
    end
    
    s1 = [timesAll;nNumAll];
    s1 = s1(:,timesAll<=T);
    
%    s1Input = [times0All;nNum0All];
    
    
    I1=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
    
    I1=I1(I1<=(Ne*2)+(Ni*2));
    
    countsAll = [countsAll, counts];
    
%     for mm=1:Ninput
%             reInput(mm,:)=hist(s1Input(1,mm-1/4<s1Input(2,:) & s1Input(2,:)<=mm+1/4),binedges);%*1e3;
%     end
%     
%     reInput = reInput * (1000/Tw_bin);

    
    
    for iCluster = 0:1
        
        %timesAll = [];
        %nNumAll = [];
        meanR = [];
        re1 = [];
        re1Inhib = [];
        
        %         for ci = ((iCluster * Ne) + 1):(Ne * (iCluster + 1))
        %             timevals = times(ci,1:ns(ci));
        %             nNum = ci * ones(1,length(timevals));
        %             timesAll = [timesAll, timevals];
        %             nNumAll = [nNumAll,nNum];
        %         end
        
        
        %timesAll = sxChunkBothHemi(1,:);
        
        %s1 = sxChunkBothHemi;
        
        %s1 = s1(:,timesAll<=10000);
        
        %I1_new=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
        
        %I1_new=I1_new(I1_new<=(Ne*2));
        
        if iCluster == 0 %grab spike time samples from a single cluster
            %Ic1=randsample(I1(I1<=Ne),Nc/2);
            Ic1 = I1(I1<=Ne);
            diffL = setdiff(1:Ne,Ic1);
            Ic1Inhib = I1(I1>(2*Ne)&I1<=(2*Ne+Ni));
            diffLInhib = setdiff(((2*Ne+1):(2*Ne+Ni)),Ic1Inhib);
            blackShade = 'k';
            redShade = [.91 .31 .22];
            
        else
            %Ic1=randsample(I1(I1>Ne),Nc/2);
            
            Ic1 = I1(I1>Ne&I1<=(Ne*2));
            diffR = setdiff((Ne+1):(Ne*2),Ic1);
            Ic1Inhib = I1(I1>(2*Ne+Ni)&I1<=(2*Ne+2*Ni));
            diffRInhib = setdiff(((2*Ne+Ni+1):(2*Ne+2*Ni)),Ic1Inhib);
            blackShade = [.63 .63 .64];
            redShade = [.99 .62 .62];
            
        end
        
        
        
        cat = 3;
        
        %for mm=1:Nc/2 %get spike counts per time window for single cluster
        parfor mm=1:length(Ic1)
            %             if ismember(mm,[diffLInhib diffRInhib])
            %                re1(mm,:) = zeros(1,size(re1,2));
            %             else
            re1(mm,:)=histc(s1(1,Ic1(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1(mm)+1/4),binedges);%*1e3;
            %             end
        end
        
        %re1_s=imfilter(re1(:,Tburn+1:end),ones(1,Tw)/Tw);
        re1_s = re1 * (1000/Tw_bin);
        
        parfor mm=1:length(Ic1Inhib)
            %             if ismember(mm,[diffLInhib diffRInhib])
            %                 re1Inhib(mm,:)=zeros(1,size(re1Inhib,2));
            %             else
            re1Inhib(mm,:)=histc(s1(1,Ic1Inhib(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1Inhib(mm)+1/4),binedges);%*1e3;
            %             end
        end
        
        %re1_sInhib=imfilter(re1Inhib(:,Tburn+1:end),ones(1,Tw)/Tw);
        re1_sInhib = re1Inhib * (1000/Tw_bin);
        
        
%         set(0,'CurrentFigure',FRFig);
%         subplot(4,1,3);
%         plot(binedges,mean(re1_s),'linewidth',3,'color',blackShade);
%         axis([plotTimeMin plotTimeMax 0 10]);
%         %title('Mean E Pop FR (Hz)')
%         ylabel('Mean FR E Pops (Hz)','fontsize',16);
%         xlabel('Time (ms)','fontsize',16);
%         hold all;
%         box off;
%         subplot(4,1,4);
%         plot(binedges,mean(re1_sInhib),'linewidth',3,'color',redShade);
%         hold all;
%         %title('Mean I Pop FR (Hz)')
%         axis([plotTimeMin plotTimeMax 0 24]);
%         ylabel('Mean FR I Pops (Hz)','fontsize',16);
%         box off;
        
        
        reMeanHemi = [reMeanHemi; mean(re1)];
        reMeanHemiInhib = [reMeanHemiInhib; mean(re1Inhib)];
        
       scBoth = [scBoth; re1; re1Inhib];
        

        
        reBoth = [reBoth; re1_s; re1_sInhib];
        %reBothInhib = [reBothInhib; re1_sInhib];
        
        re1_s_reshape = reshape(mean(re1_s'),size(mean(re1_s'),1)*size(mean(re1_s'),2),1);
        re1_s_reshape = re1_s_reshape';
        re_all = [re_all, re1_s_reshape];
        
        re1_s_reshapeInhib = reshape(mean(re1_sInhib'),size(mean(re1_sInhib'),1)*size(mean(re1_sInhib'),2),1);
        re1_s_reshapeInhib = re1_s_reshapeInhib';
        re_allInhib = [re_allInhib, re1_s_reshapeInhib];
        
        
        
        
        %         COV=cov(re1_s');
        %         Var=diag(COV);
        %         %rate1=mean(re1_s,2);
        %
        %         R = COV./sqrt(Var*Var');
        %         R_perPair = triu(R,1);
        %         R_perPair(R_perPair==0)=NaN;
        %
        %
        %         %R_perPair_Within = R_perPair(1:(N0/2),1:(N0/2))
        %
        %         R_perPairReshape = reshape(R_perPair,size(R_perPair,1).^2,1);
        %
        %         R_perPair_all = [R_perPair_all; R_perPairReshape];
        %
        %         if iCluster == 0
        %             R_perPair_allA = [R_perPair_allA; R_perPairReshape];
        %         else
        %             R_perPair_allB = [R_perPair_allB; R_perPairReshape];
        %         end
        %
        %         %R_perPair_allA = reshape(R_perPair_allA,size(R_perPair_allA,1).^2,1);
        %         %R_perPair_allB = reshape(R_perPair_allB,size(R_perPair_allB,1).^2,1);
        %
        %
        %         %meanR((iCluster+1)) = nanmean(R_perPair);
        %
    end
    %
    %     COVBoth=cov(reBoth');
    %     VarBoth=diag(COVBoth);
    %     %rate1=mean(re1_s,2);
    %
    %     RBoth = COVBoth./sqrt(VarBoth*VarBoth');
    %     R_perPairBoth = triu(RBoth,1);
    %     R_perPairBoth(R_perPairBoth==0)=NaN;
    %
    %     R_perPairWithinA = reshape(R_perPairBoth(1:(N0),1:(N0)),(N0).^2,1);
    %     R_perPairWithinB = reshape(R_perPairBoth((N0 + 1):(N0*2),(N0 + 1):(N0*2)),(N0).^2,1);
    %     R_perPairBetween = reshape(R_perPairBoth(1:(N0),(N0 + 1):(N0*2)),(N0).^2,1);
    %     R_perPairWithin = [R_perPairWithinA; R_perPairWithinB];
    %     R_perPairBetweenPad = padarray(R_perPairBetween,size(R_perPairBetween),NaN,'post');
    %     withinMinusBetween = nanmean(R_perPairWithin) - nanmean(R_perPairBetween);
    %     %R_withinMinusBetweenAll(R_Jidx,sigma0idx) = nanmean([R_withinMinusBetweenAll(R_Jidx,sigma0idx),withinMinusBetween]);
    %     %R_WithinAll(R_Jidx,sigma0idx) = nanmean([R_WithinAll(R_Jidx,sigma0idx),nanmean(R_perPairWithin)]);
    %
    %     R_perPairWithinHistAll = [R_perPairWithinHistAll; R_perPairWithin];
    %     R_perPairBetweenHistAll = [R_perPairBetweenHistAll; R_perPairBetweenPad];
    %
    %     %figure;
    %     %imagesc(R_perPairBoth);
    %    R_perPairBoth = reshape(R_perPairBoth,size(R_perPairBoth,1).^2,1);
    %     R_perPairBoth_all = [R_perPairBoth_all; R_perPairBoth];
    
    positiveSetEL = setdiff(1:2000,diffL);
    positiveSetIL = setdiff(4001:4500,diffLInhib);
    positiveSetER = setdiff(2001:4000,diffR);
    positiveSetIR = setdiff(4501:5000,diffRInhib);
    
    positiveSetAll = [positiveSetEL positiveSetIL positiveSetER positiveSetIR];
    
    [row, col] = find(reMeanHemi(1,:)>0.25 & reMeanHemi(2,:)<0.25 & reMeanHemi(2,:)>0.08);
    [row2, col2] = find(reMeanHemi(2,:)>0.25 & reMeanHemi(1,:)<0.25 & reMeanHemi(1,:)>0.08);
    
    [rowADS, colADS] = find(reMeanHemi(1,:)>-3 & reMeanHemi(2,:)>-3); %100 is amount by which current inputs are downsampled
    [rowBDS, colBDS] = find(reMeanHemi(2,:)>0.25 & reMeanHemi(1,:)<0.25 & reMeanHemi(1,:)>0.08);
    
%     set(0,'CurrentFigure',FRFig);
%     subplot(4,1,1);
%     plot(binedges(col),ones(size(col)),'.g');
%     hold all;
%     plot(binedges(col2),ones(size(col2)),'.y');
%     axis([plotTimeMin plotTimeMax 0.9 1.1])
%     set(gca,'visible','off')
%     subplot(4,1,2);
%     plot(timesAll,nNumAll,'.k');
%     ylabel('Neuron Idx','fontsize',16);
%     axis([plotTimeMin plotTimeMax 1 5000]);
    
%     figure;
%     subplot(2,1,1);
%     plot(times0All,nNum0All,'.k');
%     ylabel('Neuron Idx','fontsize',16);
%     title('Input Spikes','fontsize',20);
%     axis([plotTimeMin plotTimeMax 1 4000]);    
%     subplot(2,1,2);
%     plot(mean(reInput(1:2000,2:end)));
%     hold all;
%     plot(mean(reInput(2001:4000,2:end)));
%     xlim([plotTimeMin/Tw_bin plotTimeMax/Tw_bin]); %the alignment on this isn't correct 

    

    
    
    

    
    
    
    %  synInputMatPreppedEL_2 = synInputMat(positiveSetEL,(Tburn/downsampTime):end);
    %  synInputMatPreppedIL_2 = synInputMat(positiveSetIL,(Tburn/downsampTime):end);
    %  synInputMatPreppedER_2 = synInputMat(positiveSetER,(Tburn/downsampTime):end);
    %  synInputMatPreppedIR_2 = synInputMat(positiveSetIR,(Tburn/downsampTime):end);
    
    
    synInputMatPreppedEL = synInputMat(positiveSetEL,(Tburn/downsampTime)+1:Tw_bin/downsampTime:end);
    synInputMatPreppedIL = synInputMat(positiveSetIL,(Tburn/downsampTime)+1:Tw_bin/downsampTime:end);
    synInputMatPreppedER = synInputMat(positiveSetER,(Tburn/downsampTime)+1:Tw_bin/downsampTime:end);
    synInputMatPreppedIR = synInputMat(positiveSetIR,(Tburn/downsampTime)+1:Tw_bin/downsampTime:end);
    
    synInputMatPreppedEL_cum = synInputMat(positiveSetEL,:);
    synInputMatPreppedIL_cum = synInputMat(positiveSetIL,:);
    synInputMatPreppedER_cum = synInputMat(positiveSetER,:);
    synInputMatPreppedIR_cum = synInputMat(positiveSetIR,:);
    
    synInputMatPreppedEL_cum = downSampCumSum(synInputMatPreppedEL_cum(:,(Tburn/downsampTime)+1:end),Tw_bin/downsampTime);
    synInputMatPreppedIL_cum = downSampCumSum(synInputMatPreppedIL_cum(:,(Tburn/downsampTime)+1:end),Tw_bin/downsampTime);
    synInputMatPreppedER_cum = downSampCumSum(synInputMatPreppedER_cum(:,(Tburn/downsampTime)+1:end),Tw_bin/downsampTime);
    synInputMatPreppedIR_cum = downSampCumSum(synInputMatPreppedIR_cum(:,(Tburn/downsampTime)+1:end),Tw_bin/downsampTime);



    
    
    
    
    
    currentEL_StateA = zeros(Ne,length(colADS));
    currentIL_StateA = zeros(Ni,length(colADS));
    currentER_StateA = zeros(Ne,length(colADS));
    currentIR_StateA = zeros(Ni,length(colADS));
    
    FR_EL_StateA = currentEL_StateA;
    FR_IL_StateA = currentIL_StateA;
    FR_ER_StateA = currentER_StateA;
    FR_IR_StateA = currentIR_StateA;
    
    currentEL_StateB = zeros(Ne,length(colBDS));
    currentIL_StateB = zeros(Ni,length(colBDS));
    currentER_StateB = zeros(Ne,length(colBDS));
    currentIR_StateB = zeros(Ni,length(colBDS));
    
    FR_EL_StateB = currentEL_StateB;
    FR_IL_StateB = currentIL_StateB;
    FR_ER_StateB = currentER_StateB;
    FR_IR_StateB = currentIR_StateB;
    
    
    
    currentEL_StateA(positiveSetEL,:) = synInputMatPreppedEL_cum(:,colADS);
    currentIL_StateA((positiveSetIL-4000),:) = synInputMatPreppedIL_cum(:,colADS);
    currentER_StateA((positiveSetER-2000),:) = synInputMatPreppedER_cum(:,colADS);
    currentIR_StateA((positiveSetIR-4500),:) = synInputMatPreppedIR_cum(:,colADS);
    
    currentEL_StateB(positiveSetEL,:) = synInputMatPreppedEL_cum(:,colBDS);
    currentIL_StateB((positiveSetIL-4000),:) = synInputMatPreppedIL_cum(:,colBDS);
    currentER_StateB((positiveSetER-2000),:) = synInputMatPreppedER_cum(:,colBDS);
    currentIR_StateB((positiveSetIR-4500),:) = synInputMatPreppedIR_cum(:,colBDS);
    
    %  reBothDS = reBoth(:,1:10:end);
    reBothDS = reBoth;
    
    
    FR_EL_StateA(positiveSetEL,:) = reBothDS((positiveSetAll<=2000),colADS);
    FR_IL_StateA((positiveSetIL-4000),:) = reBothDS((positiveSetAll>4000 & positiveSetAll<=4500),colADS);
    FR_ER_StateA((positiveSetER-2000),:) = reBothDS((positiveSetAll>2000 & positiveSetAll<=4000),colADS);
    FR_IR_StateA((positiveSetIR-4500),:) = reBothDS((positiveSetAll>4500 & positiveSetAll<=5000),colADS);
    
    FR_EL_StateB(positiveSetEL,:) = reBothDS((positiveSetAll<=2000),colBDS);
    FR_IL_StateB((positiveSetIL-4000),:) = reBothDS((positiveSetAll>4000 & positiveSetAll<=4500),colBDS);
    FR_ER_StateB((positiveSetER-2000),:) = reBothDS((positiveSetAll>2000 & positiveSetAll<=4000),colBDS);
    FR_IR_StateB((positiveSetIR-4500),:) = reBothDS((positiveSetAll>4500 & positiveSetAll<=5000),colBDS);
    
    %  FR_EL_StateA = reBothDS(1:Ne,colADS);
    %  FR_IL_StateA = reBothDS((Ne*2+1):(Ne*2+Ni),colADS);
    %  FR_ER_StateA = reBothDS((Ne+1):(Ne*2),colADS);
    %  FR_IR_StateA = reBothDS((Ne*2+Ni+1):(Ne*2+Ni*2),colADS);
    %
    %  FR_EL_StateB = reBothDS(1:Ne,colBDS);
    %  FR_IL_StateB = reBothDS((Ne*2+1):(Ne*2+Ni),colBDS);
    %  FR_ER_StateB = reBothDS((Ne+1):(Ne*2),colBDS);
    %  FR_IR_StateB = reBothDS((Ne*2+Ni+1):(Ne*2+Ni*2),colBDS);
    
    
    currentEL_StateA_All = [currentEL_StateA_All currentEL_StateA];
    currentIL_StateA_All = [currentIL_StateA_All currentIL_StateA];
    currentER_StateA_All = [currentER_StateA_All currentER_StateA];
    currentIR_StateA_All = [currentIR_StateA_All currentIR_StateA];
    
    currentEL_StateB_All = [currentEL_StateB_All currentEL_StateB];
    currentIL_StateB_All = [currentIL_StateB_All currentIL_StateB];
    currentER_StateB_All = [currentER_StateB_All currentER_StateB];
    currentIR_StateB_All = [currentIR_StateB_All currentIR_StateB];
    
    FR_EL_StateA_All = [FR_EL_StateA_All FR_EL_StateA];
    FR_IL_StateA_All = [FR_IL_StateA_All FR_IL_StateA];
    FR_ER_StateA_All = [FR_ER_StateA_All FR_ER_StateA];
    FR_IR_StateA_All = [FR_IR_StateA_All FR_IR_StateA];
    
    FR_EL_StateB_All = [FR_EL_StateB_All FR_EL_StateB];
    FR_IL_StateB_All = [FR_IL_StateB_All FR_IL_StateB];
    FR_ER_StateB_All = [FR_ER_StateB_All FR_ER_StateB];
    FR_IR_StateB_All = [FR_IR_StateB_All FR_IR_StateB];
    
    
    
    StateBAll = [FR_EL_StateB_All; FR_ER_StateB_All; FR_IL_StateB_All; FR_IR_StateB_All];
    StateAAll = [FR_EL_StateA_All; FR_ER_StateA_All; FR_IL_StateA_All; FR_IR_StateA_All];
   
    
    %inputAll = [inputAll reInput(:,2:end)];   
    
end

StateAAll = StateAAll./(1000/Tw_bin);
StateBAll = StateBAll./(1000/Tw_bin);

covStateA = cov(StateAAll');
covStateB = cov(StateBAll');

% inputAll = inputAll./(1000/Tw_bin);
% covInput = cov(inputAll');
% 
% covInputTriU = triu(covInput,1);
% covInputTriU(covInputTriU==0)=NaN;
% inputHemL = covInputTriU(1:2000,1:2000);
% inputHemL = inputHemL(:);
% inputHemR = covInputTriU(2001:4000,2001:4000);
% inputHemR = inputHemR(:);
% inputBetHems = covInputTriU(1:2000,2001:4000);
% inputBetHems = inputBetHems(:);
% 
% figure;
% subplot(2,2,1);
% [n,x] = hist2Stair(diag(covInput),20,'k');
% yLims = get(gca,'YLim');
% hold on;
% plot([varEstIn varEstIn],[0 yLims(2)],'r', 'linewidth',3);
% hold on;
% plot([nanmean(diag(covInput)) nanmean(diag(covInput))], [0 yLims(2)], 'k--','linewidth',3);
% subplot(2,2,2);
% [n,x] = hist2Stair(inputBetHems,20,'k');
% yLims = get(gca,'YLim');
% hold on;
% plot([0 0],[0 yLims(2)],'r', 'linewidth',3);
% hold on;
% plot([nanmean(inputBetHems) nanmean(inputBetHems)], [0 yLims(2)], 'k--','linewidth',3);
% subplot(2,2,3);
% [n,x] = hist2Stair(inputHemL,20,'k');
% yLims = get(gca,'YLim');
% hold on;
% plot([covEstIn covEstIn],[0 yLims(2)],'r','linewidth',3);
% hold on;
% plot([nanmean(inputHemL) nanmean(inputHemL)], [0 yLims(2)], 'k--','linewidth',3);
% subplot(2,2,4);
% [n,x] = hist2Stair(inputHemR,x,'k');
% yLims = get(gca,'YLim');
% hold on;
% plot([covEstIn covEstIn],[0 yLims(2)],'r','linewidth',3);
% hold on;
% plot([nanmean(inputHemR) nanmean(inputHemR)], [0 yLims(2)], 'k--','linewidth',3);

toc



























