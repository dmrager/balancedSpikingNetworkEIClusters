tic

TwArray = [50,250];

for iTw = 1:length(TwArray) 
    
Tw = TwArray(iTw);   
    
Ne = 2000; %Cells in 1 E pop
Nc = 500; %Sample size for corr calcs, between both pops
Nc_FA = 50; %Sample size for FA
N0 = Nc/2; %Sample size for corr calcs within single pop

        T = 12000; %sim time in ms
        Tburn = 500; %0;
        time=0:1:(T-1);
        binedges=Tburn+1:Tw:T;

        %Tw=150; 
        %time=0:Tw:T;
        %Tw_FAprep=150;
        %Tburn = 500;
        
searchStr = 'spHemi_.*.csv';
searchStrNs = 'nsHemi_.*.csv';

%searchStrNs = 'nsHemi_\w*.csv';
%searchStr = 'sp0_R2_5_10s_OU_sim10.csv';
%searchStrNs = 'ns0_R2_5_10s_0U_sim10.csv';
%searchStr = 'spHemi_Balance_10s_OU_sim11.csv';
%searchStrNs = 'nsHemi_Balance_10s_OU_sim11.csv';


currentDir = pwd;
fnames = dir(fullfile(currentDir,'*.csv'));
fnames = {fnames.name}.';
FIND = @(str)cellfun(@(c)~isempty(c),regexp(fnames,str,'once'));
fnames_subset = fnames(FIND(searchStr));
fnames_subsetNs =  fnames(FIND(searchStrNs));

re_all = [];
R_perPair_all = [];
R_perPairBoth_all = [];
R_perPair_allA = [];
R_perPair_allB = [];
countsAll = [];

reMeanHemi = [];
VarBothAll = [];


simRefArray = zeros(3,length(fnames_subset));
R_JSweep = 1:0.2:2.4;
sigma0Sweep = 0:0.4:2.4;

R_withinMinusBetweenAll = NaN(length(R_JSweep),length(sigma0Sweep));
R_WithinAll = NaN(length(R_JSweep),length(sigma0Sweep));

R_perPairBetweenHistAll = [];
R_perPairWithinHistAll = [];

for iSim = 1:20%(round(length(fnames_subset)))  
         reBoth = [];
    
        fileStr = fnames_subset{iSim};
        fileStrNs = fnames_subsetNs{iSim};
        %simNumRegex = regexp(fileStr,'(^spHemi_)(\d*)(.csv)','tokens');
        %simNum = str2double(simNumRegex{1}(2));
        %sigma0idx = ceil(simNum/24);
        %R_Jidx = ceil(mod(ceil(simNum/3),8));
        %if R_Jidx == 0
            %R_Jidx = 8;
        %end
        %simRefArray(:,iSim) = [simNum; sigma0Sweep(sigma0idx); R_JSweep(R_Jidx)];
        times = csvread(fileStr);
        ns = csvread(fileStrNs);
        %subplot(5,1,iSim);
    
        timesAll = [];
        nNumAll = [];
        meanR = [];
        re1 = [];
        Y = [];
        counts = [];
    
        for ci = 1:(Ne*2)
            timevals = times(ci,1:ns(ci));
            nNum = ci * ones(1,length(timevals));
            timesAll = [timesAll, timevals];
            nNumAll = [nNumAll,nNum];
        end
    
    %nNumAll = sxChunkBothHemi(2,:);
    %timesAll = sxChunkBothHemi(1,:);
    
    %T = 10000;
    %TW_bin = 150;
    %Tburn = 0; %should probably add burn in
    
    s1 = [timesAll;nNumAll];
    s1 = s1(:,timesAll<=T);
    
    
    I1=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
    
    I1=I1(I1<=(Ne*2));
    
    %Ic1=randsample(I1(I1<=1000),Nc_FA/2); %why is this 1000? should be Ne? 
    %Ic2=randsample(I1(I1>1000),Nc_FA/2);
    %Icall = [Ic1;Ic2];
    
    %for mm=1:Nc_FA
        %counts(mm,:)=histc(s1(1,Icall(mm)-1/4<s1(2,:) & s1(2,:)<=Icall(mm)+1/4),binedges);
    %end
    
    countsAll = [countsAll, counts];
    %subplot(5,1,iSim);
    
    
    for iCluster = 0:1
        
        %timesAll = [];
        %nNumAll = [];
        meanR = [];
        re1 = [];
        
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
            Ic1=randsample(I1(I1<=Ne),Nc/2);
            %Ic1 = I1(I1<=Ne);
            diffL = setdiff(1:Ne,Ic1);
            
        else
            Ic1=randsample(I1(I1>Ne),Nc/2);
            
            %Ic1 = I1(I1>Ne);
            diffR = setdiff((Ne+1):(Ne*2),Ic1);
        end
        
        
        
        cat = 3;
        
        for mm=1:Nc/2 %get spike counts per time window for single cluster
        %for mm=1:size(Ic1)
            re1(mm,:)=histc(s1(1,Ic1(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1(mm)+1/4),binedges);%*1e3;
        end
        
        re1_s = re1;
        %re1_s=imfilter(re1(:,Tburn+1:end),ones(1,Tw)/Tw);
        %re1_s=re1_s(:,Tw/2-1:end-Tw/2); %smooth these to get FR for single cluster
        %plot(mean(re1_s));
        xlabel('Time (ms)');
        ylabel('Mean FR (Hz)');
        hold all;
        
        reMeanHemi = [reMeanHemi; mean(re1_s)];
        
        reBoth = [reBoth; re1_s];
        
        re1_s_reshape = reshape(mean(re1_s'),size(mean(re1_s'),1)*size(mean(re1_s'),2),1);
        re1_s_reshape = re1_s_reshape';
        re_all = [re_all, re1_s_reshape];
        
        
        
        
         COV=cov(re1_s');
         Var=diag(COV);
%         %rate1=mean(re1_s,2);
%         
        R = COV./sqrt(Var*Var');
        R_perPair = triu(R,1);
        
        R_perPair(R_perPair==0)=NaN;
        
        
        %R_perPair_Within = R_perPair(1:(N0/2),1:(N0/2))
        
        R_perPairReshape = reshape(R_perPair,size(R_perPair,1).^2,1);

        R_perPair_all = [R_perPair_all; R_perPairReshape];

        if iCluster == 0
            R_perPair_allA = [R_perPair_allA; R_perPairReshape];
        else
            R_perPair_allB = [R_perPair_allB; R_perPairReshape];
        end
        
%        R_perPair_allA = reshape(R_perPair_allA,size(R_perPair_allA,1).^2,1);
%       R_perPair_allB = reshape(R_perPair_allB,size(R_perPair_allB,1).^2,1);

        
        meanR((iCluster+1)) = nanmean(R_perPairReshape);
%         
    end
    
    COVBoth=cov(reBoth');
    VarBoth=diag(COVBoth);
    VarBothAll = [VarBothAll;VarBoth];
    %rate1=mean(re1_s,2);
    
    RBoth = COVBoth./sqrt(VarBoth*VarBoth');
    R_perPairBoth = triu(RBoth,1);
    R_perPairBoth(R_perPairBoth==0)=NaN;
    
    R_perPairWithinA = reshape(R_perPairBoth(1:(N0),1:(N0)),(N0).^2,1);
    R_perPairWithinB = reshape(R_perPairBoth((N0 + 1):(N0*2),(N0 + 1):(N0*2)),(N0).^2,1);
    R_perPairBetween = reshape(R_perPairBoth(1:(N0),(N0 + 1):(N0*2)),(N0).^2,1);
    R_perPairWithin = [R_perPairWithinA; R_perPairWithinB];
    R_perPairBetweenPad = padarray(R_perPairBetween,size(R_perPairBetween),NaN,'post');
    R_perPairBetweenPad = R_perPairBetweenPad(:,1);
    withinMinusBetween = nanmean(R_perPairWithin) - nanmean(R_perPairBetween);
    %R_withinMinusBetweenAll(R_Jidx,sigma0idx) = nanmean([R_withinMinusBetweenAll(R_Jidx,sigma0idx),withinMinusBetween]);
    %R_WithinAll(R_Jidx,sigma0idx) = nanmean([R_WithinAll(R_Jidx,sigma0idx),nanmean(R_perPairWithin)]);
    
    R_perPairWithinHistAll = [R_perPairWithinHistAll; R_perPairWithin];
    R_perPairBetweenHistAll = [R_perPairBetweenHistAll; R_perPairBetweenPad]; 
    
    %figure;
    %imagesc(R_perPairBoth);
    R_perPairBoth = reshape(R_perPairBoth,size(R_perPairBoth,1).^2,1);
    R_perPairBoth_all = [R_perPairBoth_all; R_perPairBoth];
    
    
    
end



R_perPairWithinHist_TwSurf(:,iTw) = R_perPairWithinHistAll;
R_perPairBetweenHist_TwSurf(:,iTw) = R_perPairBetweenHistAll;
VarAll_TwSurf(:,iTw) = VarBothAll;
re_TwSurf(:,iTw) = re_all';


end

figure;

subplot(2,1,1);
hist2Stair(R_perPairBetweenHist_TwSurf(:,1),30,'r')
hold on
hist2Stair(R_perPairWithinHist_TwSurf(:,1),30,'k')
hold on
plot([nanmean(R_perPairWithinHist_TwSurf(:,1)) nanmean(R_perPairWithinHist_TwSurf(:,1))],[0 14e4],'--k','linewidth',3)
plot([nanmean(R_perPairBetweenHist_TwSurf(:,1)) nanmean(R_perPairBetweenHist_TwSurf(:,1))],[0 14e4],'--r','linewidth',3)

subplot(2,1,2);
hist2Stair(R_perPairBetweenHist_TwSurf(:,end),30,'r')
hold on
hist2Stair(R_perPairWithinHist_TwSurf(:,end),30,'k')
hold on
plot([nanmean(R_perPairWithinHist_TwSurf(:,end)) nanmean(R_perPairWithinHist_TwSurf(:,end))],[0 14e4],'--k','linewidth',3)
plot([nanmean(R_perPairBetweenHist_TwSurf(:,end)) nanmean(R_perPairBetweenHist_TwSurf(:,end))],[0 14e4],'--r','linewidth',3)

savefig('withinVsBetweenCorr.fig')










