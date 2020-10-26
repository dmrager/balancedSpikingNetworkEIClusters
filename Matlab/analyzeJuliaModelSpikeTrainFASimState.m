segmentFA = 0; %change to 0 for linear response!

figure;

Ne = 2000;
Nc = 500;
Nc_FA = 100;
N0 = Nc;
TW_bin = 50;
T = 12000; %WAS 10000!!!
Tburn = 500; %WAS 0!!!!
binedges=Tburn+1:TW_bin:T; %was Tburn not +1?

searchStr = 'spHemi_.*.csv';
searchStrNs = 'nsHemi_.*.csv';
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
countsAllFullPop = [];


simRefArray = zeros(3,length(fnames_subset));
R_JSweep = 1:0.2:2.4;
sigma0Sweep = 0:0.4:2.4;

R_withinMinusBetweenAll = NaN(length(R_JSweep),length(sigma0Sweep));
R_WithinAll = NaN(length(R_JSweep),length(sigma0Sweep));

R_perPairBetweenHistAll = [];
R_perPairWithinHistAll = [];

s1Struct(length(fnames_subset)).s1 = NaN(2,300000);

if segmentFA == 1
    simsUsed = 1:225;
else
    simsUsed = 1:length(fnames_subset);
end

parfor iSimCount = 1:length(simsUsed)%length(fnames_subset)
        iSim = simsUsed(iSimCount);
    
         reBoth = [];
    
        fileStr = fnames_subset{iSim};
        fileStrNs = fnames_subsetNs{iSim};
        %simNumRegex = regexp(fileStr,'(^spHemi_)(\d*)(.csv)','tokens');
        %simNum = str2double(simNumRegex{1}(2));
       % sigma0idx = ceil(simNum/24);
        %R_Jidx = ceil(mod(ceil(simNum/3),8));
        %if R_Jidx == 0
            R_Jidx = 8;
       % end
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
    
    
    s1 = [timesAll;nNumAll];
    s1 = s1(:,timesAll<=T);
    s1Struct(iSimCount).s1 = s1;
    
end
    
 for iSimCount = 1:length(simsUsed)%length(fnames_subset)
     
    s1 = s1Struct(iSimCount).s1; 
     
    if iSimCount == 1
        I1=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
    
        I1=I1(I1<=(Ne*2));
        Ic1=randsample(I1(I1<=Ne),Nc_FA/2); %why is this 1000? should be Ne? 
        Ic2=randsample(I1(I1>Ne),Nc_FA/2);
        Icall = [Ic1;Ic2];
    end
    
    cat = 3;
    
        
%     for mm=1:Nc_FA
%         counts(mm,:)=histc( s1(1, Icall(mm)-1/4<s1(2,:) & s1(2,:)<=Icall(mm)+1/4) ,binedges);
%         %counts2(mm,:)=histc(s1(1,s1(2,:)==Icall(mm)),binedges);
%         cat = 3;
%         
%     end
    
     parfor mm2=1:length(I1)
          countsFullPop(mm2,:)=histc( s1(1, I1(mm2)-1/4<s1(2,:) & s1(2,:)<=I1(mm2)+1/4) ,binedges);
     end
    
    %countsAll = [countsAll, counts];
    countsAllFullPop = [countsAllFullPop, countsFullPop];
    
    cat = 3;
    
    
    
        
%         timesAll = []; nNumAll = [];

%         
%         for ci = ((iCluster * Ne) + 1):(Ne * (iCluster + 1))
%             timevals = times(ci,1:ns(ci)); nNum = ci *
%             ones(1,length(timevals)); timesAll = [timesAll, timevals];
%             nNumAll = [nNumAll,nNum];
%         end
%         
%         T = 10000;
%         time=0:1:T;
%         Tw=200;
%         time=0:Tw:T; Tw_FAprep=150;
%         Tburn = 0;
%         
%         timesAll = sxChunkBothHemi(1,:);
%         
%         s1 = sxChunkBothHemi;
%         
%         s1 = s1(:,timesAll<=10000);
%         
%         I1=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
%         
%         I1=I1(I1<=(Ne*2));
%         
%         if iCluster == 0
%             Ic1=randsample(I1(I1<=Ne),500);
%         else
%             Ic1=randsample(I1(I1>Ne),500);
%         end
%         
%         cat = 3;
%         
%         for mm=1:Nc
%             re1(mm,:)=hist(s1(1,Ic1(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1(mm)+1/4),time)*1e3;
%         end
%         
%         re1_s=imfilter(re1(:,Tburn+1:end),ones(1,Tw)/Tw);re1_s=re1_s(:,Tw/2-1:end-Tw/2);
%         plot(mean(re1_s));
%         xlabel('Time (ms)');
%         ylabel('Mean FR (Hz)');
%         hold all;
%         
%         reBoth = [reBoth; re1_s];
%         
%         re1_s_reshape = reshape(mean(re1_s'),size(mean(re1_s'),1)*size(mean(re1_s'),2),1);
%         re1_s_reshape = re1_s_reshape';
%         re_all = [re_all, re1_s_reshape];
%         
%         
%         COV=cov(re1_s');
%         Var=diag(COV);
%         rate1=mean(re1_s,2);
%         
%         R = COV./sqrt(Var*Var');
%         R_perPair = triu(R,1);
%         R_perPair(R_perPair==0)=NaN;
%         
%         
%         R_perPair_Within = R_perPair(1:(N0/2),1:(N0/2))
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
%         R_perPair_allA =
%         reshape(R_perPair_allA,size(R_perPair_allA,1).^2,1);
%         R_perPair_allB =
%         reshape(R_perPair_allB,size(R_perPair_allB,1).^2,1);
% 
%         
%         meanR((iCluster+1)) = nanmean(R_perPair);
%         
%     end
%     
%     COVBoth=cov(reBoth');
%     VarBoth=diag(COVBoth);
%     rate1=mean(re1_s,2);
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
%     R_withinMinusBetweenAll(R_Jidx,sigma0idx) = nanmean([R_withinMinusBetweenAll(R_Jidx,sigma0idx),withinMinusBetween]);
%     R_WithinAll(R_Jidx,sigma0idx) = nanmean([R_WithinAll(R_Jidx,sigma0idx),nanmean(R_perPairWithin)]);
%     
%     R_perPairWithinHistAll = [R_perPairWithinHistAll; R_perPairWithin];
%     R_perPairBetweenHistAll = [R_perPairBetweenHistAll; R_perPairBetween];
%     
%     figure; imagesc(R_perPairBoth);
%     R_perPairBoth = reshape(R_perPairBoth,size(R_perPairBoth,1).^2,1);
%     R_perPairBoth_all = [R_perPairBoth_all; R_perPairBoth];
    
    
    
end



 %Good params for 50 ms
 %[row, col] = find(mean(countsAll(1:(Nc_FA/2),:))>0.2 & mean(countsAll((Nc_FA/2+1):Nc_FA,:))<0.2 & mean(countsAll((Nc_FA/2+1):Nc_FA,:))>0.08);
 %[row2, col2] = find(mean(countsAll((Nc_FA/2+1):Nc_FA,:))>0.2 & mean(countsAll(1:(Nc_FA/2),:))<0.2 & mean(countsAll(1:(Nc_FA/2),:))>0.08);
 
 %[row, col] = find(mean(countsAll(1:(Nc_FA/2),:))>0.15 & mean(countsAll((Nc_FA/2+1):Nc_FA,:))<0.15 & mean(countsAll((Nc_FA/2+1):Nc_FA,:))>0.02);
 %[row2, col2] = find(mean(countsAll((Nc_FA/2+1):Nc_FA,:))>0.15 & mean(countsAll(1:(Nc_FA/2),:))<0.15 & mean(countsAll(1:(Nc_FA/2),:))>0.02);
 
 %meanLMinus = mean(countsAll(1:(Nc_FA/2),:)) - 0.25*std(countsAll(1:(Nc_FA/2),:));%/length(mean(countsAll(1:(Nc_FA/2),:)));
 %meanLPlus = mean(countsAll(1:(Nc_FA/2),:)) + 0.25*std(countsAll(1:(Nc_FA/2),:));%/length(mean(countsAll(1:(Nc_FA/2),:)));
 %meanRMinus = mean(countsAll((Nc_FA/2+1):Nc_FA,:)) - 0.25*std(countsAll((Nc_FA/2+1):Nc_FA,:));%/length(mean(countsAll((Nc_FA/2+1):Nc_FA,:)));
 %meanRPlus = mean(countsAll((Nc_FA/2+1):Nc_FA,:)) + 0.25*std(countsAll((Nc_FA/2+1):Nc_FA,:));%/length(mean(countsAll((Nc_FA/2+1):Nc_FA,:)));
 
%[row, col] = find(meanLMinus>meanRPlus);
%[row2, col2] = find(meanRMinus>meanLPlus);




 
 %countsAllL = countsAll(:,col);
 %countsAllR = countsAll(:,col2);
 






%  figure;
%  subplot(3,1,2);
%  plot(timesAll,nNumAll,'.k');
%  ylabel('Neuron Idx','fontsize',16);
%  axis([1 1000 1 4000]);
%  subplot(3,1,3);
%  plot(mean(countsAll(1:(Nc_FA/2),:)),'b');
%  hold on;
%  plot(mean(countsAll((Nc_FA/2+1):Nc_FA,:)),'r');
%  hold on;
%  plot(mean(countsAllFullPop(1:2000,:)),'k');
%  hold on;
%  plot(mean(countsAllFullPop(2001:end,:)),'m');
%    
%  %plot([round(4000/TW_bin) round(6000/TW_bin)],[0.5 0.5],'--k');
%  %plot([round(4000/TW_bin) round(6000/TW_bin)],[0.1 0.1],'--k');
%  plot(find(idx==1),ones(size(find(idx==1))),'.g');
%  plot(find(idx==2),ones(size(find(idx==2))),'.m');
%  axis([round(1/TW_bin) round(1000/TW_bin) 0 1.2])
%  subplot(3,1,1);
%  plot(P(:,1),'.k')
%  axis([round(1/TW_bin) round(1000/TW_bin) 0 1.1])

 
 
 
 