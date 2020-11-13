%% Fit the ploynomial FI curve

% pmu is p(mu), fitted FI curve as function of mu_sim
% Replot pmu_test = pmu to check

 FIFig = figure;
 gainFig = figure;
 %everyNeuronGainVect = NaN(5000,1); %make this not hard coded
 
 everyNeuronGainVect = [];
 
 state = 'A';
 
 switch state
     case 'A'
         titleString = "State A";
         varString = "State A";
 end
 
 
 
 

for i=1:4

    switch i
        case 1 
            currentMat = currentEL_StateA_All;
            frMat = FR_EL_StateA_All./(1000/Tw_bin);
            titleText = "E_L State A";
            gainIndsForWeightMat = positiveSetEL;
        case 3
            currentMat = currentIL_StateA_All;
            frMat = FR_IL_StateA_All./(1000/Tw_bin);
            titleText = "I_L State A";
            gainIndsForWeightMat = positiveSetIL;
        case 2 
            currentMat = currentER_StateA_All;
            frMat = FR_ER_StateA_All./(1000/Tw_bin);
            titleText = "E_R State A";
            gainIndsForWeightMat = positiveSetER;

            
        case 4
            currentMat = currentIR_StateA_All;
            frMat = FR_IR_StateA_All./(1000/Tw_bin);
            titleText = "I_R State A";
            gainIndsForWeightMat = positiveSetIR;

    end
            

    mu_sim = mean(currentMat');
    rates_sim = mean(frMat');


    %p_mu = polyfit(mu_sim,rates_sim,7); 
    %pmu = polyval(p_mu,mu_sim);
    
    
    set(0,'CurrentFigure',FIFig);
    subplot(2,2,i);
    
    slm = slmengine(mu_sim,rates_sim,'plot','off','knots',3,'increasing','on');
    pmu = slmeval(mu_sim,slm);
    

    %figure;
    %scatterhist(mu_sim,rates_sim,'Kernel','on','Location','NorthWest','Direction','out','Color','k','LineStyle',{'-'},'LineWidth',2,'Marker','.','MarkerSize',6);
    plot(mu_sim,rates_sim,'.k');
    hold on
    plot(mu_sim,pmu,'r.')
    title(titleText);
    xlabel('I');
    ylabel('F');










dp_test = slmeval(mu_sim,slm,1);

% check pmu 
% pmu_test = zeros(1,length(mu_sim));
% dp_test = zeros(1,length(mu_sim));
% for k = 1:length(mu_sim)
%     pmu_test(k) = Peval(mu_sim(k),p_mu);
%     dp_test(k) = dPdx(mu_sim(k),p_mu);
% end

%figure;
%plot(mu_sim,pmu_test,'g.')
%xlabel('I');
%ylabel('plotfitResult');

%checkBadFI = (mu_sim > -0.8) + (rates_sim < 0.2);
%dp_test(checkBadFI==2)=0;

% plot gain as function of mu
set(0,'CurrentFigure',gainFig);
subplot(2,2,i);
plot(mu_sim,dp_test,'.')
title(titleText);
xlabel('I');
ylabel('gain');




everyNeuronGainVect = [everyNeuronGainVect dp_test];




end


everyNeuronGainVectPad = [everyNeuronGainVect zeros(1,4000)];

%idxPosIntoWeights = [positiveSetEL positiveSetER positiveSetIL positiveSetIR];
%idxPosIntoWeightsPad = [idxPosIntoWeights zeros(1,4000)];

%weightMat = load('weights_NoCouple_sig00_7_1_tau0_60_freeze1(1)_1_20_20.csv');
%weightMat = h5read('weights_09_08_2020_freeze8a1ed1d0-f1bd-11ea-03cb-3d2ca16203c0.h5','/weights');

weightMatStr = 'Weights_.*.h5';
fnamesH5 = dir(fullfile(currentDir,'*.h5'));
fnamesH5 = {fnamesH5.name}.';
FIND_H5 = @(str)cellfun(@(c)~isempty(c),regexp(fnamesH5,str,'once'));
fnames_subsetWeight = fnamesH5(FIND_H5(weightMatStr));

weightMat = h5read(fnames_subsetWeight{1},'/weights');

%weightMat = weightMat(1:5000,1:5000);
%weightMat = weightMat([idxPosIntoWeights,5001:end],[idxPosIntoWeights,5001:end]);
%hjhjk
weightMatRecurrent = weightMat(1:5000,1:5000);
%weightMat = weightMat(idxPosIntoWeights,idxPosIntoWeights);

weightMatInput = weightMat(1:5000,5001:end);



 for i = 1:5000%length(idxPosIntoWeights) %Pad
     for j=1:5000%length(idxPosIntoWeights) %Pad
         if i < 4001
            synScale = 0.25; %0.25;
         else
            synScale = 0.15; %0.15;
         end
         GainMat(i,j)=weightMatRecurrent(i,j)*everyNeuronGainVect(i)*synScale; 
         if j < (size(weightMatInput,2) + 1)
            GainMatInput(i,j) = weightMatInput(i,j)*everyNeuronGainVect(i)*synScale;
         end
     end
 end
 
 
 
eigCloud = eig(GainMat);
%eigCloudInput = eig(GainMatInput);

figure;
plot(eigCloud,'.k');


stability = eigs(GainMat, 1, 'lm'); %spectral radius calc

covInTheoryBlock = covEstIn .* ones(2000,2000);
covInTheoryBlock = covInTheoryBlock + diag(varEstIn*ones(2000,1) - diag(covInTheoryBlock));
covInTheoryFull = blkdiag(covInTheoryBlock,covInTheoryBlock);





%blockEL = ones(2000,2000);
%blockER = ones(2000,2000);
%blockIL = ones(500,2000);
%blockIR = ones(500,2000);

covSimState = covStateA;

covY0 = GainMatInput*covInTheoryFull*GainMatInput';
%covY0 = GainMatInput*covInput*GainMatInput';

covY0VarReplacement = covY0 + diag(diag(covSimState)- diag(covY0));

 figure;
 subplot(2,2,2);
 imagesc(weightMatInput);
 title('Input Weights');
 subplot(2,2,3);
 imagesc(GainMatInput);
 title('G_{Input}')
 subplot(2,2,1);
 imagesc(covInTheoryFull);
 title('V');
 subplot(2,2,4);
 imagesc(covY0VarReplacement);
 title('Estimated Cov(Y0,Y0)');


covY = inv(eye(5000)-GainMat) * covY0VarReplacement * inv(eye(5000)-GainMat');
covY = covY + diag(NaN(size(covY,1),1) - diag(covY));

LHS = (eye(5000)-GainMat) * covSimState * (eye(5000)-GainMat);
RHS = GainMat * covY0VarReplacement * GainMat';


covSimState = covSimState + diag(NaN(size(covSimState,1),1) - diag(covSimState));
%covY0 = covY0 + diag(NaN(size(covY0,1),1) - diag(covY0));

figure;
subplot(2,2,1);
imagesc(weightMatRecurrent);
title('Recurrent Weights');
subplot(2,2,2);
imagesc(GainMat);
title('G_{R}');
subplot(2,2,4);
imagesc(covSimState);
cSim = caxis;
title('Sim State Cov, All');
subplot(2,2,3);
imagesc(covY);
caxis(cSim);
title('Cov Est, All');

figure;

covSimStateTriU = covSimState; %triu(covSimState,1);
%covSimStateTriU(covSimStateTriU==0)=NaN;
simState_EL_EL = covSimStateTriU(1:2000,1:2000);
subplot(4,2,1);
imagesc(simState_EL_EL(10:30,10:30));
caxis(cSim);
title('EL:EL Subset Sim')
simState_EL_EL = simState_EL_EL(:);
simState_ER_ER = covSimStateTriU(2001:4000,2001:4000);
subplot(4,2,3);
imagesc(simState_ER_ER(10:30,10:30));
caxis(cSim);
title('ER:ER Subset Sim')
simState_ER_ER = simState_ER_ER(:);
simState_ER_EL = covSimStateTriU(2001:4000,1:2000);
subplot(4,2,5);
imagesc(simState_ER_EL(10:30,10:30));
caxis(cSim);
title('ER:EL Subset Sim')
simState_ER_EL = simState_ER_EL(:);
%randER_EL = normrnd( 0 , 0.001, size(simState_ER_EL));
%simState_ER_EL = randER_EL;
simState_EL_ER = covSimStateTriU(1:2000,2001:4000);
subplot(4,2,7);
imagesc(simState_EL_ER(10:30,10:30));
caxis(cSim);
title('EL:ER Subset Sim')
simState_EL_ER = simState_EL_ER(:);
%randEL_ER = normrnd( 0 , 0.001, size(simState_ER_EL));
%simState_EL_ER = randER_EL;

simState_IL_IL = covSimStateTriU(4001:4500,4001:4500);
simState_IL_IL = simState_IL_IL(:);
simState_IR_IR = covSimStateTriU(4501:5000,4501:5000);
simState_IR_IR = simState_IR_IR(:);

covEstTriU = covY; %triu(covY,1);
%covEstTriU(covY==0)=NaN;
covEst_EL_EL = covEstTriU(1:2000,1:2000);
subplot(4,2,2);
imagesc(covEst_EL_EL(10:30,10:30));
caxis(cSim);
title('EL:EL Subset Theory')
covEst_EL_EL = covEst_EL_EL(:);
covEst_ER_ER = covEstTriU(2001:4000,2001:4000);
subplot(4,2,4);
imagesc(covEst_ER_ER(10:30,10:30));
caxis(cSim);
title('ER:ER Subset Theory')
covEst_ER_ER = covEst_ER_ER(:);
covEst_ER_EL = covEstTriU(2001:4000,1:2000);
subplot(4,2,6);
imagesc(covEst_ER_EL(10:30,10:30))
caxis(cSim);
title('ER:EL Subset Theory')
covEst_ER_EL = covEst_ER_EL(:);
covEst_EL_ER = covEstTriU(1:2000,2001:4000);
subplot(4,2,8);
imagesc(covEst_EL_ER(10:30,10:30));
caxis(cSim);
title('EL:ER Subset Theory')
covEst_EL_ER = covEst_EL_ER(:);

covEst_IL_IL = covEstTriU(4001:4500,4001:4500);
covEst_IL_IL = covEst_IL_IL(:);
covEst_IR_IR = covEstTriU(4501:5000,4501:5000);
covEst_IR_IR = covEst_IR_IR(:);

perEr_All = (covY - covSimState)./covY;

perEr_EL_EL = (covEst_EL_EL*1000 - simState_EL_EL*1000)./(simState_EL_EL*1000);
%perEr_EL_EL_NoOut = perEr_EL_EL(perEr_EL_EL > -5);
%perEr_EL_EL_NoOut = perEr_EL_EL_NoOut(perEr_EL_EL_NoOut < 5);

perEr_ER_ER = (covEst_ER_ER*1000 - simState_ER_ER*1000)./(simState_ER_ER*1000);
%perEr_ER_ER_NoOut = perEr_ER_ER(perEr_ER_ER > -5);
%perEr_ER_ER_NoOut = perEr_ER_ER_NoOut(perEr_ER_ER_NoOut < 5);

perEr_IL_IL = (covEst_EL_ER*1000 - simState_EL_ER*1000)./(simState_EL_ER*1000);
%perEr_IL_IL_NoOut = perEr_IL_IL(perEr_IL_IL > -5);
%perEr_IL_IL_NoOut = perEr_IL_IL_NoOut(perEr_IL_IL_NoOut < 5);


perEr_IR_IR = (covEst_ER_EL*1000 - simState_ER_EL*1000)./(simState_ER_EL*1000);
%perEr_IR_IR_NoOut = perEr_IR_IR(perEr_IR_IR > -5);
%perEr_IR_IR_NoOut = perEr_IR_IR_NoOut(perEr_IR_IR_NoOut < 5);



figure;
subplot(2,2,1);
[n,x] = hist2Stair(simState_EL_EL,40,'k');
hold on;
[n,x] = hist2Stair(covEst_EL_EL,x,'r');
hold on;
yLims = get(gca,'YLim');
plot([nanmean(simState_EL_EL) nanmean(simState_EL_EL)],[0 yLims(2)],'k--', 'linewidth',3);
hold on;
plot([nanmean(covEst_EL_EL) nanmean(covEst_EL_EL)], [0 yLims(2)], 'r--','linewidth',3);
title('EL:EL Sim Black Theory Red')

subplot(2,2,2);
[n,x] = hist2Stair(simState_ER_ER,40,'k');
hold on;
[n,x] = hist2Stair(covEst_ER_ER,x,'r');
hold on;
yLims = get(gca,'YLim');
plot([nanmean(simState_ER_ER) nanmean(simState_ER_ER)],[0 yLims(2)],'k--', 'linewidth',3);
hold on;
plot([nanmean(covEst_ER_ER) nanmean(covEst_ER_ER)], [0 yLims(2)], 'r--','linewidth',3);
title('ER:ER Sim Black Theory Red')


subplot(2,2,3);
[n,x] = hist2Stair(simState_EL_ER,40,'k');
hold on;
[n,x] = hist2Stair(covEst_EL_ER,x,'r');
hold on;
yLims = get(gca,'YLim');
plot([nanmean(simState_EL_ER) nanmean(simState_EL_ER)],[0 yLims(2)],'k--', 'linewidth',3);
hold on;
plot([nanmean(covEst_EL_ER) nanmean(covEst_EL_ER)], [0 yLims(2)], 'r--','linewidth',3);
title('EL:ER Sim Black Theory Red')



subplot(2,2,4);
[n,x] = hist2Stair(simState_ER_EL,40,'k');
hold on;
[n,x] = hist2Stair(covEst_ER_EL,x,'r');
hold on;
yLims = get(gca,'YLim');
plot([nanmean(simState_ER_EL) nanmean(simState_ER_EL)],[0 yLims(2)],'k--', 'linewidth',3);
hold on;
plot([nanmean(covEst_ER_EL) nanmean(covEst_ER_EL)], [0 yLims(2)], 'r--','linewidth',3);
title('ER:EL Sim Black Theory Red')



figure;
subplot(2,2,1);
[n,x] = hist2Stair(perEr_EL_EL,[-5:0.05:5],'k');
yLims = get(gca,'YLim');
hold on;
plot([nanmean(perEr_EL_EL) nanmean(perEr_EL_EL)],[0 yLims(2)],'k--', 'linewidth',3);
title('EL:EL Error')

subplot(2,2,2);
[n,x] = hist2Stair(perEr_ER_ER,[-5:0.05:5],'k');
yLims = get(gca,'YLim');
hold on;
plot([nanmean(perEr_ER_ER) nanmean(perEr_ER_ER)],[0 yLims(2)],'k--', 'linewidth',3);
title('ER:ER Error')


subplot(2,2,3);
[n,x] = hist2Stair(perEr_IL_IL,[-5:0.05:5],'k');
yLims = get(gca,'YLim');
hold on;
plot([nanmean(perEr_IL_IL) nanmean(perEr_IL_IL)],[0 yLims(2)],'k--', 'linewidth',3);
title('EL:ER Error')


subplot(2,2,4);
[n,x] = hist2Stair(perEr_IR_IR,[-5:0.05:5],'k');
yLims = get(gca,'YLim');
hold on;
plot([nanmean(perEr_IR_IR) nanmean(perEr_IR_IR)],[0 yLims(2)],'k--', 'linewidth',3);
title('ER:EL Error')
























%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function output = GainGenerator2(x,y,orderTerm)
p = polyfit(x,y,orderTerm); 
f = polyval(p,x);
TOL = 0.01;
errorvect = y-f;
maxerror = max(errorvect);
disp(maxerror)
% if maxerror > TOL
%     disp('Not a good fit')
% end


%temp22 = polyval(p,x);
%figure;plot(x,temp22,'.');
% Compute the derivative

k = polyder(p)
output = polyval(k,x);
end





function y = Peval(x,p)
    y = 0;
    for k = 1:8
        y = y + p(k) * x^(8-k);
    end
    %y = p(1)*x^7 + p(2)*x^6 +p(3)*x^5 +p(4)*x^4 +p(5)*x^3 +p(6)*x^2 +p(7)*x +p(8);
end




function y = dPdx(x,p)
    y = 0;
    for k = 1:7
        y = y + p(k) * (8-k) * x^(7-k);
    end
end





function g_mu = gamma_mu(params)

    dmu = 0.01;

    %params_th1 = [E0, mu0+dmu, tau0, Vth, Vre, tref, sigma0, sigma2];
    params_th1 = params;
    params_th1(2) = params(2) + dmu;
    params_th0 = params;

    
    [r1,~] = test_20190322_LIFth(params_th1);
    [r0,~] = test_20190322_LIFth(params_th0);
    
    g_mu = (r1-r0)/dmu;


end











