figure;
covEst_Reg = [covEst_EL_EL; covEst_ER_ER; covEst_EL_ER];
simState_Reg = [simState_EL_EL; simState_ER_ER;simState_EL_ER];
notNaN = ~isnan(covEst_Reg);
linFit = polyfit(simState_Reg(notNaN),covEst_Reg(notNaN),1);
[P,S] = polyfit(simState_Reg(notNaN),covEst_Reg(notNaN),1);
y_est = polyval(linFit,simState_Reg(notNaN));

C_simVC_theoryFigAutogen(simState_ER_ER, covEst_ER_ER, simState_ER_ER, covEst_ER_ER, simState_EL_ER,covEst_EL_ER,simState_EL_EL,covEst_EL_EL, simState_Reg(notNaN),y_est);

y = covEst_Reg(notNaN);
RSq = 1 - (S.normr/norm(y - mean(y)))^2

figure;
%compare to shuffled data
simState_RegNew = simState_Reg(notNaN);
covEst_RegNew = covEst_Reg(notNaN);
rand_pos = randperm(length(simState_RegNew)); %array of random positions
% new array with original data randomly distributed 
for k = 1:length(simState_RegNew)
    simState_shuff(k) = simState_RegNew(rand_pos(k));
end
simState_shuff = simState_shuff';
linFit = polyfit(simState_shuff,covEst_RegNew,1);
y_est = polyval(linFit,simState_shuff);

C_simVC_theoryFigShuff(simState_shuff,covEst_RegNew,simState_shuff,y_est);


[P,S] = polyfit(simState_shuff,covEst_RegNew,1);
y = covEst_RegNew;
Rsq_shuff = 1 - (S.normr/norm(y - mean(y)))^2