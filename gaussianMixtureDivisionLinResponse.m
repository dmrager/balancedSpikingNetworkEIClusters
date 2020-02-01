FR_EL_StateA_All(:,idxLowActivity)=[];
FR_ER_StateA_All(:,idxLowActivity)=[];
FR_IL_StateA_All(:,idxLowActivity)=[];
FR_IR_StateA_All(:,idxLowActivity)=[];

currentEL_StateA_All(:,idxLowActivity)=[];
currentER_StateA_All(:,idxLowActivity)=[];
currentIL_StateA_All(:,idxLowActivity)=[];
currentIR_StateA_All(:,idxLowActivity)=[];

FR_EL_StateB_All=FR_EL_StateA_All(:,cluster2Good);
FR_ER_StateB_All=FR_ER_StateA_All(:,cluster2Good);
FR_IL_StateB_All=FR_IL_StateA_All(:,cluster2Good);
FR_IR_StateB_All=FR_IR_StateA_All(:,cluster2Good);

currentEL_StateB_All=currentEL_StateA_All(:,cluster2Good);
currentER_StateB_All=currentER_StateA_All(:,cluster2Good);
currentIL_StateB_All=currentIL_StateA_All(:,cluster2Good);
currentIR_StateB_All=currentIR_StateA_All(:,cluster2Good);

FR_EL_StateA_All=FR_EL_StateA_All(:,cluster1Good);
FR_ER_StateA_All=FR_ER_StateA_All(:,cluster1Good);
FR_IL_StateA_All=FR_IL_StateA_All(:,cluster1Good);
FR_IR_StateA_All=FR_IR_StateA_All(:,cluster1Good);

currentEL_StateA_All=currentEL_StateA_All(:,cluster1Good);
currentER_StateA_All=currentER_StateA_All(:,cluster1Good);
currentIL_StateA_All=currentIL_StateA_All(:,cluster1Good);
currentIR_StateA_All=currentIR_StateA_All(:,cluster1Good);

StateBAll = [FR_EL_StateB_All; FR_ER_StateB_All; FR_IL_StateB_All; FR_IR_StateB_All];
StateAAll = [FR_EL_StateA_All; FR_ER_StateA_All; FR_IL_StateA_All; FR_IR_StateA_All];

StateAAll = StateAAll./(1000/Tw_bin);
StateBAll = StateBAll./(1000/Tw_bin);

covStateA = cov(StateAAll');
covStateB = cov(StateBAll');



















