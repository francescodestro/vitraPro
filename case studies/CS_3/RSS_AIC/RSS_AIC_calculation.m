%% Model by Destro and Braatz, 2024) 
clc, clear, close all

load('simulation_results.mat') % simulation output - model by Destro and Braatz (2024)

% Load experimental data from:
%   Rüdiger, D., Pelz, L., Hein, M.D., Kupke, S.Y. and Reichl, U., 2021. 
%   Multiscale model of defective interfering particle replication for 
%   influenza A virus infection in animal cell culture. PLoS Computational 
%   Biology, 17(9), p.e1009357.
load influenza_data.mat 

%  RSS calculation 
start_point=0.01; % neglect initial condition for RSS calculation

for i=1:12

    % STV
    STV_exp=data{i}.STV_PFU;
    t_exp=data{i}.t; 
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(STV_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(STV_exp));
    end
    STV_exp=STV_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    STV_sim=STV_PFU_sim{i};
    STV_pred=interp1(t_sim{i},STV_sim,t_exp);
    scaling_factor_per_exp=max(STV_exp);
    RSS(i,1)=sum(((STV_pred-STV_exp)/scaling_factor_per_exp).^2);

    % DIP
    DIP_exp=data{i}.DIP;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(DIP_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(DIP_exp));
    end
    DIP_exp=DIP_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    DIP_sim_RSS=DIP_sim{i};
    DIP_pred=interp1(t_sim{i},DIP_sim_RSS,t_exp);
    scaling_factor_per_exp=max(DIP_exp);
    RSS(i,2)=sum(((DIP_pred-DIP_exp)/scaling_factor_per_exp).^2);

    % STV vRNA
    STV_vRNA_exp=data{i}.FL_vRNA;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(STV_vRNA_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(STV_vRNA_exp));
    end
    STV_vRNA_exp=STV_vRNA_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    DNA1_sim=DNA1_v_sim{i}+DNA1_nv_sim{i};
    STV_vRNA_sim=DNA1_sim./(I1_sim{i}+I2_sim{i}+Co_sim{i}+T_sim{i}+NV_sim{i});
    STV_vRNA_pred=interp1(t_sim{i},STV_vRNA_sim,t_exp);
    scaling_factor_per_exp=max(STV_vRNA_exp);
    RSS(i,3)=sum(((STV_vRNA_pred-STV_vRNA_exp)/scaling_factor_per_exp).^2);

    % DIP vRNA
    DIP_vRNA_exp=data{i}.DI_vRNA;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(DIP_vRNA_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(DIP_vRNA_exp));
    end
    DIP_vRNA_exp=DIP_vRNA_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    DIP_vRNA_sim=DNA2_sim{i}./(I1_sim{i}+I2_sim{i}+Co_sim{i}+T_sim{i}+NV_sim{i});
    DIP_vRNA_pred=interp1(t_sim{i},DIP_vRNA_sim,t_exp);
    scaling_factor_per_exp=max(DIP_vRNA_exp);
    RSS(i,4)=sum(((DIP_vRNA_pred-DIP_vRNA_exp)/scaling_factor_per_exp).^2);

    % VCC
    VCC_exp=data{i}.VCC;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(VCC_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(VCC_exp));
    end
    VCC_exp=VCC_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    VCC_sim=T_sim{i}+NV_sim{i}+I1_sim{i}+I2_sim{i}+Co_sim{i};
    VCC_pred=interp1(t_sim{i},VCC_sim,t_exp);
    scaling_factor_per_exp=max(VCC_exp);
    RSS(i,5)=sum(((VCC_pred-VCC_exp)/scaling_factor_per_exp).^2);

end
RSS_destro=sum(RSS,2);

% AIC calculation
N_samples=[28, 42, 68, 51, 24, 38, 41, 41, 24, 38, 41, 41]';
N_par=11; % number of estimated parameters in Destro and Braatz, 2024 (CS3)
AIC_destro=2*N_par+N_samples.*log(RSS_destro./N_samples);

RSS_cum_destro=sum(RSS_destro);
AIC_cum_destro=2*N_par+sum(N_samples).*log(RSS_cum_destro./sum(N_samples));

%% Model by Rudiger et al., 2021) 
% Load simulation results from:
%   Rüdiger, D., Pelz, L., Hein, M.D., Kupke, S.Y. and Reichl, U., 2021. 
%   Multiscale model of defective interfering particle replication for 
%   influenza A virus infection in animal cell culture. PLoS Computational 
%   Biology, 17(9), p.e1009357.
load('sim_output_rudiger.mat')

%  RSS calculation 
% Retrieve states
SimStatesIn = p.InToPop.states;
SimValuesIn = p.InToPop.populationvalues;
SimValuesIn(isnan(SimValuesIn)) = 0;

SimTimeEx   = result.Ex.time;
SimStatesEx = [result.Ex.states, result.Ex.variables];
SimValuesEx = [result.Ex.statevalues result.Ex.variablevalues];

for i = 1 : length(SimStatesIn)
    SimIndex.(SimStatesIn{i}) = strcmp(SimStatesIn{i}, SimStatesIn);
end
for i = 1 : length(SimStatesEx)
    SimIndex.(SimStatesEx{i}) = strcmp(SimStatesEx{i}, SimStatesEx);
end
for i = 1 : length(d.In.ExpStateNames)
    ExpIndex.(d.In.ExpStateNames{i}) = strcmp(d.In.ExpStateNames{i}, d.In.ExpStateNames);
end
for i = 1 : length(d.Ex.ExpStateNames)
    ExpIndex.(d.Ex.ExpStateNames{i}) = strcmp(d.Ex.ExpStateNames{i}, d.Ex.ExpStateNames);
end

for i=1:12

    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.RvSeg1, i+1)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime);  

    % STV
    STV_exp=data{i}.STV_PFU;
    t_exp=data{i}.t; 
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(STV_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(STV_exp));
    end
    STV_exp=STV_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    STV_pred=interp1(SimTimeEx,SimValuesEx(:, SimIndex.VOnlyRel, i+1),t_exp);
    scaling_factor_per_exp=max(STV_exp);
    RSS(i,1)=sum(((STV_pred-STV_exp)/scaling_factor_per_exp).^2);

    % DIP
    DIP_exp=data{i}.DIP;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(DIP_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(DIP_exp));
    end
    DIP_exp=DIP_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    DIP_pred=interp1(SimTimeEx,SimValuesEx(:, SimIndex.DTot, i+1),t_exp); 
    scaling_factor_per_exp=max(DIP_exp);
    RSS(i,2)=sum(((DIP_pred-DIP_exp)/scaling_factor_per_exp).^2);

    % STV vRNA
    STV_vRNA_exp=data{i}.FL_vRNA;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(STV_vRNA_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(STV_vRNA_exp));
    end
    STV_vRNA_exp=STV_vRNA_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    t_exp(end)=round(t_exp(end));
    STV_vRNA_pred=interp1(SimPlotTime,SimValuesIn(1:length(SimPlotTime),...
        SimIndex.RvSeg1,i+1),t_exp);
    scaling_factor_per_exp=max(STV_vRNA_exp);
    RSS(i,3)=sum(((STV_vRNA_pred-STV_vRNA_exp)/scaling_factor_per_exp).^2);

    % DIP vRNA
    DIP_vRNA_exp=data{i}.DI_vRNA;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(DIP_vRNA_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(DIP_vRNA_exp));
    end
    DIP_vRNA_exp=DIP_vRNA_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    if length(t_exp)>1
        t_exp(end)=round(t_exp(end));
    end
    DIP_vRNA_pred=interp1(SimPlotTime,SimValuesIn(1:length(SimPlotTime),...
            SimIndex.RvSeg9,i+1),t_exp);
    scaling_factor_per_exp=max(DIP_vRNA_exp);
    RSS(i,4)=sum(((DIP_vRNA_pred-DIP_vRNA_exp)/scaling_factor_per_exp).^2);

    % VCC
    VCC_exp=data{i}.VCC;
    t_exp=data{i}.t;
    if i==4
        indices_exp=and(and(t_exp>start_point,~isnan(VCC_exp)),t_exp<40);
    else
        indices_exp=and(t_exp>start_point,~isnan(VCC_exp));
    end
    VCC_exp=VCC_exp(indices_exp);
    t_exp=t_exp(indices_exp);
    VCC_sim=SimValuesEx(:, SimIndex.T, i+1)+...
         SimValuesEx(:, SimIndex.I, i+1)+ SimValuesEx(:, SimIndex.Ia, i+1)+...
         SimValuesEx(:, SimIndex.Ta, i+1);
    VCC_pred=interp1(SimTimeEx,VCC_sim,t_exp);
    scaling_factor_per_exp=max(VCC_exp);
    RSS(i,5)=sum(((VCC_pred-VCC_exp)/scaling_factor_per_exp).^2);

end

RSS_rudiger=sum(RSS,2);

% AIC calculation
N_samples=[28, 42, 68, 51, 24, 38, 41, 41, 24, 38, 41, 41]';
N_par=19; % number of estimated parameters in Rudiger et al., 2021
AIC_rudiger=2*N_par+N_samples.*log(RSS_rudiger./N_samples);

RSS_cum_rudiger=sum(RSS_rudiger);
AIC_cum_rudiger=2*N_par+sum(N_samples).*log(RSS_cum_rudiger./sum(N_samples));

