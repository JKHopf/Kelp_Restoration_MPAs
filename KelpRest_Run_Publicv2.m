% KelpRest_Run_Publicv2.m
%
% Jess Hopf (2021) 
%
% This is the main model code for "No-take marine protected areas can
% enhance the fish and fishery outcomes of kelp forest restoration"
% Hopf, Caselle, White 
%
% This calls on a number of function files, the main model projection one
% is 'Func_ProjPop_KR4pop_BHv1.m'
%
% Model overview:
% 4 population model (two MPA/reserve, two fished). 
% Density-dep model: Beverton-Holt (intra-cohort) recrtuiment
% Fishing with fishery squeeze at time of MPA
% assumes global distrubance
% inc kelp restoration in either the MPA or fished area
% Long and short-term dynamics

%-------------------------------------------------------------------------  
    
%% Model

clear
% add file with relevant functions
addpath('.\Functions')

% Variable variables -------------------------------------    

    % number of populations/patches:
    p = 4; 

    % run-times 
    tInit = 500;  
    tPost = 150;
    
    % single run? (use if looking at short-term dynamics, will plot short-term dyanmics)
    single_run = 'n'; % 'y'; %    
    
    % Disturbance timing 
    % 2yrs for long-term runs
    % [] for short-term runs
    % 105yrs if not having disturbance and/or restoration
    yrDist_Vec = 2;% [5,15,50];%  105;%  year dist happens (aka reserve age) 0 = yr MPA est
    rDist_Vec = 2;% [1:5,105];% 105;%  length of disturbance before kelp resotration
   
	% which species
    spp = 'Kelp bass'; %  'California sheephead'; %  
    
        % Species paras:
        load SppParameterVals.mat
        SpParas = SppParas(SppParas.sp == spp,:);
        clear SppParas
    
    % Exploration parameters: 
                                    
        % Proportion of area in reserves
        Area_vec = 0.2; %0;%  

        % Proportion of area restored (needs to be <Area)
        A_r_vec = 0.2; % 0:0.05:0.2; %    

        % magnitude of disturbance (reduction in recruit survival) 
        % for juveniles affected (iDist = 1-gamma_0)
        iDist_vec = 0.5;% 0.05:0.05:0.95;%

        % Who is affected by the disturbance?
        who_dist = 'juvs'; %  'adults'; % 

        % does fish density affect kelp density (and therefore survival?) 
        fishaffect = 'no'; %'yes'; %     

        % which population does the restoration occur in?
        rest_pop = 'rest_fish';  % 'rest_res'; %    'rest_both';% 
        
        % fishing pressure (nan if not varying over fishing pressure)
        mf_vec = nan; %linspace(SpParas.m/2,SpParas.m*2,10); %    

        % override fishing pressure (mf)
            % extremes for Cali SH 
%             SpParas.mf = SpParas.m*2;%  SpParas.m/2;% 
        
    
        
% Model set-up -----------------------------    

% DD Function:
DD_func =  str2func('Func_ProjPop_KR4pop_BHv1'); 

% Beverton-Holt Para values: 
    % a = slope at/near zero/origin
    % calculuate life-time egg production
    Fun = Func_fecunds((1:SpParas.A_max)',...
          Func_Length((1:SpParas.A_max)',SpParas.L_inf,SpParas.K,SpParas.A0)*10,...
                             SpParas.c, SpParas.d);
    Fun(1:(SpParas.A_mat-1))=0;                      
    LEP = sum(cumprod([1;repmat(exp(-SpParas.m), SpParas.A_max-1, 1)]).*Fun);
    % set so that if population drops below 25% of the unfished LEP then
    % population declines
    SpParas.BHa = 1/(0.25*LEP); % 
    
    % b = max density of settlers
    SpParas.BHb = 1000; 

    
% Shape of disturbance multiplier vs. biomass density curve (proxy for fish supporting
% kelp, which increases survival)
% >1 (higher = steeper curve) 
% 1.000000001 is effectively 1 (which means no effect of fish on kelp)
    switch fishaffect
        case 'yes'
            Sg =  [0,logspace(-6,0,40)]; % 1+logspace(-9,-2,40); % 0; %
        case 'no'
            Sg = 0; 
    end

% Variable varying
    Var_vec = nan; 
    
    if length(Sg) >1
    Var_vec = Sg;
    end

    if length(Area_vec) >1
    Var_vec = Area_vec;
    end

    if length(A_r_vec) >1
    Var_vec = A_r_vec;
    end
    
    if length(mf_vec) >1
    Var_vec = mf_vec;
    end
    
    if length(iDist_vec) >1
    Var_vec = iDist_vec;
    end

% Empty matrices
NaN_Mat = NaN(length(yrDist_Vec), length(rDist_Vec));

NDpcM = NaN_Mat; RtimeAM = NaN_Mat; RateAM = NaN_Mat; 
RtimeBM = NaN_Mat; RateBM = NaN_Mat;
NDpcInM = NaN_Mat; RtimeAInM = NaN_Mat; RateAInM = NaN_Mat;
RtimeBInM = NaN_Mat; RateBInM = NaN_Mat;

NDpcR = NaN_Mat; RtimeAR = NaN_Mat; RateAR = NaN_Mat; 
RtimeBR = NaN_Mat; RateBR = NaN_Mat;
NDpcInR = NaN_Mat; RtimeAInR = NaN_Mat; RateAInR = NaN_Mat;
RtimeBInR = NaN_Mat; RateBInR = NaN_Mat;

NDpcF = NaN_Mat; RtimeAF = NaN_Mat; RateAF = NaN_Mat; 
RtimeBF = NaN_Mat; RateBF = NaN_Mat;
NDpcInF = NaN_Mat; RtimeAInF = NaN_Mat; RateAInF = NaN_Mat;
RtimeBInF = NaN_Mat; RateBInF = NaN_Mat;

NM = NaN(length(Var_vec),1); N1 = NM; N2 = NM; N3 = NM; N4 = NM;
YM = NM; Y1 = NM; Y2 = NM; Y3 = NM; Y4 = NM;
NBM = NM; NB1 = NM; NB2 = NM; NB3 = NM; NB4 = NM;
YBM = NM; YB1 = NM; YB2 = NM; YB3 = NM; YB4 = NM;
NNRM = NM; NNR1 = NM; NNR2 = NM; NNR3 = NM; NNR4 = NM;
YNRM = NM; YNR1 = NM; YNR2 = NM; YNR3 = NM; YNR4 = NM;

% Run over variable vector of choice -----------------------------------

for i = 1:length(Var_vec)
    
    % Shape of disturbance multiplier
    SpParas.Sg = Sg(1);
    
    % Proportion of area in reserves
    Area = Area_vec(1);
    
    % Proportion of area restored (needs to be <Area)
    A_r = A_r_vec(1);
    
    % Magnitdue of distrubance
    iDist = iDist_vec(1);
    
    % reassign values if varying over parameter range
    if length(Sg) >1
    SpParas.Sg = Sg(i);
    end
    
    if length(Area_vec) >1
    Area = Area_vec(i);
    end
    
    if length(A_r_vec) >1
    A_r = A_r_vec(i);
    end
    
    if length(iDist_vec) >1
    iDist = iDist_vec(i);
    end
    
    % reassing fishing pressure
    if length(mf_vec) >1
    SpParas.mf = mf_vec(i);
    end
    
    % Area vector (proportional, should sum to 1) 
    % (reserve restored, reserve, fished restored, fished)
    Area_all = [A_r, Area-A_r, A_r, 1-Area-A_r]; % [0.1, 0.1, 0.1, 0.7]; % [0.05, 0.15, 0.05, 0.75]; %  [0.2, 0, 0.2, 0.6]; % 

    
% Pre-reserv (intial conditions)-----       

% Fished initial conditions (no disturbance):
% SpParas.Sg = 1.01;
[NPre, ~, YBPre, ~] = DD_func('fished',fishaffect, SpParas,...
                                   repmat(100,SpParas.A_max*p,1),...
                                   tInit, Area_all, p,...
                                   zeros(tInit,p), zeros(tInit,p));
                               
% figure(100)
% hold on
% plot(1:tInit, sum(NPre,1))


% Post-reserve ----

% baseline equilibrium values (with reserves)
% NO disturbance
[NPostBase, NBPostBase, YPostBase, ~] = DD_func('reserve', fishaffect, SpParas,...
                             NPre(:,end,:), tPost, Area_all,...
                             p, zeros(tPost,p), zeros(tPost,p));                     

% Disturbance, no restoration
switch who_dist
    case 'juvs'
        JD1 = repmat(iDist,tPost,p);
        AD1 = zeros(tPost,p);
    case 'adults'
        JD1 = zeros(tPost,p);
        AD1 = repmat(iDist,tPost,p);
end
        
[NPostNR, NBPostNR, YPostNR, ~] = DD_func('reserve', fishaffect, SpParas,...
                             NPre(:,end,:), tPost, Area_all, p, JD1, AD1);
% 
% figure(100)
% hold on
% plot(1:tPost, sum(NPostBase,1)) 
% plot(1:tPost, sum(NPostNoRest,1))  

% With distrubance: 
    
    for x = 1:length(yrDist_Vec)
        for y = 1:length(rDist_Vec)    

 
	% year disturbance happens in (after res est)
    yrDist = yrDist_Vec(x);
	% length of the disturbance     
    rDist = rDist_Vec(y);
    
    % Disturbance affects survival:
    % set up empty matricies (resets them each time)
    LvarD = zeros(tPost,p);
    AdD = zeros(tPost,p);
    
    % set affected stage
    switch who_dist
        % Larvae affected
        case 'juvs'
            LvarD = Func_pop_rest(rest_pop, LvarD, iDist, yrDist, rDist);
        
        % Adults affected
        case 'adults'
            AdD = Func_pop_rest(rest_pop, AdD, iDist, yrDist, rDist);
    end
    
    
    
    % After MPA with disturbance
[NPostD, NBPostD, YBPostD, Rs] = DD_func('reserve', fishaffect, SpParas, ...
                             NPre(:,end,:), tPost, Area_all,...
                             p, LvarD, AdD);


% population & yeild biomass at equilibrium
t_eq = 100;
NM(i) = sum(NBPostD(:,t_eq),1);
N1(i) = sum(NBPostD(1:SpParas.A_max,t_eq),1);
N2(i) = sum(NBPostD(SpParas.A_max+1:SpParas.A_max*2,t_eq),1);
N3(i) = sum(NBPostD(2*SpParas.A_max+1:SpParas.A_max*3,t_eq),1);
N4(i) = sum(NBPostD(3*SpParas.A_max+1:SpParas.A_max*4,t_eq),1);

YM(i) = sum(YBPostD(:,t_eq),1);
Y1(i) = sum(YBPostD(1:SpParas.A_max,t_eq),1);
Y2(i) = sum(YBPostD(SpParas.A_max+1:SpParas.A_max*2,t_eq),1);
Y3(i) = sum(YBPostD(2*SpParas.A_max+1:SpParas.A_max*3,t_eq),1);
Y4(i) = sum(YBPostD(3*SpParas.A_max+1:SpParas.A_max*4,t_eq),1);  
                
% no disturbance no rest
NBM(i) = sum(NBPostBase(:,t_eq),1);
NB1(i) = sum(NBPostBase(1:SpParas.A_max,t_eq),1);
NB2(i) = sum(NBPostBase(SpParas.A_max+1:SpParas.A_max*2,t_eq),1);
NB3(i) = sum(NBPostBase(2*SpParas.A_max+1:SpParas.A_max*3,t_eq),1);
NB4(i) = sum(NBPostBase(3*SpParas.A_max+1:SpParas.A_max*4,t_eq),1);

YBM(i) = sum(YPostBase(:,t_eq),1);
YB1(i) = sum(YPostBase(1:SpParas.A_max,t_eq),1);
YB2(i) = sum(YPostBase(SpParas.A_max+1:SpParas.A_max*2,t_eq),1);
YB3(i) = sum(YPostBase(2*SpParas.A_max+1:SpParas.A_max*3,t_eq),1);
YB4(i) = sum(YPostBase(3*SpParas.A_max+1:SpParas.A_max*4,t_eq),1);

% diturbance no rest
NNRM(i) = sum(NBPostNR(:,t_eq),1);
NNR1(i) = sum(NBPostNR(1:SpParas.A_max,t_eq),1);
NNR2(i) = sum(NBPostNR(SpParas.A_max+1:SpParas.A_max*2,t_eq),1);
NNR3(i) = sum(NBPostNR(2*SpParas.A_max+1:SpParas.A_max*3,t_eq),1);
NNR4(i) = sum(NBPostNR(3*SpParas.A_max+1:SpParas.A_max*4,t_eq),1);

YNRM(i) = sum(YPostNR(:,t_eq),1);
YNR1(i) = sum(YPostNR(1:SpParas.A_max,t_eq),1);
YNR2(i) = sum(YPostNR(SpParas.A_max+1:SpParas.A_max*2,t_eq),1);
YNR3(i) = sum(YPostNR(2*SpParas.A_max+1:SpParas.A_max*3,t_eq),1);
YNR4(i) = sum(YPostNR(3*SpParas.A_max+1:SpParas.A_max*4,t_eq),1); 


if single_run == 'y'
% Plotting individual scenarios:
    tx = 1:tPost; 
    switch rest_pop
        case 'rest_res'
           lincol = '-g';
        case 'rest_fish'
            lincol = '--b';
    end 

    if rDist > 100
     lincol = '--r'; 
    end
    
    % biomass
    figure(331)
    hold on
    Func_OverTimePlot(NBPostD, tx, lincol, SpParas)
    Func_OverTimePlot(NBPostBase, tx, 'k', SpParas)
    ylabel('Biomass')
    
%     % abundance
% 	figure(332)
%     hold on
%     Func_OverTimePlot(NPostD, tx, lincol, SpParas)
%     Func_OverTimePlot(NPostBase, tx, 'k', SpParas)
%     ylabel('Abundance')
%      
    % yield biomass
	figure(333)
    hold on
    Func_OverTimePlot(YBPostD, tx, lincol, SpParas)
    Func_OverTimePlot(YPostBase, tx, 'k', SpParas)
    ylabel('Yield Bio')
      
    
end


        end 
     end
end


%% Plotting: Equilb values vs varied variable 

figure(1) 

switch rest_pop
    case 'rest_res'
       lincol = 'g';
    case 'rest_fish'
        lincol = 'r';
    case 'rest_both'
        lincol = 'm';
end

% Subplots
% need to chose the appropriate plot and associate labels

subplot(1,2,1)
% plot(Var_vec, NM./NNRM, '--', 'Color', lincol)
% plot(Var_vec./SpParas.m, NM./NNRM, '-', 'Color', lincol)
semilogx(Sg, NM./NNRM, '-', 'Color', lincol)
xline(0.1, "--k")
% grid on
% ylim([1, 1.8])
% xlim([10^-5, 10^-2])
ylabel('Metapopulation Biomass (as prop of dist & no rest state)')
% xlabel('Proportion of total area restored (area protected = 0.2)')
% xlabel('Proportion of area protected (area restored = 0.1)')
% xlabel('Fishing pressure (F/M) (area protect = 0.2, area restored = 0.2)')
% xlabel('Magnitude of disturbance (gamma) (area protect = 0.2, area restored = 0.2)')
% xlabel('Proportion of area protected (area restored = area pro)')
xlabel('Strength of fish-kelp itneraction (gamma_s)')
hold on

subplot(1,2,2)
% plot(Var_vec, YM./YNRM, '--', 'Color', lincol)
% plot(Var_vec./SpParas.m, YM./YNRM, '-', 'Color', lincol)
semilogx(Sg, YM./YNRM, '-', 'Color', lincol)
xline(0.1, "--k")
% grid on
% ylim([1, 1.8])
% xlim([10^-5, 10^-2])
ylabel('Yield Biomass (as prop of dist & no rest state)')
hold on


%% Plotting fish-to-habitat interaction

clear 

Abio = 1:10^3;
G0 = 0.25;
Gs = [0,logspace(-6,0,7)];


figure
hold on
for i = 1:length(Gs)
    G = G0 - (Gs(i).*Abio)./(1+(Gs(i).*Abio)./G0);
    
    plot(Abio,G)
    xlabel('Population biomass (B_i)')
    ylabel('\gamma_i (B_i)')
end

% add biomass range from Caselle et al 2018 (Fig.3)
xline(0.7*10^2, "--k")
xline(3.2*10^2, "--k")

