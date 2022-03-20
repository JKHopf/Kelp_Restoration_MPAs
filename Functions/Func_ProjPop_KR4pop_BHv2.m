function [N, Nbio, Ybio, Rs, S, R] = Func_ProjPop_KR4pop_BHv2(scenario, fishaffect, SpParas,...
                                      Ninit, Rtime, Area_vec, p,...
                                      Rs, As)
                                  
% Projecting the popualtion over time
% Can be used for both fully fished or spatially managed pops
%
% Jess Hopf
% May 2021
%
% Inputs:
%   - scenario = 'fished' or 'reserves'
%   - SpParas = single row table of parameters for the species
%   - Ninit = initial population size (matrix size = age-max x nsim)
%   - Rtime = the run time of the projection
%   - Area_vec = area in populations
%   - nsim = number of simulations
%   - p = number of populations
%   - Rs = matrix of recrtuiment disturbance multiplied mortality over time
%   - As = matrix of adult disturbance additional mortality over time
%

% ------------------------------------------------------------------------
% Pre-allocate matrices:
    N = NaN(SpParas.A_max * p, Rtime);
    Nbio = N;
    
% ------------------------------------------------------------------------    
% Set up demography

% Age, length, weight relationships:
    % Per capita length (cm) at start of age year (L)
    Lengths = Func_Length((1:SpParas.A_max)',...
                             SpParas.L_inf,SpParas.K,SpParas.A0);

    % Per capita weight at the length at the start of age year (W)
    Weights = Func_Weights(Lengths,SpParas.y,SpParas.z);
    
% Fecundity:
    % Per capita eggs produced           
        % Density indep fecundity (per capita) at length at the start 
        % of age year
        % lengths need to be in mm
        % divide by 2, assuming 50:50 sex ratio
        Fun = Func_fecunds((1:SpParas.A_max)',Lengths*10,...
                             SpParas.c, SpParas.d);

    % Set ages which reproduce
     Fun(1:(SpParas.A_mat-1))=0;   
     
% ------------------------------------------------------------------------
% run over time (inc time-dependent parameters)

% Intial conditions:
N(:,1,:) = Ninit;

% first year disturbance happens
Yr_distR = find(sum(Rs, 2) < p, 1);
if isempty(Yr_distR)
    Yr_distR = 0;
end
Yr_distA = find(sum(As, 2) < p, 1);
if isempty(Yr_distA)
    Yr_distA = 0;
end

% Fishing mortality:
    switch scenario
        case 'fished'
          mf = SpParas.mf;  
        case 'reserve'
          mf = SpParas.mf/(1-sum(Area_vec(1:2)));  
        case 'unfished'
          mf = 0;
    end
    
% Survival probabilities:
    % survive natural mortality 
    SrU = exp(-SpParas.m);   
    % survive natural & fishing mort 
    % (inc reallocated effort if reserve scenario)
    SrF = exp(-(SpParas.m+mf));
    
% Proportion caught
    Df = (mf./(SpParas.m+mf)).*(1-exp(-SpParas.m-mf));
   
for t = 1:Rtime-1 
% Demographic matrices:
    % Reefs that are fished 
    RF = diag([repmat(SrU,1,(SpParas.Ac-1)),...
          repmat(SrF,1,SpParas.A_max-SpParas.Ac)],-1);
    RF(1,:) = Fun;

    % Reefs that are not fished
    RR = diag(repmat(SrU,1,(SpParas.A_max-1)),-1);
    RR(1,:) = Fun;   
    
    % Combine pops
    switch scenario
        case 'fished' 
            B = blkdiag(RF,RF,RF,RF);            
        case 'reserve'
            B = blkdiag(RR,RR,RF,RF);   
        case 'unfished'
            B = blkdiag(RR,RR,RR,RR);   
    end
    
    % Transition matrix
    M = B; 
    R = NaN(1,p);
    
    % density of incoming larvae (same for both pop under well-mixed
    % assumption)
    S = sum(repmat(Fun,p,1).*N(:,t,:));
    
    % for individual populations
    for i = 1:p    
    % M matrix index values for start/end of population
    Mmin = i*SpParas.A_max-SpParas.A_max+1;
    Mmax = i*SpParas.A_max;
    
    % If fish dependent survival rates (As & Rs)
    % adult (2+ years) biomass density in population
    % make zero if Area_vec(i) = 0
    Abio = sum(N(i*SpParas.A_max-SpParas.A_max+2:i*SpParas.A_max,t,:)...
        .* Weights(2:end))./Area_vec(i);
    
    if Area_vec(i) == 0
        Abio = 0;
    end
    

    % if fish have an effect on kelp, then include functional response 
    % this only begins to apply the year after the disturbance is done
    % two diff functional forms: top is ricker-esque, bottom is BH-esque
    if strcmp(fishaffect, 'yes')
        if Rs(t,i)>0 && t>Yr_distR
%             Rs(t,i) = SpParas.Sg^(-Abio+log(Rs(t,i))/log(SpParas.Sg)); % Rs(t,i) + Rg.*Abio;
            Rs(t,i) = Rs(t,i) - (SpParas.Sg.*Abio)./(1+(SpParas.Sg.*Abio)./Rs(t,i));
        end
        if As(t,i)>0 && t>Yr_distA
%             As(t,i) = SpParas.Sg^(-Abio+log(As(t,i))/log(SpParas.Sg)); % As(t,i) + Ag.*Abio;
            As(t,i) = As(t,i) - (SpParas.Sg.*Abio)./(1+(SpParas.Sg.*Abio)./As(t,i));
        end
    end
     
    % DD p.c. survival - Beverton-Holt function (recruits affect each other)
    % includes different disturbance effects in each pop
    % changes alpha in BH 
%     R(i) = (SpParas.BHa.*(1-Rs(t,i)))./(1+(SpParas.BHa.*(1-Rs(t,i)).*S)./SpParas.BHb);
    % changes beta in BH 
%     R(i) = SpParas.BHa./(1+(SpParas.BHa.*S)./(SpParas.BHb.*(1-Rs(t,i))));
    % affects overall survival
    R(i) = (SpParas.BHa./(1+(SpParas.BHa.*S)./(SpParas.BHb))).*(1-Rs(t,i));
   
    % Multiply recruits into projection matrix
    % i.e. fecundity * area * p.c. survival                                    
    M(Mmin,:,:,:) = repmat(M(Mmin,Mmin:Mmax,:,:)...
                           .*Area_vec(i).*R(i),1,p);
            
    % multiply in disturbance effects on adults
    M(Mmin:Mmax,Mmin:Mmax)= M(Mmin:Mmax,Mmin:Mmax)+...
                            diag(diag(M(Mmin:Mmax,Mmin:Mmax),-1).*(1-As(t,i))-...
                            diag(M(Mmin:Mmax,Mmin:Mmax),-1),-1);

    end
    
    N(:,t+1,:,:) = pagemtimes(M,N(:,t,:,:));      

end

    Nbio = N.*repmat(Weights,p,Rtime);
    

   switch scenario
        case 'fished' 
            Y = repmat([zeros(SpParas.Ac-1,size(N,2)); % fished ages only
                    repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],4,1).*N; 
            Ybio = repmat([zeros(SpParas.Ac-1,size(N,2)); % fished ages only
                    repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],4,1).*Nbio; 
        case 'reserve'
            Y = [zeros(SpParas.A_max*2,size(N,2)); % pops 1 & 2, fished ages only
                 repmat([zeros(SpParas.Ac-1,size(N,2)); % pops 3 & 4, fished ages only
                    repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],2,1)].*N; 
            Ybio = [zeros(SpParas.A_max*2,size(N,2)); % pops 1 & 2, fished ages only
                    repmat([zeros(SpParas.Ac-1,size(N,2)); % pops 3 & 4, fished ages only
                        repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],2,1)].*Nbio;
        case 'unfished'
            Y = 0.*N;
            Ybio = 0.*Nbio;
    end



