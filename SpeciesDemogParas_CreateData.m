% SpeciesDemogParas_CreateData.m
%
% Jess Hopf
% June 2020
%
% This file creates .m data files holding the parameter estimates for
% target nearshore species in California.
%
% Data manually entered from:
%   - Kaplan et al. (2019) Ecol App. (Table 2)
%
% Parameter descriptions:
%   sp = species name
%   von Bertalanffy growth paras 
%       Linf = Asymptotic max length (cm)
%       K = growth para
%       A0 = age at length 0
%   y, z = weight-length paras 
%   A_mat = age at maturity
%   Ac = entry to fishery age (age of first capture)
%   A_max = max age
%   L_mat = length at maturity
%   Lc = length of first capture
%   m = natural mortality rate
%   mf = fishing mortality rate (estimated)
%
%  c, d = fecundity at length length para*

% ------------------------------------------------------------------------

clear

sp = ["California sheephead";"Kelp bass"];
  
A_max = [53;33];

% L_mat = [24.00;22.30];
     
% Lc = [30.48;30.48];

m = [0.3;0.18];  % SH: mean from Caselle et al (2011) 

mf = [0.25;0.12];
 
L_inf = [83.86;69.80]; % SH: Alonzo et al (2004)
     
K = [0.068;0.06]; % SH: Alonzo et al (2004)
 
A0 = [-1.47;-3.5]; % SH: solved VBF with L_0 = 8cm, Alonzo et al (2004)
  
A_mat = [4;4];

Ac = [6;6];  % SH: Alonzo et al (2004), length @ 1st capture = 13 inches/33cm

y = [0.27;0.0273].*0.0001; % SH: Alonzo et al (2004)

z = [2.86;3.27]; % SH: Alonzo et al (2004)

c = [7.04; exp(-5.57)]; % KB: Oda et al 1993, SH: Easter et al 2020
    
d = [2.95; 2.93]; % KB: Oda et al 1993, SH: Easter et al 2020
 
SppParas = table(sp,L_inf,K,A0,y,z,A_mat,Ac,A_max,m,mf,c,d);

save SppParameterVals.mat SppParas
 