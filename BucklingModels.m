% Copyright (C) 2019 Petros Lazaridis


% Author: Petros Lazaridis
% Created: 2019-11-23

function [sigmaT, sigmaDM, f_sc] = BucklingModels (eps, L, D, eps_y, eps_sh,...
                                                  eps_u, fy, fu)
% All units must be in mm and MPa.     
% fy   : Yield strength (in MPa) 
% fu   : Ultimate strength (in MPa)
% eps_y : Yield strain 
% eps_sh : Strain at once strain hardening 
% eps_u : Ultimate strain 
% L    : Unsupported length (mm)
% D    : Diameter of circular section (mm)  
  P = 1;
  Es = fy/eps_y;
  rb = L./D*sqrt(fy/100);
  
  sigmaT = interp1([0  eps_y eps_sh eps_u], [0  fy fy fu], eps, 'linear');
  condMainEq = eps > eps_sh & eps < eps_u;
  sigmaT(condMainEq) = fu + (fy-fu)*((eps_u-eps(condMainEq))/(eps_u - eps_sh)).^P;
  
  eps_i0 = (55 - 2.3*rb) * eps_y;
  eps_i = eps_i0;
  alpha2 = (1.1 - 0.016*rb);

  %-------------- Refined Dhakal & Maekawa (RDM) Buckling Model ----------------
  % Yildir Akkaya Sehran Guner and Frank J.Vecchio, Constitutive Model for Ine-
  % lastic Behaviour of Reinforcing Bars, ACI Structural Journal, p. 195-204,
  % March 2019
  %-----------------------------------------------------------------------------
  sigma = sigmaT;
  fracLD_min = 5;
  rb_min = fracLD_min * sqrt(fy/100);
  eps_imax = (55 - 2.3*rb_min) * eps_y;
  cond1fromEq3 = eps_i0 < eps_u & eps_u < eps_imax;
  eps_i(cond1fromEq3) = eps_i0(cond1fromEq3).*eps_u./eps_imax;
  eps_i(eps_i < 7*eps_y) = 7*eps_y;
  
  f_it = interp1(eps, sigma, eps_i, 'linear');
  alpha1 = .8 + 1.8*(fu/fy)*D/L;
  alpha = .75*alpha1.*alpha2;
  alpha(eps_i > eps_sh) = alpha1(eps_i > eps_sh).*alpha2(eps_i > eps_sh);
  eq5LastCond = eps_u < eps_imax & eps_i == 7*eps_y;
  alpha(eq5LastCond) = .75*alpha2(eq5LastCond).*(f_it(eq5LastCond)/fy);
  
  f_i = alpha*fy; % Intermediate point stress
  
  f_i(f_i < .2*fy) = .2*fy;
  f_i(f_i > f_it) = f_it;
  f_sc = sigma;
  condEq7 = eps > eps_y & eps < eps_i;
  f_sc(condEq7) = f_sc(condEq7).*(1- (1-f_i/f_it)*(eps(condEq7)-eps_y)/(eps_i-eps_y));
  
  eps_ii = eps_i + .25*f_i/.02/Es;
  cond1Eq8 = eps_i < eps & eps< eps_ii;
  f_sc(cond1Eq8) = f_i - .02*Es*(eps(cond1Eq8) - eps_i);
  cond2Eq8 = eps_ii < eps & eps_u >= eps ;
  f_sc(cond2Eq8) = .75*f_i - .01*Es*(eps(cond2Eq8)-eps_ii);
  f_sc(eps > eps_i0 & f_sc < .2*fy) = .2*fy;

  
  %f_i is equivalent to sigmaStar (intemediate stress)
  %f_it is equivalent to sigmaStar_l
  
  %------------------- Dhakal & Maekawa (DM) Buckling Model --------------------
  % Rajesh Prasad Dhakal and Koichi Maekawa, Modelling for Postyield Buckling of 
  % Reinforcement, Journal of Structural Engineering, vol.128, p. 1139-1147,
  % Sept 2002
  %-----------------------------------------------------------------------------
  sigma(condMainEq) = fu + (fy-fu)*((eps_u-eps(condMainEq))/(eps_u - eps_sh)).^1;
  alpha0_1 = 0.75 + (eps_u - eps_sh)/300/eps_y;
  alpha0_1(alpha0_1 > fu/1.5/fy) = fu/1.5/fy;
  alpha0_1(alpha0_1 < .75) = .75;
  alpha0_1(alpha0_1 > 1) = 1;
  rb_min = 8;
  eps_i0(eps_i0 < 7*eps_y) = 7*eps_y;
  sigmaStar_l = interp1(eps, sigma, eps_i0, 'linear');
  sigmaStar = alpha0_1*alpha2.*sigmaStar_l;
  sigmaStar(sigmaStar < 0.2*fy) = 0.2*fy;
  fracSigmaStars = sigmaStar./sigmaStar_l;
  cond1Eq1 = eps > eps_y & eps < eps_i0';
  sigmaDM = sigma;
  
  sigmaDM(cond1Eq1) = sigmaDM(cond1Eq1).*( 1- (1-fracSigmaStars)*(eps(cond1Eq1)- eps_y)./(eps_i0- eps_y)');
  sigmaDM(eps > eps_i0) = sigmaStar - 0.02*Es*(eps(eps > eps_i0) - eps_i0);
  sigmaDM(eps > eps_i0 & sigmaDM < .2*fy) = .2*fy;

end
