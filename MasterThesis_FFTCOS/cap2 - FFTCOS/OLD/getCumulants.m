function resu = getCumulants(model,t,r,q,param)
% Provides a library of cumulants for various models. Those quantities are
% necessary to truncate the integration range [a,b] in order to compute the 
% inverse Fourier integral using the cosine series expansion.
% 
% Literature: see the paper by F. Fang and C. W. Oosterlee [2008],
% «A Novel Pricing Method for European Options Based on 
% Fourier-Cosine Series Expansions».
%
% INPUTS
% x
% y
%
% OUTPUT
% resu 
%

ME1 = MException('InconsistentInput:InvalideNumberOfArguments',...
                 'Invalid number of input arguments');

if strcmp(model,'BlackScholes')
    % List of parameters
    sigma = param(1);
    % Give parameters to characteristic function
    resu = getCumulants_bs(t,r,q,sigma);

elseif strcmp(model,'Merton')
    % List of parameters
    sigma = param(1);
    muJ = param(2);
    sigmaJ = param(3);
    lambda = param(4);
    % Give parameters to characteristic function
    resu = getCumulants_merton(t,r,q,sigma,muJ,sigmaJ,lambda);

elseif strcmp(model,'Kou')
    % List of parameters
    sigma = param(1);
    p1 = param(2);
    eta1 = param(3);
    eta2 = param(4);
    lambda = param(5);
    % Give parameters to characteristic function
    resu = getCumulants_kou(t,r,q,sigma,p1,eta1,eta2,lambda);

elseif strcmp(model,'VG')
    % List of parameters
    sigmaVG = param(1);
    beta = param(2);
    theta = param(3);
    S0 = param(4);
    % Give parameters to characteristic function
    resu = getCumulants_vg(t,r,q,sigmaVG,beta,theta,S0);

elseif strcmp(model,'Heston')
    % List of parameters
    kappa = param(1);
    theta = param(2);
    eta = param(3);
    rho = param(4);
    v0 = param(5);
    % Give parameters to characteristic function
    resu = getCumulants_heston(t,r,q,kappa,theta,eta,rho,v0);

elseif strcmp(model,'Bates')
    % List of parameters
    kappa = param(1);
    theta = param(2);
    eta = param(3);
    rho = param(4);
    v0 = param(5);
    muJ = param(6);
    sigmaJ = param(7);
    lambda = param(8);
    % Give parameters to characteristic function
    resu = getCumulants_bates(t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda);

% elseif strcmp(model,'WMSV') TODO
    % List of parameters
    % TODO
    % Give parameters to characteristic function
    % TODO

else
    warning('Method not found')
end

end

%-------------------------------------------------------------------------%
% Here we effectively compute the cumulants of various model
%  
% Below we can add more functions
%-------------------------------------------------------------------------%

% Black Scholes
function c = getCumulants_bs(t,r,q,sigma)
    c = zeros(1,4);
    % First cumulant
    c(1) = (r - q - 0.5*sigma^2) * t; 
    % Second cumulant
    c(2) = sigma^2 * t;
    % High order cumulants
    c(3) = 0.0; c(4) = 0.0;
end

% Merton Jump diffusion
function c = getCumulants_merton(t,r,q,sigma,muJ,sigmaJ,lambda)
    c = zeros(1,4);
    % Term for E[exp(J) - 1]
    EX1 = exp(muJ + 0.5 * sigmaJ * sigmaJ) - 1.0;
    % Drift correction term
    mu = (r - q) - 0.5 * sigma^2 - lambda * EX1;
    % First cumulant
    c(1) = t * (mu + lambda*muJ);
    % Second cumulant
    c(2) = t * (sigma^2 + lambda*muJ^2 + sigmaJ^2*lambda);
    % High order cumulants
    c(3) = 0.0;
    c(4) = t * lambda * (muJ^4 + 6*sigmaJ^2*muJ^2 + 3*sigmaJ^4*lambda);
end

% Kou 
function c = getCumulants_kou(t,r,q,sigma,p1,eta1,eta2,lambda)
    c = zeros(1,4);
    p2 = 1 - p1;
    % Term for E[J]
    EX1 = (p1 * eta1)/(eta1 - 1.0) + (p2 * eta2)/(eta2 + 1.0);
    % Drift correction term
    mu = (r - q) - 0.5 * sigma^2 - lambda * EX1; 
    % First cumulant
    c(1) = t * (mu + lambda*p1*(1.0/eta1) - lambda*p2*(1.0/eta2));
    % Second cumulant
    c(2) = t * (sigma^2 + 2*lambda*p1*(1.0/eta1^2)...
               + 2*lambda*p2*(1.0/eta2^2));
    % High order cumulants
    c(3) = 0.0;
    c(4) = 24 * t * lambda * (p1/(eta1^4) + p2/(eta2^4));
end

% Variance Gamma [VG]
function c = getCumulants_vg(t,r,q,sigmaVG,beta,theta)
    c = zeros(1,4);
    % Correction for martingality
    omega = (1/beta) * log(1 - theta*beta - 0.5*sigmaVG^2*beta);
    % First cumulant
    c(1) = (r - q - omega + theta) * t;
    % Second cumulant
    c(2) = (sigmaVG^2 + beta*theta^2) * t;
    % High order cumulants
    c(3) = 0.0; 
    c(4) = 3 * (sigmaVG^4*beta + 2*theta^4*beta^3 ...
               + 4*sigmaVG^2*theta^2*beta^2) * t;
end

% Heston Stochastic Volatility 
function c = getCumulants_heston(t,r,q,kappa,theta,eta,rho,v0)
    c = zeros(1,4); 
    % First cumulant
    c(1) = (r-q)*t + (theta - v0)/(2*kappa) * (1 - exp(-kappa*t)) - 0.5*theta*t;
    % Second cumulant
    c(2) = 1/8*kappa^3 *... 
        ( kappa*eta*t*exp(-kappa*t) * (v0 - theta) * (8*kappa*rho - 4*eta)...
          + 8*kappa*rho*eta * (1 - exp(-kappa*t)) * (2*theta - v0)...
          + 2*kappa*theta*t * (-4*kappa*rho*eta + eta^2 + 4*kappa^2) ...
          + eta^2 * ( (theta - 2*v0) * exp(-2*kappa*t) + theta * (6*exp(-kappa*t) - 7) + 2*v0 ) ...
          + 8*kappa^2 * (v0 - theta) * (1 - exp(-kappa*t)) ) ;
    % High order cumulants
    c(3) = 0.0; c(4) = 0.0;
end

% Bates Jump diffusion stochastic volatility
function c = getCumulants_bates(t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda)
    c = zeros(1,4);
    %-----Heston cumulants-----%
    c1_hes = (r-q)*t + (theta - v0)/(2*kappa) * (1 - exp(-kappa*t)) - 0.5*theta*t;
    c2_hes = 1/8*kappa^3 *... 
        ( kappa*eta*t*exp(-kappa*t) * (v0 - theta) * (8*kappa*rho - 4*eta)...
          + 8*kappa*rho*eta * (1 - exp(-kappa*t)) * (2*theta - v0)...
          + 2*kappa*theta*t * (-4*kappa*rho*eta + eta^2 + 4*kappa^2) ...
          + eta^2 * ( (theta - 2*v0) * exp(-2*kappa*t) + theta * (6*exp(-kappa*t) - 7) + 2*v0 ) ...
          + 8*kappa^2 * (v0 - theta) * (1 - exp(-kappa*t)) ) ;
    %-----Jump cumulants-----%
    % Term for E[exp(J) - 1]
    EX1 = exp(muJ + 0.5 * sigmaJ * sigmaJ) - 1.0;
    % Correction term
    omega = lambda * EX1;
    % First cumulant
    c1_jump = t * (-1.0*omega + lambda*muJ);
    % Second cumulant
    c2_jump = t * (lambda*muJ^2 + sigmaJ^2*lambda);
    %-----Bates cumulants-----%
    % First cumulants
    c(1) = c1_hes + c1_jump;
    % Second cumulant
    c(2) = c2_hes + c2_jump;
    % High order cumulants
    c(3) = 0.0; c(4) = 0.0;
end
