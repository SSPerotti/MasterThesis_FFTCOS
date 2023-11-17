function resu = getCumulants(model,t,varargin)
% GETCUMULANTS(model,t,varargin);
%
% Provides a library of cumulants for various models. Those quantities are
% necessary to truncate the integration range [a,b] in order to compute the 
% inverse Fourier integral using the cosine series expansion.
% 
% Literature: see the paper by F. Fang and C. W. Oosterlee [2008],
%             «A Novel Pricing Method for European Options Based on 
%              Fourier-Cosine Series Expansions».
%
% Reference:  see Table 5.2 (page 223) in Financial Modelling. Theory, 
%             Implementation and Practice with MATLAB Source by Jorg Kienitz
%             and Daniel Wetterau.
%
% INPUTS
% model:       Chose model of the underlying dynamics. Model accepts as a 
%              parameter a string that represents the name of the model. 
%              The models implemented are: 'BlackScholes', 'Merton',
%              'Kou', 'VG', 'Heston', 'Bates', 'WMSV'.
% t:           Time-to-maturity
% 
% FUNCTIONS
% getCumulants_bs:      Compute the Black Scholes [1973] model cumulants.
%                       The only parameter is: 
%                           - sigma = volatility.
%
% getCumulants_merton:  Compute the Merton Jump diffusion [1976] model cumulants
%                       The parameters are
%                           - sigma = volatility,
%                           - muJ = jump size average,
%                           - sigmaJ = jump standard deviation,
%                           - lambda = intensity of the Poisson process.
%  
% getCumulants_kou:     Compute the Kou [2002] model cumulants. 
%                       The parameters are:
%                           - sigma = volatility,
%                           - p1 = probability of positive jump,
%                           - eta1 = positive jump size,
%                           - eta2 = negative jump size,
%                           - lambda = intensity of the Poisson process.
%
% getCumulants_vg:      Compute the Variance Gamma model cumulants. 
%                       The parameters are:
%                           - sigmaVG = TODO
%                           - beta = TODO
%                           - theta = TODO
%
% getCumulants_heston:  Compute the Heston [1993] model cumulants.
%                       The parameters are:
%                           - kappa = speed of mean reversion, 
%                           - theta = long run average, 
%                           - eta = vol-of-vol, 
%                           - rho = correlation, 
%                           - v0 = initial variance value. 
%
% getCumulants_bates:   Compute the Bates [1996] model cumulants. 
%                       The parameters for the Heston part are: 
%                           - kappa = speed of mean reversion, 
%                           - theta = long run average, 
%                           - eta = vol-of-vol, 
%                           - rho = correlation, 
%                           - v0 = initial variance value. 
%                       Instead for the jump part we have: 
%                           - muJ = jump size average, 
%                           - sigmaJ = jump standard deviation, 
%                           - lambda = intensity of the Poisson process.
%
% getCumulants_wmsv:    Compute the Wishart Heston Volatility model cumulants via
%                       Finite Difference.
%                       The model parameters are:
%                         - M11, M12, M21, M22 = elements of (2x2) negative
%                                                semidefinite matrix in order to
%                                                ensure mean-reversion
%                         - R11, R12, R21, R22 = elements of (2x2) matrix that
%                                                represents correlation between
%                                                matrix brownian motions
%                         - Q11, Q12, Q21, Q22 = elements of (2x2) matrix
%                         - Sigma11, Sigma12, 
%                           Sigma21, Sigma22   = elements of (2x2) positive
%                                                semidefinite matrix of initial
%                                                state of variance process
%                         - beta = Gindikin parameter has the same role of
%                                  Feller's conditions in the square-root process.
%
% OUTPUT
% resu:        Return a vector that contain the cumulants c1, c2 and c4. 
%
% ------------------------------------------------------------------------

ME1 = MException('InconsistentInput:InvalideNumberOfArguments',...
                 'Invalid number of input arguments');

if strcmp(model,'BlackScholes')
    if nargin == 5
        funObj = @getCumulants_bs;
    else
        throw(ME1);
    end
elseif strcmp(model,'Merton')
    if nargin == 8
        funObj = @getCumulants_merton;
    else
        throw(ME1);
    end
elseif strcmp(model,'Kou')
    if nargin == 9
        funObj = @getCumulants_kou;
    else
        throw(ME1)
    end
elseif strcmp(model,'VG')
    if nargin == 7
        funObj = @getCumulants_vg;
    else
        throw(ME1)
    end
elseif strcmp(model,'Heston')
    if nargin == 9
        funObj = @getCumulants_heston;
    else
        throw(ME1);
    end
elseif strcmp(model,'Bates')
    if nargin == 12
        funObj = @getCumulants_bates;
    else
        throw(ME1)
    end
elseif strcmp(model,'WMSV') 
    if nargin == 21
        funObj = @getCumulants_wmsv;
    else
        throw(ME1)
    end
end

resu = feval(funObj,t,varargin{:});

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
    % Recover high order cumulants via finite difference
    i = complex(0,1); 
    h = 1e-2;
    psi = @(u) log(getCharacteristicFunction('Heston',-i*u,t,r,q, ...
                   kappa,theta,eta,rho,v0));
    fh = psi(0);
    f1h = psi(1.0*h);  
    f_1h = psi(-1.0*h);
    f2h = psi(2.0*h);  
    f_2h = psi(-2.0*h);
    % High order cumulants
    c(3) = 0.0;
    c(4) = (1.0*f2h - 4.0*f1h + 6.0*fh - 4.0*f_1h + 1.0*f_2h) / (1.0 * h^4);
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
    % Recover high order cumulants via finite difference
    i = complex(0,1); 
    h = 1e-2;
    psi = @(u) log(getCharacteristicFunction('Bates',-i*u,t,r,q, ...
                   kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda));
    fh = psi(0);
    f1h = psi(1.0*h);  
    f_1h = psi(-1.0*h);
    f2h = psi(2.0*h);  
    f_2h = psi(-2.0*h);
    % High order cumulants
    c(3) = 0.0;
    c(4) = (1.0*f2h - 4.0*f1h + 6.0*fh - 4.0*f_1h + 1.0*f_2h) / (1.0 * h^4);
end


function c = getCumulants_wmsv(t,r,q, ...
                               m11, m12, m21, m22, ...
                               r11, r12, r21, r22, ...
                               q11, q12, q21, q22, ...
                               Sig11, Sig12, Sig21, Sig22, ...
                               beta)
    i = complex(0,1); 
    h = 1e-2;
    psi = @(u) log(getCharacteristicFunction('WMSV',-i*u,t,r,q, ...
                                              m11, m12, m21, m22, ...
                                              r11, r12, r21, r22, ...
                                              q11, q12, q21, q22, ...
                                              Sig11, Sig12, Sig21, Sig22, ...
                                              beta));
    fh = psi(0);
    f1h = psi(1.0*h);  
    f_1h = psi(-1.0*h);
    f2h = psi(2.0*h);  
    f_2h = psi(-2.0*h);
    % First cumulant
    c(1) = (1.0*f1h - 1.0*f_1h) / (2.0 * h);
    % Second cumulant
    c(2) = (1.0*f1h - 2.0*fh + 1.0*f_1h) / h^2;
    % High order cumulants
    c(3) = 0.0;
    c(4) = (1.0*f2h - 4.0*f1h + 6.0*fh - 4.0*f_1h + 1.0*f_2h) / (1.0 * h^4);
end
