function resu = getCharacteristicFunction(model,u,t,r,q,varargin)
% GETCHARACTERISTICFUNCTION(model,u,t,r,q,varargin) 
%
% It provides a library of some characteristic functions for various model.
%
% INPUTS
% model:       Chose model of the underlying dynamics. Model accepts as a 
%              parameter a string that represents the name of the model. 
%              The models implemented are: 'BlackScholes', 'Merton',
%              'Kou', 'VG', 'Heston', 'Bates', 'WMSV'.
% u:           Valuation point of the characteristic function
% t:           Time-to-maturity
% r:           Risk free interest rate
% q:           Dividend yield
% 
% FUNCTIONS
% chf_bs:      Compute the Black Scholes [1973] model characteristic
%              function. The only parameter is: 
%                 - sigma = volatility.
%
% chf_merton:  Compute the Merton Jump diffusion [1976] model 
%              characteristic function. The parameters are
%                 - sigma = volatility,
%                 - muJ = jump size average,
%                 - sigmaJ = jump standard deviation,
%                 - lambda = intensity of the Poisson process.
%  
% chf_kou:     Compute the Kou [2002] model characteristic function. The 
%              parameters are:
%                 - sigma = volatility,
%                 - p1 = probability of positive jump,
%                 - eta1 = positive jump size,
%                 - eta2 = negative jump size,
%                 - lambda = intensity of the Poisson process.
%
% chf_vg:      Compute the Variance Gamma model characteristic function. 
%              The parameters are:
%                 - sigmaVG = TODO
%                 - beta = TODO
%                 - theta = TODO
%
% chf_heston:  Compute the Heston [1993] model characteristic function. 
%              The parameters are:
%                 - kappa = speed of mean reversion, 
%                 - theta = long run average, 
%                 - eta = vol-of-vol, 
%                 - rho = correlation, 
%                 - v0 = initial variance value. 
%
% chf_bates:   Compute the Bates [1996] model characteristic function. 
%              The parameters for the Heston part are: 
%                 - kappa = speed of mean reversion, 
%                 - theta = long run average, 
%                 - eta = vol-of-vol, 
%                 - rho = correlation, 
%                 - v0 = initial variance value. 
%              Instead for the jump part we have: 
%                 - muJ = jump size average, 
%                 - sigmaJ = jump standard deviation, 
%                 - lambda = intensity of the Poisson process.
%
% chf_wmsv:    Compute the Wishart Heston Volatility model characteristic
%              function (see Da Fonseca, Grasselli and Tebaldi [2008]). 
%              The model parameters are:
%                 - M11, M12, M21, M22 = elements of (2x2) negative
%                                        semidefinite matrix in order to
%                                        ensure mean-reversion
%                 - R11, R12, R21, R22 = elements of (2x2) matrix that
%                                        represents correlation between
%                                        matrix brownian motions
%                 - Q11, Q12, Q21, Q22 = elements of (2x2) matrix
%                 - Sigma11, Sigma12, 
%                   Sigma21, Sigma22   = elements of (2x2) positive
%                                        semidefinite matrix of initial
%                                        state of variance process
%                 - beta = Gindikin parameter has the same role of
%                          Feller's conditions in the square-root process.
%
% OUTPUT
% resu:        Return the characteristic function in the form Psi(u,t)
%              for a generic model. The initial condition i*u*X(t_0) 
%              where X = log(S_t) is not included. 
%
% ------------------------------------------------------------------------

ME1 = MException('InconsistentInput:InvalideNumberOfArguments',...
                 'Invalid number of input arguments');

if strcmp(model,'BlackScholes')
    if nargin == 6
        funObj = @chf_bs;
    else
        throw(ME1);
    end
elseif strcmp(model,'Merton')
    if nargin == 9
        funObj = @chf_merton;
    else
        throw(ME1);
    end
elseif strcmp(model,'Kou')
    if nargin == 10
        funObj = @chf_kou;
    else
        throw(ME1)
    end
elseif strcmp(model,'VG')
    if nargin == 8
        funObj = @chf_vg;
    else
        throw(ME1)
    end
elseif strcmp(model,'Heston')
    if nargin == 10
        funObj = @chf_heston;
    else
        throw(ME1)
    end
elseif strcmp(model,'Bates')
    if nargin == 13
        funObj = @chf_bates;
    else
        throw(ME1)
    end
elseif strcmp(model,'WMSV') 
    if nargin == 22
        funObj = @chf_wmsv;
    else
        throw(ME1)
    end
else
    warning('Method not found')
end


resu = feval(funObj,u,t,r,q,varargin{:});

end

%-------------------------------------------------------------------------%
% Here we effectively compute the characteristic function of various model
%  
% Below we can add more functions
%-------------------------------------------------------------------------%

% Black Scholes
function cf = chf_bs(u,t,r,q,sigma)
    i = complex(0,1);
    % Drift correction term
    mu = (r - q) - 0.5*sigma^2;
    % ChF for the Black Scholes model
    cf = exp(i*u*mu*t - 0.5*sigma^2*u.^2*t);  
end

% Merton Jump diffusion
function cf = chf_merton(u,t,r,q,sigma,muJ,sigmaJ,lambda)
    i = complex(0,1);
    % Term for E[exp(J) - 1]
    EX1 = exp(muJ + 0.5 * sigmaJ * sigmaJ) - 1.0;
    % Term for E[exp(iuJ) - 1]
    EX2 = exp(i*u*muJ - 0.5 * sigmaJ * sigmaJ * u.^2) - 1.0;
    % Drift correction term
    mu = (r - q) - 0.5 * sigma^2 - lambda * EX1;
    % ChF for the Merton model
    cf = exp(i*u*mu*t - 0.5*sigma*sigma*u.^2 + lambda*t*EX2);
end

% Kou 
function cf = chf_kou(u,t,r,q,sigma,p1,eta1,eta2,lambda)
    i = complex(0,1);
    p2 = 1 - p1;
    % Term for E[J]
    EX1 = (p1 * eta1)/(eta1 - 1.0) + (p2 * eta2)/(eta2 + 1.0);
    % Term for E[iuJ]
    EX2 = (p1 * eta1) ./ (eta1 - i*u) + (p2 * eta2) ./ (i*u + eta2);
    % Drift correction term
    mu = (r - q) - 0.5 * sigma^2 - lambda * EX1; 
    % Chf for the Kou model
    cf = exp(i*u*mu*t - 0.5*sigma*sigma*u.^2 + lambda*t*EX2);
end

% Variance Gamma [VG]
function cf = chf_vg(u,t,r,q,sigmaVG,beta,theta)
    i = complex(0,1);
    % Correction for martingality
    omega = (1/beta) * log(1 - theta*beta - 0.5*sigmaVG^2*beta);
    % Drift correction term
    mu = (r - q) + omega;
    % Levy process part
    psi = (1 - i*u*theta*beta + 0.5*sigmaVG^2*beta*u.^2) .^ (-t/beta); 
    % Chf for Variance Gamma model
    cf = exp(i*u*mu*t) .* psi;
end

% Heston Stochastic Volatility (Little Trap, see Albrecher et. al)
function cf = chf_heston(u,t,r,q,kappa,theta,eta,rho,v0)
    i = complex(0,1);
    % Functions D and g
    D = sqrt((-kappa + i*rho*eta.*u).^2 + (u.^2 + i*u)*eta^2 );
    g = (kappa - i*rho*eta*u - D) ./ (kappa - i*rho*eta*u + D);
    % Complex-valued functions A and C
    C = (1/eta^2) * (1-exp(-D*t))./(1-g.*exp(-D*t)) .* (kappa - eta*rho*i*u - D);
    A = i*u*(r-q)*t + kappa*theta/eta^2 * (kappa - eta*rho*i*u - D)*t... 
        - kappa*theta/eta^2 * 2*log((1-g.*exp(-D*t))./(1-g));
    % ChF for the Heston model
    cf = exp(A + C * v0);
end

% Bates Jump diffusion stochastic volatility
function cf = chf_bates(u,t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda)
    i = complex(0,1);
    % Functions D and g
    D = sqrt((-kappa + i*rho*eta.*u).^2 + (u.^2 + i*u)*eta^2 );
    g = (kappa - i*rho*eta*u - D) ./ (kappa - i*rho*eta*u + D);
    % Complex-valued functions A and C
    C = (1/eta^2) * (1-exp(-D*t))./(1-g.*exp(-D*t)) .* (kappa - eta*rho*i*u - D);
    A = i*u*(r-q)*t + kappa*theta/eta^2 * (kappa - eta*rho*i*u - D)*t... 
        - kappa*theta/eta^2 * 2*log((1-g.*exp(-D*t))./(1-g));
    % Adjustment for the Bates model
    A = A - lambda*i*u*t*(exp(muJ+1/2*sigmaJ^2)-1) + lambda*t*(exp(i*u*muJ-1/2*sigmaJ^2*u.^2)-1);
    % ChF for the Bates model
    cf = exp(A + C * v0);
end

% Heston Wishart Stochastic Volatility
function cf = chf_wmsv(u,t,r,q, ...
                       m11, m12, m21, m22, ...
                       r11, r12, r21, r22, ...
                       q11, q12, q21, q22, ...
                       Sig11, Sig12, Sig21, Sig22, ...
                       beta)
    M = [m11 m12; m21 m22];
    R = [r11 r12; r21 r22];
    Q = [q11 q12; q21 q22];
    Sigma0 = [Sig11, Sig12; Sig21 Sig22];
    i = complex(0,1);
    cf = zeros(1,length(u));   
    for k = 1:length(u) 
        w = u(1,k);
        MATEXP = expm(t*[M, -2.0*(Q')*Q; ...
                         0.5*((i*w)^2 - (i*w))*eye(2), ...
                         -1.0*((M') + 2.0*(i*w)*(R')*Q)]);
        A21 = [MATEXP(3:4,1:2)];
        A22 = [MATEXP(3:4,3:4)];
        % Computation of matrix function A(tau)
        A = (A22)^(-1) * A21;
        % Computation of scalar function C(tau)
        C = -0.5*beta * trace(logm(A22) + ((M') + 2.0*(i*w)*R*Q)*t) ...
            + (i*w)*(r-q)*t;
        % Characteristic Function for WMSV model
        cf(1,k) = exp(trace(A*Sigma0) + C);
    end
end