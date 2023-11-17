%************************************************************************%
% Here we implement we recover the density for various model             %
% using COS inversion                                                    %
%************************************************************************%
clear variables; close all; clc

model = 'WMSV';

% Market parameters
S0 = 100;            % Initial stock price (not used for density)
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield

t = 1.0;             % Time to maturity

if strcmp(model,'BlackScholes')
    % Parameters
    sigma = 0.25;        % Volatility
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigma);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigma);

elseif strcmp(model,'Merton')
    % Parameters
    sigma = 0.25;        % Volatility
    muJ = 0.02;          % Jump mean
    sigmaJ = 0.08;       % Jump standard deviation
    lambda = 8;          % Intensity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigma,muJ,sigmaJ,lambda);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigma,muJ,sigmaJ,lambda);

elseif strcmp(model,'Kou')
    % Parameters
    sigma = 0.25;        % Volatility
    p1 = 0.4;            % Probability of positive jump
    eta1 = 0.02;         % Positive jump size  
    eta2 = 0.08;         % Negative jump size
    lambda = 1;          % Intensity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigma,p1,eta1,eta2,lambda);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigma,p1,eta1,eta2,lambda);

elseif strcmp(model,'VG')
    % Parameters
    sigmaVG = 0.25;      % Volatility
    beta = 0.2;          % 
    theta = -0.14;       % 
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigmaVG,beta,theta);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigmaVG,beta,theta);

elseif strcmp(model,'Heston')
    % Parameters
    v0 = 0.0225;         % Initial variance  
    kappa = 1.5768;      % Speed of mean reversion 
    theta = 0.0398;      % Long-run average 
    eta = 0.5751;        % vol-of-vol 
    rho = -0.5711;       % Correlation 
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,kappa,theta,eta,rho,v0);
    % Compute cumulants
    c = getCumulants(model,t,r,q,kappa,theta,eta,rho,v0);

elseif strcmp(model,'Bates')
    % Parameters
    v0 = 0.0225;         % Initial variance  
    kappa = 1.5768;      % Speed of mean reversion 
    theta = 0.0398;      % Long-run average 
    eta = 0.5751;        % vol-of-vol 
    rho = -0.5711;       % Correlation 
    muJ = 0.02;          % Jump mean
    sigmaJ = 0.08;       % Jump standard deviation
    lambda = 8;          % Intensity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda);
    % Compute cumulants
    c = getCumulants(model,t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda);

elseif strcmp(model,'WMSV')
    % Parameters
    % Elements of (2x2) negative semidefinite matrix M in order to
    % ensure typical mean reversion
    m11 = -3.0; m12 = 0.0; m21 = 0.0; m22 = -3.0;
    % Elements of (2x2) correlation matrix
    r11 = -0.7; r12 = 0.0; r21 = 0.0; r22 = -0.7; 
    q11 = 0.25; q12 = 0.0; q21 = 0.0; q22 = 0.25; 
    % Elements of (2x2) positive semidefinite matrix. It represents
    % variance process initial state
    Sig11 = 0.01; Sig12 = 0.0; 
    Sig21 = 0.0; Sig22 = 0.01;  
    % Gindikin condition (same as role of Feller condition in scalar case)
    beta = 3;    
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                        m11, m12, m21, m22, ...
                                        r11, r12, r21, r22, ...
                                        q11, q12, q21, q22, ...
                                        Sig11, Sig12, Sig21, Sig22, ...
                                        beta);
    % Compute cumulants (via finite difference)
    c = getCumulants(model,t,r,q, ...
                     m11, m12, m21, m22, ...
                     r11, r12, r21, r22, ...
                     q11, q12, q21, q22, ...
                     Sig11, Sig12, Sig21, Sig22, ...
                     beta);

else
    disp('Method not found')
end

% Domain for the density f(x)
x = linspace(-1.5,1.5,1000);

% Define the COS method integration range
L = 12; N = 8192;
a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) );
b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) );  

% Reference value
fx = getRecoveredDensityByCOS(cf,x,N,a,b);

% Plot density recovery
figure(1)
plot(x,fx,'linewidth',1.5,'Marker','o')
x_txt = '$\bar{x}=\log(S_t)$';
y_txt = "$\phi(\bar{y}|\bar{x})$";
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
title_txt = sprintf('%s Density recovery in t = %4.1f with COS',model,t);
title(title_txt)
grid on

