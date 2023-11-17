%% Density recovery with COS method
clear variables
close all 
clc

format long e 

model = 'Merton';

if strcmp(model,'BlackScholes')
    % Parameters
    S0 = 100;            % Initial stock price
    r = 0.05;            % Risk free rate 
    q = 0;               % Dividend yield
    sigma = 0.25;        % Volatility
    t = 0.1;             % Time to maturity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigma); % model,u,t,r,q,param
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigma); % model,t,r,q,param

elseif strcmp(model,'Merton')
    % Parameters
    S0 = 100;            % Initial stock price (not used for density)
    r = 0.05;            % Risk free rate 
    q = 0;               % Dividend yield
    sigma = 0.25;        % Volatility
    muJ = 0.10;          % Jump mean
    sigmaJ = 0.25;       % Jump standard deviation
    lambda = 8;          % Intensity
    t = 0.1;             % Time to maturity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,[sigma,muJ,sigmaJ,lambda]);
    % Compute cumulants
    c = getCumulants(model,t,r,q,[sigma,muJ,sigmaJ,lambda]);

elseif strcmp(model,'Kou')
    % Parameters
    S0 = 100;            % Initial stock price (not used for density)
    r = 0.05;            % Risk free rate 
    q = 0;               % Dividend yield
    sigma = 0.25;        % Volatility
    p1 = 0.4;            % Probability of positive jump
    eta1 = 10;           % Positive jump size  
    eta2 = 5;            % Negative jump size
    lambda = 8;          % Intensity
    t = 0.1;             % Time to maturity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,[sigma,p1,eta1,eta2,lambda]);
    % Compute cumulants
    c = getCumulants(model,t,r,q,[sigma,p1,eta1,eta2,lambda]);

elseif strcmp(model,'VG')
    % Parameters
    S0 = 100;            % Initial stock price (not used for density)
    r = 0.05;            % Risk free rate 
    q = 0;               % Dividend yield
    sigmaVG = 0.25;      % Volatility
    beta = 0.2;          % TODO
    theta = -0.14;       % TODO
    t = 0.1;             % Time to maturity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,[sigmaVG,beta,theta,S0]);
    % Compute cumulants
    c = getCumulants(model,t,r,q,[sigmaVG,beta,theta]);
    %#- use L = 150 for Variance Gamma

elseif strcmp(model,'Heston')
    % Parameters
    S0 = 100;            % Initial stock price (not used for density)
    r = 0.05;            % Risk free rate 
    q = 0;               % Dividend yield
    v0 = 0.0225;         % Initial variance  
    kappa = 1.5768;      % Speed of mean reversion 
    theta = 0.0398;      % Long-run average 
    eta = 0.5751;        % vol-of-vol 
    rho = -0.5711;       % Correlation 
    t = 0.1;             % Time to maturity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,[kappa,theta,eta,rho,v0]);
    % Compute cumulants
    c = getCumulants(model,t,r,q,[kappa,theta,eta,rho,v0]);

elseif strcmp(model,'Bates')
    % Parameters
    S0 = 100;            % Initial stock price (not used for density)
    r = 0.05;            % Risk free rate 
    q = 0;               % Dividend yield
    v0 = 0.0225;         % Initial variance  
    kappa = 1.5768;      % Speed of mean reversion 
    theta = 0.0398;      % Long-run average 
    eta = 0.5751;        % vol-of-vol 
    rho = -0.5711;       % Correlation 
    muJ = 0.02;          % Jump mean
    sigmaJ = 0.08;       % Jump standard deviation
    lambda = 8;          % Intensity
    t = 0.1;             % Time to maturity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,[kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda]);
    % Compute cumulants
    c = getCumulants(model,t,r,q,[kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda]);

elseif strcmo(model,'WMSV')
    % Parameters
    % Characteristic function
    % Compute cumulants
    disp('TODO')

else
    disp('Method not found')
end


% Domain for the density f(x)
x = linspace(-1.0,1.0,1000);

% Define the COS method integration range
L = 10; 
a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4)) ) ); 
% a = - L * sqrt(t)   % alternative
b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4)) ) ); 
% b = + L * sqrt(t);  % alternative

% Reference value (N = 2^16)
f_XExact = getRecoveredDensityByCOS(cf,x,2^16,a,b);

% Plot density recovery
figure(1)
plot(x,f_XExact,'linewidth',1.5,'Marker','o')
xlabel('$\bar{x}=\log(S_t)$',Interpreter='latex',FontSize=12)
ylabel("$\phi(\bar{y}|\bar{x})$",Interpreter='latex',FontSize=12)
title_txt = sprintf('model = %s Density recovery in t =%4.1f with COS method',model,t);
title(title_txt)
grid on
