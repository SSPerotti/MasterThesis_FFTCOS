% ------------------------------------------------------------------------
% Here we implement COS method to price european and digital call (put)
% option for various model
% ------------------------------------------------------------------------
clear variables; close all ; clc

model = 'WMSV';
product = 'European';   % 'European' or 'Digital'
type = 'Call';          % 'Call' or 'Put'

% Market parameters
S0 = 100;            % Initial stock price
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
    muJ = 0.10;          % Jump mean
    sigmaJ = 0.25;       % Jump standard deviation
    lambda = 8;          % Intensity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                        sigma,muJ,sigmaJ,lambda);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigma,muJ,sigmaJ,lambda);

elseif strcmp(model,'Kou')
    % Parameters
    sigma = 0.25;        % Volatility
    p1 = 0.4;            % Probability of positive jump
    eta1 = 10;           % Positive jump size  
    eta2 = 5;            % Negative jump size
    lambda = 8;          % Intensity
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                        sigma,p1,eta1,eta2,lambda);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigma,p1,eta1,eta2,lambda);

elseif strcmp(model,'VG')
    % Parameters
    sigmaVG = 0.25;      % Volatility
    beta = 0.2;          % 
    theta = -0.14;       % 
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigmaVG,beta,theta,S0);
    % Compute cumulants
    c = getCumulants(model,t,r,q,sigmaVG,beta,theta);
    %#- use L = 150 for Variance Gamma

elseif strcmp(model,'Heston')
    % Parameters
    v0 = 0.0225;         % Initial variance  
    kappa = 1.5768;      % Speed of mean reversion 
    theta = 0.0398;      % Long-run average 
    eta = 0.5751;        % vol-of-vol 
    rho = -0.5711;       % Correlation 
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                        kappa,theta,eta,rho,v0);
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
    cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                        kappa,theta,eta,rho,v0, ...
                                        muJ,sigmaJ,lambda);
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
                                        beta ...
                                        );
    % Compute cumulants (via finite difference)
    c = getCumulants(model,t,r,q, ...
                     m11, m12, m21, m22, ...
                     r11, r12, r21, r22, ...
                     q11, q12, q21, q22, ...
                     Sig11, Sig12, Sig21, Sig22, ...
                     beta ...
                     );

else
    disp('Method not found')
end

% Define the COS method settings
N = 2^12; 
L = 10; 
K = (40:1:160)';


% Timing results
NoOfIterations = 1;
tic
for i = 1:NoOfIterations
    COSprice = getOptionPriceByCOS(cf,c,product,type,S0,r,t,K,N,L);
end
time_elapsed = toc;
sprintf('It took %f seconds to price',time_elapsed/NoOfIterations)

% Plot the results
figure(1)
plot(K,COSprice,'LineWidth',2.0)
grid on
xlabel('Strike, $K$','Interpreter','latex',FontSize=12)
ylabel('Option price','Interpreter','latex',Fontsize=12)
legend('COS Price')

