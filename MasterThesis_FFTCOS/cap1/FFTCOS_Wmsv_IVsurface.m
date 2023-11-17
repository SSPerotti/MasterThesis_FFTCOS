clear variables; close all; clc

format short

model = 'WMSV';
product = 'European';
type = 'Call';

% Parameters
S0 = 100;            % Initial stock price (not used for density)
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield
t = 1.0;             % Time to maturity

% Elements of (2x2) negative semidefinite matrix M in order to
% ensure typical mean reversion
m11 = -3.0; m12 = 0.0; 
m21 = 0.0;  m22 = -3.0;

% Elements of (2x2) correlation matrix
r11 = -0.7; r12 = 0.0; 
r21 = 0.0;  r22 = -0.7; 

% Elements of (2x2) matrix Q
q11 = 0.25; q12 = 0.0; 
q21 = 0.0;  q22 = 0.25; 

% Elements of (2x2) positive semidefinite matrix. 
% It represents variance process initial state
Sig11 = 0.01; Sig12 = 0.0; 
Sig21 = 0.0;  Sig22 = 0.01;  

% Gindikin condition (same as role of Feller condition in scalar case)
beta = 3;  

% Range of strike prices
K = linspace(40,180,25)';
figure(1)
hold on;
grid on;

% COS method settings
N = 4096;
L = 8;
TMat = linspace(90/360,720/360,24);
IV = zeros(length(TMat),length(K));

for idx = 1:1:length(TMat)
    Ttemp = TMat(idx);
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,Ttemp,r,q, ...
                                        m11, m12, m21, m22, ...
                                        r11, r12, r21, r22, ...
                                        q11, q12, q21, q22, ...
                                        Sig11, Sig12, Sig21, Sig22, ...
                                        beta);
    valCOS = getOptionPriceByCOS(cf,product,type,S0,r,t,K,N,L);
    figure(1)
    plot(K,valCOS)
    hold on
    for idy = 1:length(K)
        IV(idx,idy) = impliedVola(S0,K(idy),r,Ttemp,valCOS(idy),q)*100;
    end
end

figure(2)
surf(log(S0./K),TMat,IV)
xlabel('Moneyness, $S/K$','Interpreter','latex',FontSize=12)
ylabel('Maturity, $T$','Interpreter','latex',FontSize=12)
zlabel('Implied Volatility, $\sigma(T,K)$','Interpreter','latex',FontSize=12)