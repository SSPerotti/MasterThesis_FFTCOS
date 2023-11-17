clear variables; close all; clc

model = 'Bates';

% Option parameters
S0 = 150;
product = 'European'; 
type = 'Call';    

% Heston model part
kappa = 1.2;
theta = 0.05;
eta = 0.8;
rho = -0.75;
v0 = 0.05;
r = 0;
q = 0;

% Bates model
muJ = 0.0;
sigmaJ = 0.0;
lambda = 0.2;


% Range of strike prices
K = linspace(100,200,15)';
figure(1)
hold on;
grid on;
xlabel('strike, $K$','Interpreter','latex',FontSize=12)
ylabel('Implied Volatility, $\sigma(T,K)$','Interpreter','latex',FontSize=12)

% COS method settings
N = 4096;
L = 8;
TMat = linspace(1/12,2,20);
IV = zeros(length(TMat),length(K));

for idx = 1:1:length(TMat)
    Ttemp = TMat(idx);
    cf = @(u) getCharacteristicFunction(model,u,Ttemp,r,q, ...
                                        kappa, theta, eta, rho, v0, ...
                                        muJ, sigmaJ, lambda);
    valCOS = getOptionPriceByCOS(cf,product,type,S0,r,Ttemp,K,N,L);
    figure(1)
    plot(K,valCOS)
    hold on
    for idy = 1:length(K)
        IV(idx,idy) = ImpliedVolatility(type,valCOS(idy),K(idy), ...
                                        Ttemp,S0,r,0.3)*100;
    end
end

figure(2)
surf(log(S0./K),TMat,IV)
ylabel('Maturity, $\tau$','Interpreter','latex',FontSize=12)
xlabel('Log-Moneyness','Interpreter','latex',FontSize=12)
zlabel('Implied Volatility, $\sigma(T,K)$', ...
       'Interpreter','latex',FontSize=12)


% ------------------------------------------------------------------------
% Closed-form expression of European call/put option 
% with Black-Scholes formula
% ------------------------------------------------------------------------
function resu = getEuropeanOptionPriceByBS(type,S0,K,sigma,T,r)
    d1 = (log(S0 ./ K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    if strcmp(type,'Call') 
        resu = normal(d1,0,1) * S0 - normal(d2,0,1) .* K * exp(-r * T);
    elseif strcmp(type,'Put')
        resu = normal(-d2,0,1) .* K*exp(-r*T) - normal(-d1,0,1)*S0;
    end
end


function impliedVol = ImpliedVolatility(type,marketPrice,K,T,S0,r,initialVol)
    func = @(sigma) (getEuropeanOptionPriceByBS(type,S0,K,sigma,T,r) ...
                     - marketPrice).^1.0;
    impliedVol = fzero(func,initialVol);
end

% Erf(x) function
function resu = normal(x,mu,sigma)
    z = (x - mu) ./ sigma;
    resu = 0.5 * erfc(-z ./ sqrt(2));
end