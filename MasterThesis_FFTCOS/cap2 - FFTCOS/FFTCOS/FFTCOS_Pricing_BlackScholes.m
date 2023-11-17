%************************************************************************%
% Here we implement COS method to price european and digital call (put)  %
% option in the Black-Scholes [1973] model                               %
%************************************************************************%
clear variables; close all; clc

format long e

model = 'BlackScholes';
product = 'European';  % 'European' or 'Digital'
type = 'Call';         % 'Call' or 'Put'

% Parameters
S0 = 100;            % Initial stock price 
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield
sigma = 0.25;        % Volatility 
t = 0.1;             % Time to maturity 


% COS settings
L = 10; N = 8192;
K = (80:1:120)'; % vector of strikes


% The Chf does not include the coefficient "+iuX(t_0)" 
% as this coefficient is included internally in the evaluation
cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigma);

% Compute the cumulants for COS method
c = getCumulants(model,t,r,q,sigma);

% Timing results
NoOfIterations = 1;
tic
for i = 1:NoOfIterations
    COSprice = getOptionPriceByCOS(cf,c,product,type,S0,r,t,K,N,L);
end
time_elapsed = toc;
sprintf('It took %f seconds to price',time_elapsed/NoOfIterations)


% Plot the results
if strcmp(product,'European')
    BSprice = getEuropeanOptionPriceByBS(type,S0,K,sigma,t,r);
elseif strcmp(product,'Digital')
    BSprice = getDigitalOptionPriceByBS(type,S0,K,sigma,t,r);
end

figure(1)
clf(1)
plot(K,COSprice,'LineWidth',2.0)
hold on
plot(K,BSprice,'LineWidth',1.0,'Color','r',Marker='o')
grid on
xlabel('Strike, $K$','Interpreter','latex',FontSize=14)
ylabel('Option price','Interpreter','latex',Fontsize=14)
legend('COS price','BS price','Location','best')

% Error computation
error = zeros(1,length(K));

for i = 1:length(K)
    error(i) = abs(COSprice(i) - BSprice(i));
    sprintf('Abs error for strike %.2f is equal to %.2E',K(i),error(i))
end


% =======================================================================
% COS price vs BLACK SCHOLES price comparison
% =======================================================================
strike = 80;  
Npow = 2:1:6; 
Ngrid = floor(2.^Npow);

% Closed-form price expression
if strcmp(product,'European')
    valBS = getEuropeanOptionPriceByBS(type,S0,strike,sigma,t,r);
elseif strcmp(product,'Digital')
    valBS = getDigitalOptionPriceByBS(type,S0,strike,sigma,t,r);
end

% Timing results
NoOfIterations = 1;
error = zeros(length(Ngrid),1);
idx = 1;

for j = Ngrid
    tic
    for i = 1:NoOfIterations
        valCOS = getOptionPriceByCOS(cf,c,product,type,S0,r,t,strike,j,L);
    end
    error(idx) = abs(valBS - valCOS);
    time_elapsed = toc/i*1000;
    sprintf('For N = %.0f it took %f seconds to price',j,time_elapsed)
    sprintf('For N = %.0f the error is equal to %e',j,error(idx))
    idx = idx +1;
end


% Plot the results
figure(2); 
clf(2);
semilogy(Npow,error,'LineWidth',1.5,'LineStyle','-.', ...
         'Marker','*','MarkerSize',10);
hold on; grid on;
xlabel('$2^{N}$','Interpreter','latex',FontSize=14)
ylabel('Abs. Error','Interpreter','latex',FontSize=14)
title_txt = ['$|C_{\mathtt{bs}}^{\mathtt{eu}}(t,S) ' ...
             ' - C_{\mathtt{cos}}^{\mathtt{eu}}(t,S)|$'];
title(title_txt,'interpreter','latex',FontSize=14);



% ------------------------------------------------------------------------
% Option Price surface
% ------------------------------------------------------------------------
tGrid = 0:1/52:1;
COSmat = zeros(length(K),length(tGrid));

for j = 1:length(tGrid)
        COSmat(:,j) = getOptionPriceByCOS(cf,c,product,type, ...
                                          S0,r,tGrid(j),K,N,L);
end

[Tmat,Kmat] = meshgrid(tGrid,K);

figure(3)
clf(3)
surf(Tmat,Kmat,COSmat)
xlabel('Tempo, $\tau$','Interpreter','latex',FontSize=14)
ylabel('Strike, $K$','Interpreter','latex',Fontsize=14)

if strcmp(product,'European')
    if strcmp(type,'Call')
        zOptvaluetxt = 'Call europea: $C^{\mathtt{eu}}(t,S)$';
    elseif strcmp(type,'Put')
        zOptvaluetxt = 'Put europea: $P^{\mathtt{eu}}(t,S)$';
    end
elseif strcmp(product,'Digital')
    if strcmp(type,'Call')
        zOptvaluetxt = 'Call digitale: $C^{\mathtt{dig}}(t,S)$';
    elseif strcmp(type,'Put')
        zOptvaluetxt = 'Put digitale: $P^{\mathtt{dig}}(t,S)$';
    end
end
title(zOptvaluetxt,'interpreter','latex',FontSize=14);

if strcmp(type,'Call')
    view([-127.5 30])  % -127.5 call 
elseif strcmp(type,'Put')
    view([-307.5 30])  % -307.5 put
end



%************************************************************************%
% Here we have closed formula for Black Scholes model                    %
%                                                                        %
% - For European option we have                                          %
%       getEuropeanOptionPriceByBlackScholes(type,S0,K,sigma,T,r)        %
%                                                                        %
% - For Digital option we have                                           %
%       getDigitalOptionPriceByBlackScholes(type,S0,K,sigma,T,r)         %
%************************************************************************%

function resu = getEuropeanOptionPriceByBS(type,S0,K,sigma,t,r)
    d1 = (log(S0 ./ K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t));
    d2 = d1 - sigma * sqrt(t);
    if strcmp(type,'Call') 
        resu = normal(d1,0,1) * S0 - normal(d2,0,1) .* K * exp(-r * t);
    elseif strcmp(type,'Put')
        resu = normal(-d2,0,1) .* K*exp(-r*t) - normal(-d1,0,1)*S0;
    end
end

function resu = getDigitalOptionPriceByBS(type,S_0,K,sigma,t,r)
    d1 = (log(S_0 ./ K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t));
    d2 = d1 - sigma * sqrt(t);
    if strcmp(type,'Call') 
        resu = K * exp(-r*t) .* normal(d2,0,1) ;
    elseif strcmp(type,'Put')
        resu = K * exp(-r*t) .* (1-normal(d2,0,1));
    end
end

% Erf(x) function
function resu = normal(x,mu,sigma)
    z = (x - mu) ./ sigma;
    resu = 0.5 * erfc(-z ./ sqrt(2));
end
