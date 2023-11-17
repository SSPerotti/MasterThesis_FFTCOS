%% Pricing of European call and put options with the COS method
clear variables
close all
clc

format long

model = 'BlackScholes';
product = 'Digital';
type = 'Put';

% Parameters
S0 = 100;            % Initial stock price 
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield
sigma = 0.25;        % Volatility
t = 0.1;             % Time to maturity


% COS settings
L = 8;
K = (80:1:120)'; % vector of strikes
N = 2^16;


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
    BSprice = getEuropeanOptionPriceByBlackScholes(type,S0,K,sigma,t,r);
elseif strcmp(product,'Digital')
    BSprice = getDigitalOptionPriceByBlackScholes(type,S0,K,sigma,t,r);
end

figure(1)
clf(1)
plot(K,COSprice,'LineWidth',2.0)
hold on
plot(K,BSprice,'LineWidth',1.0,'Color','r',Marker='o')
grid on
xlabel('Strike, $K$','Interpreter','latex',FontSize=12)
ylabel('Option price','Interpreter','latex',Fontsize=12)
legend('COS price','BS price','Location','best')

% Error computation
error = zeros(1,length(K));

for i = 1:length(K)
    error(i) = abs(COSprice(i)-BSprice(i));
    sprintf('Abs error for strike %.2f is equal to %.2E',K(i),error(i))
end

% ------------------------------------------------------------------------
%                        Option Price surface
% ------------------------------------------------------------------------
tGrid = 0:1/52:1;
COS_price = zeros(length(K),length(tGrid));

for j = 1:length(tGrid)
        COS_price(:,j) = getOptionPriceByCOS(cf,c,product,type,S0,r,tGrid(j),K,N,L);
end

figure(2)
[X,Y] = meshgrid(tGrid,K);
surf(X,Y,COS_price)
xlabel('Tempo, $t$','Interpreter','latex',FontSize=12)
ylabel('Strike, $K$','Interpreter','latex',Fontsize=12)
zlabel('Option price','Interpreter','latex',Fontsize=12)
if strcmp(type,'Call')
    view([-127.5 30]) % -127.5 call | -307.5 put
elseif strcmp(type,'Put')
    view([-307.5 30])
end
% colorbar


% ------------------------------------------------------------------------
%                        Error for various N
% ------------------------------------------------------------------------
KK = 80;
Npow = 2:1:6;
Ngrid = 2.^Npow;

% Closed-form price expression
if strcmp(product,'European')
    valBS = getEuropeanOptionPriceByBlackScholes(type,S0,KK,sigma,t,r);
elseif strcmp(product,'Digital')
    valBS = getDigitalOptionPriceByBlackScholes(type,S0,KK,sigma,t,r);
end

% Timing results
NoOfIterations = 1;
errorVec = zeros(length(Ngrid),1);
idx = 1;

for j = Ngrid
    tic
    for i = 1:NoOfIterations
        valCOS = getOptionPriceByCOS(cf,c,product,type,S0,r,t,KK,j,L);
    end
    errorVec(idx)= abs(valBS - valCOS);
    time_elapsed = toc/i*1000;
    sprintf('For N= %.0f it took %f seconds to price',j,time_elapsed)
    sprintf('For N= %.0f the error is equal to %e',j,errorVec(idx))
    idx = idx +1;
end


% Plot the results
figure(3); 
clf(3);
% plot(Ngrid,errorVec,'LineWidth',1.5,'LineStyle','-.');
semilogy(Npow,errorVec,'LineWidth',1.5,'LineStyle','-.');
hold on; grid on;
xlabel('$2^{N}$','Interpreter','latex',FontSize=12)
ylabel('$\log_{10}(|\mathrm{errore}|)$','Interpreter','latex',FontSize=12)


% ------------------------------------------------------------------------
% Here we have closed formula for Black Scholes model
%
% - For European option we have
%       getEuropeanOptionPriceByBlackScholes(type,S0,K,sigma,T,r)
% 
% - For Digital option we have
%       getDigitalOptionPriceByBlackScholes(type,S0,K,sigma,T,r)
% ------------------------------------------------------------------------

function value = getEuropeanOptionPriceByBlackScholes(type,S0,K,sigma,t,r)
    d1 = (log(S0 ./ K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t));
    d2 = d1 - sigma * sqrt(t);
    if strcmp(type,'Call') 
        value = normcdf(d1) * S0 - normcdf(d2) .* K * exp(-r * t);
    elseif strcmp(type,'Put')
        value = normcdf(-d2) .* K*exp(-r*t) - normcdf(-d1)*S0;
    end
end

% Closed-form expression of cash-or-nothing option with the Black-Scholes formula
function value = getDigitalOptionPriceByBlackScholes(type,S_0,K,sigma,t,r)
    d1 = (log(S_0 ./ K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t));
    d2 = d1 - sigma * sqrt(t);
    if strcmp(type,'Call') 
        value = K * exp(-r*t) .* normcdf(d2) ;
    elseif strcmp(type,'Put')
        value = K * exp(-r*t) .* (1-normcdf(d2));
    end
end
