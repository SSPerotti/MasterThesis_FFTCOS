% ------------------------------------------------------------------------
% Here we implement COS method to price european and digital call (put)
% option in the Bates [1996] Stochastic Volatility Jump Diffusion model
% ------------------------------------------------------------------------
clear variables; close all; clc

format long e

model = 'Bates';
product = 'European'; % 'European' or 'Digital'
type = 'Call';        % 'Call' or 'Put'

% Market Parameters
S0 = 100;            % Initial stock price
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield

% Bates parameters
v0 = 0.0175;         % Initial variance  
kappa = 1.5768;      % Speed of mean reversion 
theta = 0.0398;      % Long-run average 
eta = 0.5751;        % Vol-of-vol 
rho = -0.5711;       % Correlation 
muJ = 0.02;          % Jump mean  0.02
sigmaJ = 0.08;       % Jump standard deviation 0.08
lambda = 1;          % Intensity 8

t = 1.0;             % Time to maturity


% COS settings
L = 10;
K = (80:1:120)'; % vector of strikes
N = 2^14;


% The Chf does not include the coefficient "+iuX(t_0)" 
% as this coefficient is included internally in the evaluation
cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                    kappa,theta,eta,rho,v0, ...
                                    muJ,sigmaJ,lambda ...
                                    );

% Compute the cumulants for COS method
c = getCumulants(model,t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda);

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
xlabel('Strike, $K$','Interpreter','latex',FontSize=14)
ylabel('Option price','Interpreter','latex',Fontsize=14)
legend('COS Price')



% ------------------------------------------------------------------------
% COS price vs CONVOLUTION price comparison
% ------------------------------------------------------------------------
strike = 100;
Lcos = 6;
Npow = 3:1:10;    
Ngrid = floor(2.^Npow);

% CONV price
integrationBound = 16;
alpha = 0.5; 
optCONV = getCallPriceByCONV(integrationBound,alpha,S0,strike,t,r,cf);

% Timing results
NoOfIterations = 1;
error = zeros(length(Ngrid),1);
idx = 1;

for j = 1:length(Ngrid)
    tic
    for i = 1:NoOfIterations
        optCOS = getOptionPriceByCOS(cf,c,product,type,S0,r,t, ...
                                     strike,Ngrid(j),Lcos);
    end
    error(idx) = abs(optCONV - optCOS);
    time_elapsed = toc/i*1000;
    sprintf('For N = %.0f it took %f msec to price',Ngrid(j),time_elapsed)
    sprintf('For N = %.0f the error is equal to %e',Ngrid(j),error(idx))
    idx = idx + 1;
end

% Plot the results
figure(2); 
clf(2);
semilogy(Npow,error,'LineWidth',1.5,'LineStyle','-.', ...
         'Marker','*','MarkerSize',10);
hold on; grid on;
xlabel('$2^{N}$','Interpreter','latex',FontSize=14)
ylabel('Abs. Error','Interpreter','latex',FontSize=14)
title_txt = ['$|C_{\mathtt{conv}}^{\mathtt{eu}}(t,S,v) ' ...
             ' - C_{\mathtt{cos}}^{\mathtt{eu}}(t,S,v)|$'];
title(title_txt,'interpreter','latex',FontSize=14);



%-------------------------------------------------------------------------
% Option price, Delta and Gamma surface for Bates model
%-------------------------------------------------------------------------
tGrid = 0:1/52:1;
COSmat = zeros(length(K),length(tGrid));
DELTA = zeros(length(K),length(tGrid));
GAMMA = zeros(length(K),length(tGrid));

for j = 1:length(tGrid)
    chf = @(u) getCharacteristicFunction(model,u,tGrid(j),r,q, ...
                                         kappa,theta,eta,rho,v0, ...
                                         muJ,sigmaJ,lambda ...
                                         );
    ck = getCumulants(model,tGrid(j),r,q, ...
                      kappa,theta,eta,rho,v0, ...
                      muJ,sigmaJ,lambda ...
                      );
    [COSmat(:,j),DELTA(:,j),GAMMA(:,j)] = getOptionPriceDeltaGammaByCOS(...
                        chf,ck,product,type,S0,r,tGrid(j),K,N,L);
end

[Tmat,Kmat] = meshgrid(tGrid,K);

% Option price
figure(3)
clf(3);
surf(Tmat,Kmat,COSmat)
xlabel('Tempo, $\tau$','Interpreter','latex',FontSize=14)
ylabel('Strike, $K$','Interpreter','latex',Fontsize=14)

if strcmp(product,'European')
    if strcmp(type,'Call')
        zOptvaluetxt = 'Call europea: $C^{\mathtt{eu}}(t,S,v)$';
    elseif strcmp(type,'Put')
        zOptvaluetxt = 'Put europea: $P^{\mathtt{eu}}(t,S,v)$';
    end
elseif strcmp(product,'Digital')
    if strcmp(type,'Call')
        zOptvaluetxt = 'Call digitale: $C^{\mathtt{dig}}(t,S,v)$';
    elseif strcmp(type,'Put')
        zOptvaluetxt = 'Put digitale: $P^{\mathtt{dig}}(t,S,v)$';
    end
end
title(zOptvaluetxt,'interpreter','latex',FontSize=14);

if strcmp(type,'Call')
    view([-127.5 30]) % -127.5 call 
elseif strcmp(type,'Put')
    view([-307.5 30]) % -307.5 put
end

% Delta surface
figure(4) 
clf(4);
surf(Tmat,Kmat,DELTA)
xlabel('Tempo, $\tau$','Interpreter','latex',FontSize=14)
ylabel('Strike, $K$','Interpreter','latex',FontSize=14)
deltaTxt = '$\mathrm{Delta}: \Delta = \partial_{s} C(t,S,v)$';
title(deltaTxt,'interpreter','latex',FontSize=14);

% Gamma surface
figure(5)
clf(5);
surf(Tmat,Kmat,GAMMA)
xlabel('Tempo, $\tau$','Interpreter','latex',FontSize=14)
ylabel('Strike, $K$','Interpreter','latex',Fontsize=14)
gammaTxt = '$\mathrm{Gamma}: \Gamma = \partial_{ss} C(t,s,v)$';
title(gammaTxt,'interpreter','latex',FontSize=14);



