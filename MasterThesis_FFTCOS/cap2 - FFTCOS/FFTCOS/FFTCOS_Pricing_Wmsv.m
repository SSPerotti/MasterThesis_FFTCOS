%************************************************************************%
% Here we implement COS method to price european and digital call (put)
% option in the Wishart Heston model [WMSV]
%
% Literature: see Da Fonseca, Grasselli and Tebaldi [2008]
%             "A Multifactor Volatility Heston Model"
%************************************************************************%
clear variables; close all; clc

format long

model = 'WMSV';
product = 'European';  % 'European' or 'Digital'
type = 'Call';         % 'Call' or 'Put'
t = 1.0;               % Time to maturity

% Market Parameters
S0 = 100;            % Initial stock price
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield

% WMSV parameters
% Matrix M
m11 = -0.5856; m12 = 0.0; 
m21 = 0.0;     m22 = -1.1856;   
% Matrix R
r11 = -0.2123; r12 = 0.0; 
r21 = 0.0;     r22 = -0.7562;     
% Matrix Q
q11 = 0.1752; q12 = 0.0; 
q21 = 0.0;    q22 = 0.1567;    
% Matrix Sigma0
Sig11 = 0.0192; Sig12 = 0.0; 
Sig21 = 0.0;    Sig22 = 0.1394;  
% Gindikin condition
beta = 1.1423; 

% COS settings
L = 10;
N = 8196;
K = (80:1:120)'; % vector of strikes

% Characteristic function 
cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                    m11, m12, m21, m22, ...
                                    r11, r12, r21, r22, ...
                                    q11, q12, q21, q22, ...
                                    Sig11, Sig12, Sig21, Sig22, ...
                                    beta);

% Compute the cumulants for COS method
c = getCumulants(model,t,r,q, ...
                 m11, m12, m21, m22, ...
                 r11, r12, r21, r22, ...
                 q11, q12, q21, q22, ...
                 Sig11, Sig12, Sig21, Sig22, ...
                 beta);

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


% ========================================================================
% COS price vs CARR MADAN price comparison
% ========================================================================
strike = 100;
Lcos = 6;
Npow = 3:1:10;    
Ngrid = floor(2.^Npow);

% Carr-Madan reference value 
optCM = getCallPriceByCarrMadan(S0,strike,r,t,cf);

% Timing results
NoOfIterations = 1;
error = zeros(length(Ngrid),1);
idx = 1;

for j = 1:length(Ngrid)
    tic
    for i = 1:NoOfIterations
        optCOS = getOptionPriceByCOS(cf,c,product,type,S0,r,...
                                     t,strike,Ngrid(j),Lcos);
    end
    error(idx) = abs(optCM - optCOS);
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
title_txt = ['$|C_{\mathtt{carr}}^{\mathtt{eu}}(t,S,v) ' ...
             ' - C_{\mathtt{cos}}^{\mathtt{eu}}(t,S,v)|$'];
title(title_txt,'interpreter','latex',FontSize=14);
%}


% ========================================================================
% Option price, Delta and Gamma surface for WMSV model
% ========================================================================
tGrid = 0:1/52:1;
COSmat = zeros(length(K),length(tGrid));
DELTA = zeros(length(K),length(tGrid));
GAMMA = zeros(length(K),length(tGrid));

for j = 1:length(tGrid)
    chf = @(u) getCharacteristicFunction(model,u,tGrid(j),r,q, ...
                                         m11, m12, m21, m22, ...
                                         r11, r12, r21, r22, ...
                                         q11, q12, q21, q22, ...
                                         Sig11, Sig12, Sig21, Sig22, ...
                                         beta);
    ck = getCumulants(model,t,r,q, ...
                      m11, m12, m21, m22, ...
                      r11, r12, r21, r22, ...
                      q11, q12, q21, q22, ...
                      Sig11, Sig12, Sig21, Sig22, ...
                      beta);
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