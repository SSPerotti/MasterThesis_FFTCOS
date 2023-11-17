%************************************************************************%
% Here we study the effects of the jump part of parameter vector         %
% that parametrizes Bates model                                          %
%                                                                        %
%************************************************************************%
clear variables; close all; clc

model = 'Bates';

% Market Parameters
S0 = 100;            % Initial stock price
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield

% Bates parameters
v0 = 0.0225;         % Initial variance  
kappa = 1.5768;      % Speed of mean reversion 
theta = 0.0398;      % Long-run average 
eta = 0.5751;        % Vol-of-vol 
rho = -0.5711;       % Correlation 
muJ = 0.02;          % Jump mean
sigmaJ = 0.08;       % Jump standard deviation
lambda = 8;          % Intensity

t = 0.1;             % Time to maturity

x_txt = '$\bar{x}=\log(S_t)$';
y_txt = '$\phi_{\mathtt{Bates}}(\bar{y}|\bar{x})$';

% COS method settings
N = 4096;
L = 10;


% ===========================================================
% Effect of [muJ] on Bates density
% ===========================================================
muJGrid = [0.02 0.04 0.06 0.10];
time = 2.0;
x = linspace(-1.0,1.0,1000);

figure(1)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(muJGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,theta,eta,rho,v0, ...
                                        muJGrid(j),sigmaJ,lambda ...
                                        );
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q, ...
                     kappa,theta,eta,rho,v0, ...
                     muJGrid(j),sigmaJ,lambda ...
                     );
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) );
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\mu_{J}=',num2str(muJGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \mu_J on Bates density')


% ===========================================================
% Effect of [sigmaJ] on Bates density
% ===========================================================
lambdaGrid = [0.05 0.10 0.15 0.25];
time = 2.0;
x = linspace(-1.0,1.0,1000);

figure(2)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(lambdaGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,theta,eta,rho,v0, ...
                                        muJ,lambdaGrid(j),lambda ...
                                        );
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q, ...
                     kappa,theta,eta,rho,v0, ...
                     muJ,lambdaGrid(j),lambda ...
                     );
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) );
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\sigma_{J}^{2}=',num2str(lambdaGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \sigma_J on Bates density')


% ===========================================================
% Effect of [lambda] on Bates density
% ===========================================================
lambdaGrid = [4 8 10 15];
time = 2.0;
x = linspace(-1.0,1.0,1000);

figure(3)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(lambdaGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,theta,eta,rho,v0, ...
                                        muJ,sigmaJ,lambdaGrid(j) ...
                                        );
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q, ...
                     kappa,theta,eta,rho,v0, ...
                     muJ,sigmaJ,lambdaGrid(j) ...
                     );
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) );
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\lambda=',num2str(lambdaGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \lambda on Bates density')