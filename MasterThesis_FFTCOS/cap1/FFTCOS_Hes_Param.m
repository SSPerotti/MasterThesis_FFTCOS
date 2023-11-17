%************************************************************************%
% Here we study the effects of the parameter vector                      %
% that parametrizes Heston model                                         %
%                                                                        %
%************************************************************************%
clear variables; close all; clc

model = 'Heston';

% Market parameters
S0 = 100;            % Initial stock price
r = 0;               % Risk free rate
q = 0;               % Dividend yield

% Heston parameters
v0 = 0.0175;         % Initial variance  
kappa = 1.5768;      % Speed of mean reversion 
theta = 0.0398;      % Long-run average 
eta = 0.5751;        % vol-of-vol 
rho = -0.5711;       % Correlation 

t = 0.1;             % Time to maturity

x_txt = '$\bar{x}=\log(S_t)$';
y_txt = '$\phi_{\mathtt{Hes}}(\bar{y}|\bar{x})$';

% COS method settings
N = 4096;
L = 8;

% ===========================================================
% Effect of [kappa] on Heston density
% ===========================================================
kappaGrid = [0.70 0.90 1.35 1.75];
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

for j = 1:length(kappaGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappaGrid(j),theta,eta,rho,v0);
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q,kappaGrid(j),theta,eta,rho,v0);
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\kappa=',num2str(kappaGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \kappa on Heston density')


% ===========================================================
% Effect of [theta] on Heston density
% ===========================================================
thetaGrid = [0.03 0.7 0.15 0.25];
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

for j = 1:length(thetaGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,thetaGrid(j),eta,rho,v0);
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q,kappa,thetaGrid(j),eta,rho,v0);
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\theta=',num2str(thetaGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \theta on Heston density')


% ===========================================================
% Effect of [eta] on Heston density
% ===========================================================
etaGrid = [0.10 0.30 0.50 0.90];
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

for j = 1:length(etaGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,theta,etaGrid(j),rho,v0);
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q,kappa,theta,etaGrid(j),rho,v0);
    % Define the COS method integration range 
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\eta=',num2str(etaGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \eta on Heston density')


% ===========================================================
% Effect of [rho] on Heston density
% ===========================================================
rhoGrid = [0.30 0 -0.30 -0.70];
time = 2.0;
x = linspace(-1.0,1.0,1000);

figure(4)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(etaGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,theta,eta,rhoGrid(j),v0);
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q,kappa,theta,eta,rhoGrid(j),v0);
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$\rho=',num2str(rhoGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of \rho on Heston density')


% ===========================================================
% Effect of [v0] on Heston density
% ===========================================================
v0Grid = [0.01 0.05 0.09 0.16];
time = 2.0;
x = linspace(-1.0,1.0,1000);

figure(5)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(v0Grid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        kappa,theta,eta,rho,v0Grid(j));
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q,kappa,theta,eta,rho,v0Grid(j));
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$v_{0}=',num2str(v0Grid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of v_0 on Heston density')