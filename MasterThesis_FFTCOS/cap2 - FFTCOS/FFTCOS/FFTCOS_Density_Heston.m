%************************************************************************%
% Here we recover the density for Heston [1993] Stochastic Volatility    %
%                                                                        % 
%************************************************************************%
clear variables; close all; clc

format long e

model = 'Heston';

% Market parameters
S0 = 100;            % Initial stock price (not used for density)
r = 0;               % Risk free rate
q = 0;               % Dividend yield

% Heston parameters
v0 = 0.0175;         % Initial variance  
kappa = 1.5768;      % Speed of mean reversion 
theta = 0.0398;      % Long-run average 
eta = 0.5751;        % vol-of-vol 
rho = -0.5711;       % Correlation 

t = 0.1;             % Time to maturity

% Check if Feller condition is satisfied
if 2*kappa*theta >= eta*eta
    disp("Feller condition is satisfied")
else
    disp("Feller condition is not satisfied")
end

% Domain for the density f(x)
% use [-0.5, 0.5] for t = 0.1 
% use [-1.0, 1.0] for t = 1.0
% use [-5.0, 5.0] for t = 10
x = linspace(-0.5,0.5,1000); 


% ===========================================================
%                   Model density recovery 
% ===========================================================
% Characteristic function 
cf = @(u) getCharacteristicFunction(model,u,t,r,q, ...
                                    kappa,theta,eta,rho,v0);

% Cumulants needed for the integration range
c = getCumulants(model,t,r,q,kappa,theta,eta,rho,v0); 

% Define the COS method integration range
L = 12; N = 8192;
a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 

% Reference value (N = 2^13)
f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);

figure(1)
clf(1)
plot(x,f_XExact,'linewidth',1.5,'Marker','o')
x_txt = '$\bar{x}=\log(S_t)$';
xlabel(x_txt,Interpreter='latex',FontSize=14)
y_txt = '$\phi_{\mathtt{Hes}}(\bar{y}|\bar{x})$';
ylabel(y_txt,Interpreter='latex',FontSize=14)
title_txt = sprintf(['%s Density recovery ' ...
                     'in t = %4.1f with COS'],model,t);
title(title_txt)
grid on


% ===========================================================
%             Effect of [N] in density recovery 
% ===========================================================
% Iterate over different numbers of expansion terms
Npow = 3:0.5:10;    
Ngrid = floor(2.^Npow);
NoOfIterations = 1;
idx = 1;
error = zeros(length(Ngrid),1);

figure(2)
clf(2)
hold on; grid on
xlabel(x_txt,Interpreter='latex',FontSize=14)
ylabel(y_txt,Interpreter='latex',FontSize=14)
legend_txt = cell(length(Ngrid),1);
time_txt = cell(length(Ngrid),1);

for N = Ngrid
    % Density from the COS method
    tic
    for j = 1:NoOfIterations % 10^3 experiments
        f_X = getRecoveredDensityByCOS(cf,x,N,a,b);
    end
    timeEnd = toc/j * 1000;
    toc
    % Error
    error(idx) = max(f_XExact-f_X);
    sprintf('For N = %.0f the error is equal to %e',N,error(idx))
    %{
    plot(x,f_X,'linewidth',1.5,'LineStyle',LINESTYLE(idx), ...
         'Marker',MARKERS(idx),'MarkerIndices',1:50:length(x))
    %}
    plot(x,f_X,'linewidth',1.5)
    legend_txt{idx} = sprintf('N = %i',N);
    time_txt{idx} = sprintf('time neeeded for evaluation = %.2f miliseconds',timeEnd);
    idx = idx + 1;
end

legend(legend_txt)
% error
% time_txt
%{
title_txt = sprintf(['Effect of N in %s density ' ...
                     'recovery in t =%4.1f with COS'],model,t);
title(title_txt)
%}



% ===========================================================
%             Error decay in density recovery
% ===========================================================
figure(3)
clf(3)
semilogy(Npow,error,'LineWidth',1.5,'LineStyle','-.', ...
         'Marker','*','MarkerSize',10);
grid on
xlabel('$2^{N}$','Interpreter','latex',FontSize=14)
ylabel('Abs. Error','Interpreter','latex',FontSize=14)
title_txt = sprintf('Decadimento dell''errore assoluto per t = %4.1f',t);
title(title_txt,'interpreter','latex',FontSize=12);


% ===========================================================
%                   Density time evolution
% ===========================================================
tGrid = [0.1 1 10];
xGrid = linspace(-1.5,1.5,1000);

figure(4)
clf(4)
hold on; grid on 
MARKERS = {'o','*','square'};
xlabel('$\bar{x}=\log(S_t)$', ...
       Interpreter='latex',FontSize=14)
ylabel("$\phi_{\mathtt{Hes}}(\bar{y}|\bar{x})$", ...
       Interpreter='latex',FontSize=14)
idx = 1;
legend_txt = {3};

for j = 1:length(tGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,tGrid(j),r,q, ...
                                        kappa,theta,eta,rho,v0);
    % Cumulants needed for the integration range
    c = getCumulants(model,tGrid(j),r,q,kappa,theta,eta,rho,v0);
    % Define the COS method integration range
    L = 12; 
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) );
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,xGrid,2^13,a,b);
    % Plot density time evolution
    plot(xGrid,f_XExact,'linewidth',1.5,'Marker',MARKERS(j))
    legend_txt{idx} = strcat('$$t =',num2str(tGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')