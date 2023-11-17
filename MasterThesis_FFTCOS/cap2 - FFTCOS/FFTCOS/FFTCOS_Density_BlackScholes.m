%************************************************************************%
% Here we recover the density for Black Scholes [1973] model             %
%                                                                        %
%************************************************************************%
clear variables; close all; clc

format long e

model = 'BlackScholes';

% Parameters
S0 = 100;            % Initial stock price (not used for density)
r = 0.05;            % Risk free rate 
q = 0;               % Dividend yield
t = 0.1;             % Time to maturity
sigma = 0.25;        % Volatility


% Domain for the density f(x)
% use [-1.0, 1.0] for t = 0.1 
% use [-1.0, 1.0] for t = 1.0
% use [-5.0, 5.0] for t = 10
x = linspace(-0.5,0.5,1000); 


% ===========================================================
%                   Model density recovery 
% ===========================================================
% Characteristic function 
cf = @(u) getCharacteristicFunction(model,u,t,r,q,sigma);

% Cumulants needed for the integration range
c = getCumulants(model,t,r,q,sigma);

% Define the COS method integration range
L = 12; N = 8192;
a = c(1) - L * sqrt( c(2) + sqrt(c(4)) ); 
b = c(1) + L * sqrt( c(2) + sqrt(c(4)) ); 

% Reference value (N = 2^16)
f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);

figure(1)
plot(x,f_XExact,'linewidth',1.5,'Marker','o')
x_txt = '$\bar{x}=\log(S_t)$';
y_txt = '$\phi_{\mathtt{BS}}(\bar{y}|\bar{x})$';
title_txt = sprintf(['%s Density recovery ' ...
                     'in t = %4.1f with COS'],model,t);

xlabel(x_txt,Interpreter='latex',FontSize=14)
ylabel(y_txt,Interpreter='latex',FontSize=14)
title(title_txt)
grid on


% ===========================================================
%             Effect of [N] in density recovery 
% ===========================================================
% Iterate over different numbers of expansion terms
Npow = 2:1:6;
Ngrid = 2.^Npow;
idx = 1;
error = zeros(length(Ngrid),1);

figure(2)
hold on; grid on
xlabel(x_txt,Interpreter='latex',FontSize=14)
ylabel(y_txt,Interpreter='latex',FontSize=14)
legend_txt = cell(length(Ngrid),1);
time_txt = cell(length(Ngrid),1);

for N = Ngrid
    % Density from the COS method
    tic
    for j = 1:1 % 10^3 experiments
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
time_txt
%{
title_txt = sprintf(['Effect of N in %s density ' ...
                     'recovery in t = %4.1f with COS'],model,t);
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
xlabel('$\bar{x}=\log(S_t)$',Interpreter='latex',FontSize=14)
ylabel("$\phi_{\mathtt{BS}}(\bar{y}|\bar{x})$",Interpreter='latex',FontSize=14)
idx = 1;
legend_txt = {3};

for j = 1:length(tGrid)
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,tGrid(j),r,q,sigma);
    % Cumulants needed for the integration range
    c = getCumulants(model,tGrid(j),r,q,sigma);
    % Define the COS method integration range
    L = 12; N = 8192;
    a = c(1) - L * sqrt( c(2) + sqrt(c(4)) );
    b = c(1) + L * sqrt( c(2) + sqrt(c(4)) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,xGrid,N,a,b);
    % Plot density time evolution
    plot(xGrid,f_XExact,'linewidth',1.5,'Marker',MARKERS(j))
    legend_txt{idx} = strcat('$$t =',num2str(tGrid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')