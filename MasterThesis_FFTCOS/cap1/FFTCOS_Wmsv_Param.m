%************************************************************************%
% Here we study the effects of the parameter vector                      %
% that parametrizes WMSV model                                           %
%                                                                        %
%************************************************************************%
clear variables; close all; clc

model = 'WMSV';

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

x_txt = '$\bar{x}=\log(S_t)$';
y_txt = '$\phi_{\mathtt{WMSV}}(\bar{y}|\bar{x})$';

% COS method settings
N = 4096;
L = 8;


% ===========================================================
% Effect of [M_11] on WMSV density
% ===========================================================
m11grid = [-3.0 -1.5 -0.5 -0.05]; 
time = 1.0;
x = linspace(-1.0,1.0,1000);

figure(5)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(m11grid)
    m11Temp = m11grid(j);
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        m11Temp, m12, m21, m22, ...
                                        r11, r12, r21, r22, ...
                                        q11, q12, q21, q22, ...
                                        Sig11, Sig12, Sig21, Sig22, ...
                                        beta);
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q, ...
                     m11Temp, m12, m21, m22, ...
                     r11, r12, r21, r22, ...
                     q11, q12, q21, q22, ...
                     Sig11, Sig12, Sig21, Sig22, ...
                     beta);
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$M_{11}=',num2str(m11grid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of M11 on WMSV density')


% ===========================================================
% Effect of [M_22] on WMSV density
% ===========================================================
m22grid = [-3.0 -1.5 -0.5 -0.05]; 
time = 1.0;
x = linspace(-1.0,1.0,1000);

figure(6)
hold on; grid on 
MARKERS = {'o','*','square','diamond'};
LINESTYLE = {'-','-.','-.','-'};
xlabel(x_txt,Interpreter='latex',FontSize=12)
ylabel(y_txt,Interpreter='latex',FontSize=12)
idx = 1;
legend_txt = {4};

for j = 1:length(m22grid)
    m22Temp = m22grid(j);
    % Characteristic function
    cf = @(u) getCharacteristicFunction(model,u,time,r,q, ...
                                        m11, m12, m21, m22Temp, ...
                                        r11, r12, r21, r22, ...
                                        q11, q12, q21, q22, ...
                                        Sig11, Sig12, Sig21, Sig22, ...
                                        beta);
    % Cumulants needed for the integration range
    c = getCumulants(model,time,r,q, ...
                     m11, m12, m21, m22Temp, ...
                     r11, r12, r21, r22, ...
                     q11, q12, q21, q22, ...
                     Sig11, Sig12, Sig21, Sig22, ...
                     beta);
    % Define the COS method integration range
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4))) ); 
    % Exact density
    f_XExact = getRecoveredDensityByCOS(cf,x,N,a,b);
    % Plot density time evolution
    plot(x,f_XExact,'linewidth',1.5,'LineStyle',LINESTYLE(j), ...
         'Marker',MARKERS(j),'MarkerIndices',1:20:length(x))
    legend_txt{idx} = strcat('$$M_{22}=',num2str(m22grid(j)),'$$');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')
title('Effect of M22 on WMSV density')