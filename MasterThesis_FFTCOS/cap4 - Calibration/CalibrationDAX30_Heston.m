clear variables; close all; clc

format short

global MktIV;
global K;
global r;
global q;
global S0;
global T;
global model;

model = 'Heston';

% Load script with market data
marketData = readmatrix("DAX30marketdata.xlsx");

MktIV = marketData(2:8,2:end) ./ 100;
T = marketData(1,2:end) ./ 252;
K = marketData(2:8,1);
K = repmat(K,1,length(T));
r = 5.31/100;   % T-bill 3m USA observed on 20 october 2023
q = 0;
S0 = 14889.46;  % DAX30 spot price observed on 20 october 2023


% Initial Parameter Guess for Heston Model
x0 = [5.6325, ...  % Kappa (speed of mean reversion) 
      0.0068, ...  % Theta (variance at infinity)
      0.4305, ...  % Eta (vol-of-vol)
      -0.3782, ... % Rho (correlation)
      0.0405 ...   % v0 (initial variance)
     ];

options = optimoptions('lsqnonlin', ...
                       'Algorithm','levenberg-marquardt', ...
                       'Diagnostics','on', ...
                       'MaxFunctionEvaluations',500, ...
                       'MaxIterations',500, ...
                       'FunctionTolerance',1e-10, ...
                       'StepTolerance',1e-12, ...
	                   'Display','iter');

% Bound constraints (Lower and Upper bounds on Parameters)
lb = [1e-3, 1e-3, 1e-3, -0.99, 1e-3];
ub = [10, 1.0, 0.45, 0.99, 1.0];

% Linear inequality constraints
A = []; b = [];

% Linear constraints
Aeq = []; beq = [];

% Calibration start
tic
[x,resnorm,residual,exitflag, ...
 output,lambda,jacobian] = lsqnonlin(@sigmaImpDiff,x0,lb,ub,options);
t1 = toc;

fprintf('----------------------------------------------------------------\n');
fprintf('  kappa     theta       eta        rho         v0      EstTime\n');
fprintf('----------------------------------------------------------------\n');
fprintf('%7.4f %10.4f %10.4f %10.4f %10.4f %10.3f\n',x,t1);

fprintf('\n')

kappa = x(1);
theta = x(2);
eta = x(3);
rho = x(4);
v0 = x(5);

price_1 = zeros(size(K,1),length(T));

for j = 1:length(T)
    price_1(:,j) = getEuropeanOptionSmileByCOS('Heston', ...
                       [kappa,theta,eta,rho,v0], ...
                       'Call', ...
                       S0, r, q, T(j), K(:,j), ...
                       2^13, ...
                       12);
end

fprintf('price =\n')
disp(price_1)

ModelIV = zeros(size(K,1),length(T));

for j = 1 : length(T)
    for i = 1 : size(K,1)
        ModelIV(i,j) = impliedVola(S0, K(i,j), r, ...
                                   T(j), price_1(i,j), q);
        if(isnan(ModelIV(i,j)) == 1)
            ModelIV = 1e-2;
        end
    end
end

IVmse = abs(sum(sum((ModelIV - MktIV)))) / (length(T)*size(K,1));
fprintf('IVmse = %10.7f\n',IVmse)


% =======================================================================
% Implied volatility for DAX30 Index - Heston
% =======================================================================
X = 1:7;
fig = figure('PaperSize',[21 29.7]);

for t = 1:length(T)
	subplot(3,2,t)
	plot(X',MktIV(:,t),'kx',X',ModelIV(:,t),'ro','MarkerSize',10)
    ylim([0.10 0.25])
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',{'148','156','163','170', ...
                          '178','185','193'})
	title(['DAX30 - Maturity ' num2str(T(t)*252) ' giorni'])
	legend('Market IV', 'Heston IV','Orientation','horizontal')
    legend boxoff
end


% =======================================================================
% Implied volatility DAX30 surface model vs market - Heston
% =======================================================================
figure(2)
surf(K,repmat(T, size(K,1),1), ModelIV) % MKtIV
hold on
scatter3(K,repmat(T, size(K,1),1), MktIV, 'ro', 'filled') % ModelIV
zlim([0.10 0.25])
% shading interp;
xlabel('Strike Price')
ylabel('Tempo (anni)')
zlabel('Implied Volatility [%]')
title('Heston DAX30 Implied Volatility')


% =======================================================================
% DAX30 implied risk neutral distribution - Heston
% =======================================================================
i = complex(0,1);
S0 = S0/100;  % Scale spot price 
y = linspace(120,200,1500);
fXhes = zeros(length(T),length(y));

% COS settings
N = 4096;
L = 10;

figure(3)
xlabel('DAX30 (EUR)'); ylabel('PdF')
hold on; grid on
idx = 1;
legend_txt = {6};

for j = 1:length(T)
    a = 0.0 - L * sqrt(T(j));
    b = 0.0 + L * sqrt(T(j));
    cf = @(u) getCharacteristicFunction(model,u,T(j),r,q, ...
                                        [kappa,theta,eta,rho,v0]) ...
                                        .* exp(i*u*log(S0));
    fXhes(j,:) = (1./y) .* getRecoveredDensityByCOS(cf,log(y),N,a,b);
    plot(y,fXhes(j,:),'LineWidth',1.5)
    legend_txt{idx} = strcat('$$t =',num2str(T(j)*252),'$$ giorni');
    idx = idx + 1;
end

legendObj = legend(legend_txt);
set(legendObj,'interpreter','latex')

% Non-linear constraints
function [c,ceq] = nonlinconstraint(x)
    c(1) = 2 * x(1) * x(2) - x(3) * x(3); % Feller condition
    ceq = [];
end