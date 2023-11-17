clear variables; close all; clc

format short

global MktIV;
global K;
global r;
global q;
global S0;
global T;
global model;

model = 'WMSVtrig';

% Load script with market data
marketData = readmatrix("SPXmarketdata.xlsx");

MktIV = marketData(2:8,2:end) ./ 100;
T = marketData(1,2:end) ./ 252;
K = marketData(2:8,1);
K = repmat(K,1,length(T));
r = 5.31/100;   % T-bill 3m USA at 10/20/2023
q = 0;
S0 = 4278;


% Initial Parameter Guess for WMSV Model
x0 = [-0.1948, 0, -1.4967, ... % matrix M
      -0.3422, 0, -0.6210, ... % matrix R
       0.0076, 0, -0.0290, ... % matrix Q
       0.0010, 1e-3, 0.0238, ... % matrix Sigma
      1.9544
     ];          

options = optimoptions('fmincon', ...
                       'Algorithm','trust-region-reflective', ...
                       'Diagnostics','on', ...
                       'MaxFunctionEvaluations',300, ...
                       'MaxIterations',500, ...
                       'FunctionTolerance',1e-10, ...
                       'StepTolerance', 1e-12, ...
	                   'Display','iter');

% Constraints (Lower and Upper bounds on Parameters)
lb = [-1.50, -1.50, -1.50, ... % matrix M
      -0.99, -0.99, -0.99, ... % matrix R
      -0.10, -0.10, -0.10, ... % matrix Q
      1e-4, 1e-4, 1e-4, ... % matrix Sigma
      1.0
     ];

ub = [0.25, 0.25, 0.25, ... % matrix M
      0.10, 0.10, 0.10, ... % matrix R
      0.95, 0.95, 0.95, ... % matrix Q
      0.35, 0.35, 0.35, ... % matrix Sigma
      3.0
     ];

% Linear inequality constraints
A = []; b = [];

% Linear constraints
Aeq = []; beq = [];

% Calibration start
tic
[x,resnorm,residual,exitflag, ...
 output,lambda,jacobian] = fmincon(@sigmaImpDiff,x0,lb,ub, ...
                                     A,b, ...
                                     Aeq,beq, ...
                                     @nonlinconstraint, ...
                                     options);
t1 = toc;

fprintf('=================================\n');
fprintf('   M_11     | %10.4f\n',x(1))
fprintf('---------------------------------\n');
fprintf('   M_12     | %10.4f\n',x(2))
fprintf('---------------------------------\n');
fprintf('   M_22     | %10.4f\n',x(3))
fprintf('---------------------------------\n');
fprintf('   R_11     | %10.4f\n',x(4))
fprintf('---------------------------------\n');
fprintf('   R_12     | %10.4f\n',x(5))
fprintf('---------------------------------\n');
fprintf('   R_22     | %10.4f\n',x(6))
fprintf('---------------------------------\n');
fprintf('   Q_11     | %10.4f\n',x(7))
fprintf('---------------------------------\n');
fprintf('   Q_12     | %10.4f\n',x(8))
fprintf('---------------------------------\n');
fprintf('   Q_22     | %10.4f\n',x(9))
fprintf('---------------------------------\n');
fprintf('   Sigma_11 | %10.4f\n',x(10))
fprintf('---------------------------------\n');
fprintf('   Sigma_12 | %10.4f\n',x(11))
fprintf('---------------------------------\n');
fprintf('   Sigma_22 | %10.4f\n',x(12))
fprintf('---------------------------------\n');
fprintf('   Gindikin | %10.4f\n',x(13))
fprintf('=================================\n');
fprintf('   EstTime  | %10.4f\n',t1)
fprintf('=================================\n');
fprintf('\n')

M11 = x(1);  M12 = x(2);  M22 = x(3);
R11 = x(4);  R12 = x(5);  R22 = x(6);
Q11 = x(7);  Q12 = x(8);  Q22 = x(9);
Sigma11 = x(9);  Sigma12 = x(10);  Sigma22 = x(12);
beta = x(13);

price_1 = zeros(size(K,1),length(T));

for j = 1:length(T)
    price_1(:,j) = getEuropeanOptionSmileByCOS(model, ...
                   [M11,M12,M22, ... % matrix M
                    R11,R12,R22, ... % matrix R
                    Q11,Q12, Q22, ... % matrix Q
                    Sigma11,Sigma12,Sigma22, ... % matrix Sigma
                    beta % Gindikin
                    ], ...
                   'Call', ...
                   S0, r, q, T(j), K(:,j), ...
                   2^13, ...
                   12);
end

fprintf('price=\n')
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

IVrmse = abs(sum(sum((ModelIV - MktIV)))) / (length(T)*size(K,1));
fprintf('IVRMSE = %10.7f\n',IVrmse)


% =======================================================================
% Implied volatility for SPX500 Index - [WMSVtrig]
% =======================================================================
X = 1:7;
fig = figure('PaperSize',[21 29.7]);

for t = 1:length(T)
	subplot(3,2,t)
	plot(X',MktIV(:,t),'kx',X',ModelIV(:,t),'ro','MarkerSize',10)
    ylim([0.05 0.35])
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',{'4278','4491.9','4705.8','4919.7', ...
                          '5133.6','5347.5','5561.4'})
	title(['SPX - Maturity ' num2str(T(t)*252) ' giorni'])
	legend('Market IV', 'WMSV IV','Orientation','horizontal')
    legend boxoff
end


% =======================================================================
% Implied volatility SPX500 surface model vs market - [WMSVtrig]
% =======================================================================
figure(4)
surf(K,repmat(T, size(K,1),1), MktIV)
hold on
scatter3(K,repmat(T, size(K,1),1), ModelIV, 'ro', 'filled')
zlim([0.05 0.35])
xlabel('Strike Price')
ylabel('Tempo (anni)')
title('SPX Implied Volatility')

% Non-linear constraints
function [c,ceq] = nonlinconstraint(x)
    c(1) = x(11) - sqrt(x(10) * x(12));
    c(2) = abs(x(2)) - sqrt(x(1) * x(3));
    c(3) = abs(x(8)) - sqrt(x(7) * x(9));
    ceq = [];
end