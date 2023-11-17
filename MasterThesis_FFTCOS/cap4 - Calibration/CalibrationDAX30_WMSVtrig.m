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
marketData = readmatrix("DAX30marketdata.xlsx");

MktIV = marketData(2:8,2:end) ./ 100;
T = marketData(1,2:end) ./ 252;
K = marketData(2:8,1);
K = repmat(K,1,length(T));
r = 5.31/100;   % T-bill 3m USA at 10/20/2023
q = 0;
S0 = 14889.46;


% Initial Parameter Guess for WMSVtrig Model
x0 = [-0.8468, 0, -0.3051, ... % matrix M
      -0.5082, 0, -0.5874, ... % matrix R
       0.0918, 0, 0.1196, ... % matrix Q
       0.0310, 1e-3, 0.0010, ... % matrix Sigma
       1.0499 ... % Gindikin
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

% Linear equality constraints
Aeq = []; beq = [];

% Calibration start
tic
[x,fval,residual,exitflag] = fmincon(@sigmaImpDiff,x0,lb,ub, ...
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
                   [x(1),x(2),x(3), ... % matrix M
                    x(4),x(5),x(6), ... % matrix R
                    x(7),x(8),x(9), ... % matrix Q
                    x(10),x(11),x(12), ... % matrix Sigma
                    x(13) % Gindikin
                    ], ...
                   'Call', ...
                   S0, r, q, T(j), K(:,j), ...
                   2^13, ... % COS No. series terms
                   12); % COS tolerance
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
% Implied volatility for DAX30 Index - [WMSVtrig]
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
	legend('Market IV', 'WMSV IV','Orientation','horizontal')
    legend boxoff
end


% =======================================================================
% Implied volatility DAX30 surface model vs market - [WMSVtrig]
% =======================================================================
figure(4)
surf(K,repmat(T, size(K,1),1), MktIV)
hold on
scatter3(K,repmat(T, size(K,1),1), ModelIV, 'ro', 'filled')
zlim([0.10 0.25])
xlabel('Strike Price')
ylabel('Tempo (anni)')
title('DAX30 Implied Volatility')

% Non-linear constraints
function [c,ceq] = nonlinconstraint(x)
    c(1) = x(11) - sqrt(x(10) * x(12));
    c(2) = abs(x(2)) - sqrt(x(1) * x(3));
    c(3) = abs(x(8)) - sqrt(x(7) * x(9));
    ceq = [];
end