clear variables; close all; clc

% Load script with market data
marketDataDAX = readmatrix("DAX30marketdata.xlsx");

MktIV_DAX = marketDataDAX(2:8,2:end) ./ 100;
T_DAX = marketDataDAX(1,2:end) ./ 252;
K_DAX = marketDataDAX(2:8,1);
K_DAX = repmat(K_DAX,1,length(T_DAX));
r = 5.31/100;   % T-bill 3m USA at 10/20/2023
q = 0;
S0_DAX = 14889.46;

figure(1)
clf(1)
surf(K_DAX,repmat(T_DAX, size(K_DAX,1),1), MktIV_DAX)
zlim([0.10 0.25])
xlabel('Strike Price')
ylabel('Tempo (anni)')
title('DAX30 Implied Volatility')


% Load script with market data
marketDataSPX = readmatrix("SPXmarketdata.xlsx");

MktIV_SPX = marketDataSPX(2:8,2:end) ./ 100;
T_SPX = marketDataSPX(1,2:end) ./ 252;
K_SPX = marketDataSPX(2:8,1);
K_SPX = repmat(K_SPX,1,length(T_SPX));
r = 5.31/100;   % T-bill 3m USA at 10/20/2023
q = 0;
S0_SPX = 4278;

figure(2)
clf(2)
surf(K_SPX,repmat(T_SPX, size(K_SPX,1),1), MktIV_SPX)
zlim([0.10 0.25])
xlabel('Strike Price')
ylabel('Tempo (anni)')
title('SPX Implied Volatility')
