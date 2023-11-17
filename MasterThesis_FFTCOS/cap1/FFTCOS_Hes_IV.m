clear variables; close all; clc

model = 'Heston';

% Option parameters
T = 1;
S0 = 100;
product = 'European';
type = 'Call';

% Heston model parameters
kappa = 0.1;
theta = 0.1;
eta = 0.1;
rho = -0.75;
v0 = 0.05;
r = 0.05;
q = 0.0;

% Range of strike prices
K = linspace(40,200,25)';

% COS method settings
N = 4096;
L = 8;


% ============================================================
% Effect of kappa on Heston Implied Volatility
% ============================================================
kappaVec = [0.1, 0.5, 1.0, 2.0];
ivM = zeros(length(K),length(kappaVec));
argLegend = cell(5,1);
idx = 1;
for i=1:length(kappaVec)
    kappaTemp = kappaVec(i);
    % Compute ChF for the Heston model
    cf = @(u) getCharacteristicFunction(model,u,T,r,q, ...
                                        kappaTemp,theta,eta,rho,v0);
    % The COS method
    optPrice = getOptionPriceByCOS(cf,product,type,S0,r,T,K,N,L);
    for j = 1:length(K)
        strike = K(j);
        call = optPrice(j);
        ivM(j,i) = ImpliedVolatility(type,call,strike,T,S0,r,0.3)*100;
    end
    argLegend{idx} = sprintf('\\kappa = %.1f',kappaTemp);
    idx = idx + 1;
end

MakeFigure(K,ivM,argLegend,'Effect of \kappa on implied volatility')


% ============================================================
% Effect of theta on Heston Implied Volatility
% ============================================================
thetaVec = [0.03, 0.1, 0.2, 0.3];
ivM = zeros(length(K),length(thetaVec));
argLegend = cell(5,1);
idx = 1;

for i=1:length(thetaVec)
    thetaTemp = thetaVec(i);
    % Compute ChF for the Heston model
    cf = @(u) getCharacteristicFunction(model,u,T,r,q, ...
                                        kappa,thetaTemp,eta,rho,v0);
    % The COS method
    optPrice = getOptionPriceByCOS(cf,product,type,S0,r,T,K,N,L);
    for j = 1:length(K)
        strike = K(j);
        call = optPrice(j);
        ivM(j,i) = ImpliedVolatility(type,call,strike,T,S0,r,0.3)*100;
    end
    argLegend{idx} = sprintf('\\theta = %.2f',thetaTemp);
    idx = idx + 1;
end

MakeFigure(K,ivM,argLegend,'Effect of \theta on implied volatility')


% ============================================================
% Effect of eta (vol-of-vol) on Heston Implied Volatility
% ============================================================
etaVec = [0.1, 0.3, 0.5,0.9];
ivM = zeros(length(K),length(etaVec));
argLegend = cell(5,1);
idx = 1;
for i=1:length(etaVec)
    etaTemp = etaVec(i);
    % Compute ChF for the Heston model
    cf = @(u) getCharacteristicFunction(model,u,T,r,q, ...
                                        kappa,theta,etaTemp,rho,v0);
    % The COS method
    optPrice = getOptionPriceByCOS(cf,product,type,S0,r,T,K,N,L);
    for j = 1:length(K)
        strike = K(j);
        call = optPrice(j);
        ivM(j,i) = ImpliedVolatility(type,call,strike,T,S0,r,0.3)*100;
    end
    argLegend{idx} = sprintf('\\eta = %.2f',etaTemp);
    idx = idx + 1;
end

MakeFigure(K,ivM,argLegend,'Effect of \eta on implied volatility')
% MakeFigure(K,ivM,argLegend)


% ============================================================
% Effect of rho on Heston Implied Volatility
% ============================================================
rhoVec = [0.30, -0.30, -0.50, -0.70];
ivM = zeros(length(K),length(rhoVec));
argLegend = cell(5,1);
idx = 1;

for i=1:length(rhoVec)
    rhoTemp = rhoVec(i);
    % Compute ChF for the Heston model
    cf = @(u) getCharacteristicFunction(model,u,T,r,q, ...
                                        kappa,theta,eta,rhoTemp,v0);
    % The COS method
    optPrice = getOptionPriceByCOS(cf,product,type,S0,r,T,K,N,L);
    for j = 1:length(K)
        strike = K(j);
        call = optPrice(j);
        ivM(j,i) = ImpliedVolatility(type,call,strike,T,S0,r,0.3)*100;
    end
    argLegend{idx} = sprintf('\\rho = %.2f',rhoTemp);
    idx = idx + 1;
end

MakeFigure(K,ivM,argLegend,'Effect of \rho on implied volatility')
% MakeFigure(K,ivM,argLegend)


%************************************************************************%
% Closed-form expression of European call/put option                     %
% with Black-Scholes formula                                             %
%************************************************************************%
function resu = getEuropeanOptionPriceByBS(type,S0,K,sigma,T,r)
    d1 = (log(S0 ./ K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    if strcmp(type,'Call') 
        resu = normal(d1,0,1) * S0 - normal(d2,0,1) .* K * exp(-r * T);
    elseif strcmp(type,'Put')
        resu = normal(-d2,0,1) .* K*exp(-r*T) - normal(-d1,0,1)*S0;
    end
end


function impliedVol = ImpliedVolatility(type,marketPrice,K,T,S0,r,initialVol)
    func = @(sigma) (getEuropeanOptionPriceByBS(type,S0,K,sigma,T,r) ...
                     - marketPrice).^1.0;
    impliedVol = fzero(func,initialVol);
end


function resu = normal(x,mu,sigma)
    z = (x - mu) ./ sigma;
    resu = 0.5 * erfc(-z ./ sqrt(2));
end


function MakeFigure(X1, YMatrix1, argLegend, titleIn)
% function MakeFigure(X1, YMatrix1, argLegend)
    %CREATEFIGURE(X1,YMATRIX1)
    % X1: vector of x data
    % YMATRIX1: matrix of y data
    % Auto-generated by MATLAB on 16-Jan-2012 15:26:40
    % Create figure
    figure1 = figure('InvertHardcopy','off',...
    'Colormap',[0.061875 0.061875 0.061875; 0.06875 0.06875 0.06875; ...
                0.075625 0.075625 0.075625; 0.0825 0.0825 0.0825; ...
                0.089375 0.089375 0.089375; 0.09625 0.09625 0.09625; ...
                0.103125 0.103125 0.103125; 0.11 0.11 0.11; ...
                0.146875 0.146875 0.146875; 0.18375 0.18375 0.18375; ... 
                0.220625 0.220625 0.220625; 0.2575 0.2575 0.2575; ...
                0.294375 0.294375 0.294375; 0.33125 0.33125 0.33125; ...
                0.368125 0.368125 0.368125; 0.405 0.405 0.405; ...
                0.441875 0.441875 0.441875; 0.47875 0.47875 0.47875; ...
                0.515625 0.515625 0.515625; 0.5525 0.5525 0.5525; ...
                0.589375 0.589375 0.589375; 0.62625 0.62625 0.62625; ...
                0.663125 0.663125 0.663125; 0.7 0.7 0.7; ... 
                0.711875 0.711875 0.711875; 0.72375 0.72375 0.72375; ...
                0.735625 0.735625 0.735625; 0.7475 0.7475 0.7475; ... 
                0.759375 0.759375 0.759375; 0.77125 0.77125 0.77125; ...
                0.783125 0.783125 0.783125; 0.795 0.795 0.795; ...
                0.806875 0.806875 0.806875; 0.81875 0.81875 0.81875; ...
                0.830625 0.830625 0.830625; 0.8425 0.8425 0.8425; ...
                0.854375 0.854375 0.854375; 0.86625 0.86625 0.86625; ...
                0.878125 0.878125 0.878125; 0.89 0.89 0.89; ...
                0.853125 0.853125 0.853125; 0.81625 0.81625 0.81625; ...
                0.779375 0.779375 0.779375; 0.7425 0.7425 0.7425; ...
                0.705625 0.705625 0.705625; 0.66875 0.66875 0.66875; ...
                0.631875 0.631875 0.631875; 0.595 0.595 0.595; ...
                0.558125 0.558125 0.558125; 0.52125 0.52125 0.52125; ...
                0.484375 0.484375 0.484375; 0.4475 0.4475 0.4475; ...
                0.410625 0.410625 0.410625; 0.37375 0.37375 0.37375; ...
                0.336875 0.336875 0.336875; 0.3 0.3 0.3; ...
                0.28125 0.28125 0.28125; 0.2625 0.2625 0.2625; ...
                0.24375 0.24375 0.24375; 0.225 0.225 0.225; ...
                0.20625 0.20625 0.20625; 0.1875 0.1875 0.1875; ...
                0.16875 0.16875 0.16875; 0.15 0.15 0.15],...
    'Color',[1 1 1]);
    % Create axes
    %axes1 = axes('Parent',figure1,'Color',[1 1 1]);
    axes1 = axes('Parent',figure1);
    grid on
    % Uncomment the following line to preserve the X-limits of the axes
    % xlim(axes1,[45 160]);
    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim(axes1,[19 26]);
    % Uncomment the following line to preserve the Z-limits of the axes
    % zlim(axes1,[-1 1]);
    box(axes1,'on');
    hold(axes1,'all');
    % Create multiple lines using matrix input to plot
    % plot1 = plot(X1,YMatrix1,'Parent',axes1,'MarkerEdgeColor',[0 0 0],...
    % 'LineWidth',1,...
    % 'Color',[0 0 0]);
    plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',1.5);
    set(plot1(1),'Marker','o','DisplayName',argLegend{1});
    set(plot1(2),'Marker','x','LineStyle','-.','DisplayName',argLegend{2});
    set(plot1(3),'Marker','square','LineStyle','-.','DisplayName',argLegend{3});
    set(plot1(4),'Marker','diamond','DisplayName',argLegend{4});
    % Create xlabel
    xlabel('strike, $K$','Interpreter','latex',FontSize=12);
    % Create ylabel
    ylabel('Implied volatility $\sigma(T,K)$','Interpreter','latex',FontSize=12);
    % Create title
    title(titleIn);
    % Create legend
    legend1 = legend(axes1,'show');
    set(legend1,'Color',[1 1 1]);
end