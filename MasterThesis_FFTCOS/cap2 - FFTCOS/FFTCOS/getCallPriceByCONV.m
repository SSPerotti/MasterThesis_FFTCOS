function resu = getCallPriceByCONV(intBound,alpha,S0,K,T,rf,cf)
% GETCALLPRICEBYCONV(intBound,alpha,S0,K,T,rf,cf)
% 
% Implement a basic function to calculate the price of an European call 
% option using the CONV method
%
% Literature: TODO
%
% Reference: Lecture Notes on Computational Finance, Alessandro Gnoatto,
%            summer 2015
%
% INPUT
% inteBound:     TODO
% alpha:         TODO
% S0:            Initial stock price
% K:             Strike Price
% T:             Time to maturity
% rf:            Risk free rate
% cf:            Characteristic function
%
% OUTPUT
% resu:      Call option value at t_0
%
% ------------------------------------------------------------------------

i = complex(0,1);

% N gridpoints for the FFT we choose a power of 2 as in Carr-Madan
numberOfPoints = 4096; % 8192 or 4096

% Risk free time maturity
rdt = rf * T;

% Grid spacing
Delta_y = intBound / numberOfPoints;
Delta_x = Delta_y;
Delta_u = 2 * pi / intBound;

% Grids
gridIndex = (0:numberOfPoints-1)';
minusOnePowerP = (-1) .^ gridIndex;

% Grid for the log initial stock price
x = (gridIndex .* Delta_x) - (numberOfPoints/2 * Delta_x);

% Grid for the log terminal stock price (the initial factor is just the
% rescale factor \epsilon_y from eq (44) in Lord et al. [2008])
y = log(K / S0) + (gridIndex .* Delta_y) ...
    - (numberOfPoints/2 * Delta_y);

% Grid for the Fourier frequency
u = (gridIndex .* Delta_u) - (numberOfPoints/2 * Delta_u);

% Payoff
payoff = max((S0*exp(y) - K), 0);

% Dampened option value
dampenedPayoff = payoff .* exp(alpha .* y);

% Weights for the trapezoidal rule
integrationWeights = ones(numberOfPoints,1);
integrationWeights(1) = 0.5;
integrationWeights(numberOfPoints) = 0.5;

transformedDampenedPayoff = ifft( ...
                                  (minusOnePowerP) .* integrationWeights ...
                                  .* dampenedPayoff ...
                                );

integrand = exp(i .* gridIndex .* (y(1) - x(1)) .* Delta_u) ...
            .* cf(-(u - (i * alpha))) ...
            .* transformedDampenedPayoff;

% Call prices for different values of the log initial stock price
callPrices = abs( ...
                 exp( -rdt - (alpha .* x) + (i .* u .* (y(1) - x(1))) ...
                )...
             .* (minusOnePowerP) .* fft(integrand));

% Final value: the initial stock is in the middle of the grid
resu = callPrices(numberOfPoints/2 + 1, 1); 

end

