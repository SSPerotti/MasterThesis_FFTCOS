function [price,delta,gamma] = getOptionPriceDeltaGammaByCOS(cf,c,product,type,S0,r,t,K,N,L)
% GETOPTIONPRICEDELTAGAMMABYCOS(cf,c,product,type,S0,r,t,K,N,L)
% 
% Insert detailded description here
%
% Literature: see the paper by F. Fang and C. W. Oosterlee [2008],
%             «A Novel Pricing Method for European Options Based on 
%              Fourier-Cosine Series Expansions».
%
% INPUT
% cf:         Characteristic function
% c:          Model cumulants
% product:    Type 'European' or 'Digital'
% type:       Type 'Call' or 'Put'
% S0:         Initial stock price
% r:          Risk free rate
% t:          Time to maturity
% K:          Vector of strike prices
% N:          Number of expansion terms (typ.: a power of 2)
% L:          Size of truncation domain (typ.: L=8 or L=10)
%
% OUTPUT
% price:      Option value at t_0
% delta:      Delta value at t_0
% gamma:      Gamma value at t_0
%
% -----------------------------------------------------------------------

i = complex(0,1);

% Center the grid
x0 = log(S0 ./ K);

% Truncation domain
a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4)) ) ); % a = - L * sqrt(t);
b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4)) ) ); % b = + L * sqrt(t);
% Alternative for truncation domain
% a = - L * sqrt(t);
% b = + L * sqrt(t);

% Row vector, index for expansion terms
k = 0:N-1; 

% ChF cosine arguments
u = k * pi / (b - a); 

% Hk coefficients for payoff function
H_k = getPayoffCosineCoefficients(product,type,a,b,k);

temp = (cf(u) .* H_k).';
temp(1) = 0.5 * temp(1);      % Adjust the first element by 1/2
mat = exp(i * (x0 - a) * u);  % Matrix-vector manipulations 

mat1 = exp(-i * (x0 - a) * u) .* (i * u); 
mat2 = exp(-i * (x0 - a) * u) .* ((i * u).^2 - (i * u)); 

% Final output
price = exp(-r * t) * K .* real(mat * temp);
delta = exp(-r * t) * K .* real(mat1 * temp) .* (1 ./ S0);
gamma = exp(-r * t) * K .* real(mat2 * temp) .* (1 ./ S0^2);

end


function H_k = getPayoffCosineCoefficients(product,type,a,b,k)
    if strcmp(product,'European')
        if strcmp(type,'Call')
            c = 0;
            d = b;
            [Chi_k,Psi_k] = getChiPsi(a,b,c,d,k);
            if a < b && b < 0.0
                H_k = zeros([length(k),1]);
            else
                H_k = 2.0 / (b - a) * (Chi_k - Psi_k);
            end
        elseif strcmp(type,'Put')
            c = a;
            d = 0.0;
            [Chi_k,Psi_k] = getChiPsi(a,b,c,d,k);
            H_k = 2.0 / (b - a) * (- Chi_k + Psi_k);
        end
    
    elseif strcmp(product,'Digital')
        if strcmp(type,'Call')
            c = 0;
            d = b;
            [~,Psi_k] = getChiPsi(a,b,c,d,k);
            if a < b && b < 0.0
                H_k = zeros([length(k),1]);
            else
                H_k = 2.0 / (b - a) * Psi_k;
            end
        elseif strcmp(type,'Put')
            c = a;
            d = 0.0;
            [~,Psi_k] = getChiPsi(a,b,c,d,k);
            H_k = 2.0 / (b - a) * Psi_k;
        end
    end
end


function [Chi_k, Psi_k] = getChiPsi(a,b,c,d,k)
    % Evaluation of Psi
    Psi_k = sin(k*pi * (d - a)/(b - a)) ...
            - sin(k*pi * (c - a)/(b - a));
    Psi_k(2:end) = Psi_k(2:end) * (b - a) ./ (k(2:end) * pi); % for k != 0
    Psi_k(1) = d - c; % for k = 0
    % Evaluation of Chi
    Chi_k = 1.0 ./ (1.0 + (k*pi / (b - a)).^2);
    expr1 = cos(k*pi * (d - a)/(b - a))*exp(d) ...
            - cos(k*pi * (c - a)/(b - a))*exp(c);
    expr2 = k*pi/(b - a) .* sin(k*pi * (d - a)/(b - a))*exp(d) ...
            - k*pi/(b - a) .* sin(k*pi * (c - a)/(b - a))*exp(c);
    Chi_k = Chi_k .* (expr1 + expr2);
end
