function value = getEuropeanOptionSmileByCOS(model,param,type,S0,r,q,t,K,N,L)
    % cf     - Characteristic function
    % param  - vector of params
    % type    - 'Call' or 'Put'
    % S0     - Initial stock price
    % r      - risk free rate
    % t      - Time to maturity
    % K      - Vector of strike prices
    % N      - Number of expansion terms
    % L      - Size of truncation domain (typ.:L=8 or L=10)
    i = complex(0,1);
    x0 = log(S0 ./ K);
    % Truncation domain
    a = 0.0 - L * sqrt(t);
    b = 0.0 + L * sqrt(t);
    k = 0:N-1; % Row vector, index for expansion terms
    u = k * pi / (b - a); % ChF cosine arguments
    H_k = getPayoffCosineCoefficients(type,a,b,k);
    cf = @(u) getCharacteristicFunction(model,u,t,r,q,param);
    temp = (cf(u) .* H_k).';
    temp(1) = 0.5 * temp(1); % Adjust the first element by 1/2
    mat = exp(i * (x0 - a) * u); % Matrix-vector manipulations
    % Final output
    value = exp(-r * t) * K .* real(mat * temp);
    for j = 1:length(K)
        if value(j) < 0
            value(j) = 1e-20;
        end
    end
end


function H_k = getPayoffCosineCoefficients(type,a,b,k)
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
end

function [Chi_k, Psi_k] = getChiPsi(a,b,c,d,k)
    % Evaluation of Psi
    Psi_k = sin(k*pi * (d - a)/(b - a)) - sin(k*pi * (c - a)/(b - a));
    Psi_k(2:end) = Psi_k(2:end) * (b - a) ./ (k(2:end) * pi); % for k != 0
    Psi_k(1) = d - c; % for k = 0
    % Evaluation of Chi
    Chi_k = 1.0 ./ (1.0 + (k*pi / (b - a)).^2);
    expr1 = cos(k*pi * (d - a)/(b - a))*exp(d) - cos(k*pi * (c - a)/(b - a))*exp(c);
    expr2 = k*pi/(b - a) .* sin(k*pi * (d - a)/(b - a))*exp(d) - k*pi/(b - a) .* sin(k*pi * (c - a)/(b - a))*exp(c);
    Chi_k = Chi_k .* (expr1 + expr2);
end
