function value = getOptionPriceByCOS(cf,c,product,type,S0,r,t,K,N,L)
    % cf - Characteristic function
    % c  - Model cumulants
    % CP - C for call and P for put
    % S0 - Initial stock price
    % r  - risk free rate
    % t  - Time to maturity
    % K  - Vector of strike prices
    % N  - Number of expansion terms
    % L  - Size of truncation domain (typ.:L=8 or L=10)
    i = complex(0,1);
    x0 = log(S0 ./ K);
    % Truncation domain
    a = c(1) - L * sqrt( abs(c(2)) + sqrt( abs(c(4)) ) ); % a = 0.0 - L * sqrt(t);
    b = c(1) + L * sqrt( abs(c(2)) + sqrt( abs(c(4)) ) ); % b = 0.0 + L * sqrt(t);
    k = 0:N-1; % Row vector, index for expansion terms
    u = k * pi / (b - a); % ChF cosine arguments
    H_k = getPayoffCosineCoefficients(product,type,a,b,k);
    temp = (cf(u) .* H_k).';
    temp(1) = 0.5 * temp(1); % Adjust the first element by 1/2
    mat = exp(i * (x0 - a) * u); % Matrix-vector manipulations
    % Final output
    value = exp(-r * t) * K .* real(mat * temp);
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

end
