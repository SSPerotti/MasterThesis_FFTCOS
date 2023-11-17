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
