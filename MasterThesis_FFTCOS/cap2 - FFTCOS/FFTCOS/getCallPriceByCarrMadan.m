function resu = getCallPriceByCarrMadan(S0,K,r,T,cf)
% GETCALLPRICEBYCARRMADAN(S0,K,r,T,cf) 
%
% Here we implement a basic function to calculate the price of a 
% European call using the Carr-Madan method
%
% Literature: TODO
%
% Reference: Lecture Notes on Computational Finance, Alessandro Gnoatto,
%            summer 2015
%
% INPUT
% S0:        Initial stock price
% K:         Strike price
% r:         Risk free rate
% T:         Time to maturity
% cf:        Characteristic function
%
% OUTPUT
% resu:      Call option value at t_0
%
% ------------------------------------------------------------------------

i = complex(0,1);

% Dampening factor
alpha = 1.5;

N = 8196;   % Fourier grid (a power of 2)
eta = 0.25;
lambda = (2*pi) / (N*eta);  % Equation 23 in Carr-Madan
u = 0 : eta : ((N-1)*eta);   
b = (N * lambda)/2;         % Equation 20 in Carr-Madan

ncf = exp(-r*T) .* (cf(u-(alpha+1)*i) .* exp(i*(u-(alpha+1)*i)*log(S0))); 
dcf = (alpha^2 + alpha - u.^2 + i * (2*alpha + 1) .* u);
rcf = ncf ./ dcf;

% Summation term in equation 22 in Carr-Madan
tmp = rcf .* exp(i*u*b) * eta;

% Applying Simpson's rule - formula 24 in Carr-Madan
jvec = 1:N;
tmp = (tmp/3) .* (3 + (-1) .^ jvec - ((jvec - 1)==0));

ft = fft(tmp,N);
ft = ft';
% Log strikes ordinates as in equation 19
kv = -b:lambda:((N-1)*lambda-b);
kv = kv';
% Call price vector resulting in equation 22
cpvec = (exp(-alpha*kv).*ft)/pi; 

% Perform interpolation so to give the price
kvec = exp(kv);
resu = zeros(length(K),1);

%{
for j = 1:length(K)
    [i,value] = find(kvec>K(j),1,'first');
    xp = [kvec(i-1) kvec(i)];
    yp = [real(cpvec(i-1)) real(cpvec(i))];
    resu(j) = interp1(xp,yp,K(j));
end
%}

for j = 1:length(K)
    for i = 1:N
        if kvec(i) > K(j)
            break;
        end
    end
    i = find(kvec>K(j),1,'first');
    xp = [kvec(i-1) kvec(i)];
    yp = [real(cpvec(i-1)) real(cpvec(i))];
    resu(j) = max(real(interp1(xp,yp,K(j))), 0);
end

end
