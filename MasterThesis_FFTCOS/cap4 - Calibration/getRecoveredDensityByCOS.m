function resu = getRecoveredDensityByCOS(cf,x,N,a,b)
% GETRECOVEREDDENSITYBYCOS(cf,x,N,a,b)
%
% Compute the cosine series expansion of the probability density function.
% 
% Literature: see the paper by F. Fang C. W. Oosterlee [2008] "A Novel 
%             Pricing Method for European Options based on Fourier-Cosine Series 
%             Expansions".
% 
% INPUTS 
% cf :     Characteristic function (given as a function)
% y :      Evaluation point
% N :      Number of grid points
% a :      Lower bound for inverse Fourier integral  
% b :      Upper bound for inverse Fourier integral 
%
% OUTPUT
% resu :   Probability density function approximation 
%
% -----------------------------------------------------------------------

i = complex(0,1); 
k = 0:N-1;
% Cosine argument
u = k * pi / (b - a); 
% Compute A_{k} cosine series coefficients approximation
A_k = 2 / (b - a) * real(cf(u) .* exp(-i * u * a)); 
% Correction for the first term
A_k(1) = A_k(1) * 0.5; 
% Final calculation
resu = A_k * cos(u' * (x - a)); 
end