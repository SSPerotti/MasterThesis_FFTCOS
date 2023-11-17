function [Lt]=computeLtW(x,Sigma0,w,A_i,A_j,M,Q,Gindikin,R,t,T,r_i,r_j)
% This function computes the Laplace Transform of the vector process.

% fix a dimension index
d = size(Sigma0,1);

exponent = (T-t).*[M + w.*(Q')*R*(A_i-A_j), -2*(Q')*Q; 
                0.5.*(w^2-w).*(A_i-A_j)^2, -((M') + w.*(A_i-A_j)*(R')*Q)];

BB = expm(exponent);

BB_22 = BB(d+1:2*d,d+1:2*d);
BB_21 = BB(d+1:2*d,1:d);


cal_A = w.*(r_i-r_j).*(T-t) ...
        -0.5*Gindikin*trace(logm(BB_22) + (M'+w.*(A_i-A_j)*(R')*Q)*(T-t));

cal_b = (BB_22^(-1))*BB_21;

Lt = exp(w.*x + cal_A + trace(cal_b*Sigma0));
