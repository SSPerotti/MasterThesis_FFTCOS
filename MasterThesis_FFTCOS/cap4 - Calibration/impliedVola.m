function resu = impliedVola(S0,K,r,T,mkt_price,q)
% IMPLIEDVOLA(S0,K,r,T,mkt_price,q)
% It provides a bootstrapper for implied volatility. 
% 
% Implied volatility is the value for sigma in the Black-Scholes formula 
% for which the Black-Scholes formula returns the market price mkt_price 
%
% INPUT
% S:            Initial stock price
% K:            Strike price
% r:            Risk free rate (in FX market domestic risk free rate r_d)
% T:            Maturity
% mkt_price:    Market price
% q:            Dividend yield (in FX market foreign risk free rate r_f)
%
% OUTPUT
% resu:         Implied volatility
%
% ------------------------------------------------------------------------

    options = optimset('fzero');
    options = optimset(options,'TolX',1e-10,'Display','off');
    
    [x,~,exitflag] = fzero(@difference, 0.5, options, ...
                           S0, K, r, T, mkt_price, q);
    
    if(exitflag ~= 1)
        disp('Did not find value')
        resu = 0;
    else
        resu = x;
    end

end

function resu = difference(sigma,S,K,r,T,mkt_price,q)
    bs_price = blackScholes(sigma,S,K,r,T,q);
    resu = mkt_price - bs_price;
end

function resu = blackScholes(sigma,S,K,r,T,q)
    d1 = (log(S./K) + (r - q + 0.5*sigma^2).*T) ./ (sigma.*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    resu = S*exp(-r*T)*normal(d1,0,1) - K*exp(-r*T)*normal(d2,0,1);
end

function resu = normal(x,mu,sigma)
    z = (x - mu) ./ sigma;
    resu = 0.5 * erfc(-z ./ sqrt(2));
end







