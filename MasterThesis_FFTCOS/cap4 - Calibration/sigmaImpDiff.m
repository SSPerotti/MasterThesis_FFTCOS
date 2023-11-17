function resu = sigmaImpDiff(x)

global MktIV;
global K;
global r;
global q;
global S0;
global T;
global model;

price_1 = zeros(size(K,1),length(T));

for j = 1 : length(T)
    if strcmp(model,'Heston')
        price_1(:,j) = getEuropeanOptionSmileByCOS(model, ...
                       [x(1),x(2),x(3),x(4),x(5)], ... % Heston params
                       'Call', ... 
                       S0, r, q, T(j), K(:,j), ...
                       2^13, ... % COS No. series terms
                       12); % COS tolerance

    elseif strcmp(model,'Bates')
        price_1(:,j) = getEuropeanOptionSmileByCOS(model, ...
                       [x(1),x(2),x(3),x(4), ... % Bates params
                        x(5),x(6),x(7),x(8) ...
                       ], ... 
                       'Call', ...
                       S0, r, q, T(j), K(:,j), ...
                       2^13, ... % COS No. series terms
                       12); % COS tolerance
    elseif strcmp(model,'WMSV')
        price_1(:,j) = getEuropeanOptionSmileByCOS(model, ...
                       [x(1),x(2),x(3),x(4), ... % matrix M
                        x(5),x(6),x(7),x(8), ... % matrix R
                        x(9),x(10),x(11),x(12), ... % matrix Q
                        x(13),x(14),x(15),x(16), ... % matrix Sigma
                        x(17) % Gindikin
                        ], ...
                       'Call', ...
                       S0, r, q, T(j), K(:,j), ...
                       2^13, ... % COS No. series terms
                       12); % COS tolerance
    elseif strcmp(model,'WMSVdiag')
        price_1(:,j) = getEuropeanOptionSmileByCOS(model, ...
                       [x(1),x(2), ... % matrix M
                        x(4),x(3), ... % matrix R
                        x(5),x(6), ... % matrix Q
                        x(7),x(8), ... % matrix Sigma
                        x(9) % Gindikin
                        ], ...
                       'Call', ...
                       S0, r, q, T(j), K(:,j), ...
                       2^13, ... % COS No. series terms
                       12); % COS tolerance
    elseif strcmp(model,'WMSVtrig')
        price_1(:,j) = getEuropeanOptionSmileByCOS(model, ...
                       [x(1),x(2),x(3), ... % matrix M
                        x(4),x(5),x(6), ... % matrix R
                        x(7),x(8), ... % matrix Q
                        x(9),x(10),x(11),x(12), ... % matrix Sigma
                        x(13) % Gindikin
                        ], ...
                       'Call', ...
                       S0, r, q, T(j), K(:,j), ...
                       2^13, ... % COS No. series terms
                       12); % COS tolerance
    end
end


ModelIV = zeros(size(K,1),length(T));

for j = 1 : length(T)
    for i = 1 : size(K,1)
        ModelIV(i,j) = impliedVola(S0, K(i,j), r, ...
                                   T(j), price_1(i,j), q);
    end
end

% resu = ModelIV(:) - MktIV(:);
resu = (ModelIV(:) - MktIV(:)) ./ MktIV(:) ;

end

