function resu = getCharacteristicFunction(model,u,t,r,q,param)

ME1 = MException('InconsistentInput:InvalideNumberOfArguments',...
                 'Invalid number of input arguments');

if strcmp(model,'BlackScholes')
    if length(param) == 1
    % List of parameters
    sigma = param(1);
    % Give parameters to characteristic function
    resu = chf_bs(u,t,r,q,sigma);
    else
        throw(ME1)
    end

elseif strcmp(model,'Merton')
    if length(param) == 4
    % List of parameters
    sigma = param(1);
    muJ = param(2);
    sigmaJ = param(3);
    lambda = param(4);
    % Give parameters to characteristic function
    resu = chf_merton(u,t,r,q,sigma,muJ,sigmaJ,lambda);
    else 
        throw(ME1)
    end

elseif strcmp(model,'Kou')
    if length(param) == 5
    % List of parameters
    sigma = param(1);
    p1 = param(2);
    eta1 = param(3);
    eta2 = param(4);
    lambda = param(5);
    % Give parameters to characteristic function
    resu = chf_kou(u,t,r,q,sigma,p1,eta1,eta2,lambda);
    else 
        throw(ME1)
    end

elseif strcmp(model,'VG')
    if length(param) == 3
    % List of parameters
    sigmaVG = param(1);
    beta = param(2);
    theta = param(3);
    % S0 = param(4);
    % Give parameters to characteristic function
    resu = chf_vg(u,t,r,q,sigmaVG,beta,theta);
    % resu = chf_vg(u,t,r,q,sigmaVG,beta,theta,S0);
    else
        throw(ME1)
    end

elseif strcmp(model,'Heston')
    if length(param) == 5
    % List of parameters
    kappa = param(1);
    theta = param(2);
    eta = param(3);
    rho = param(4);
    v0 = param(5);
    % Give parameters to characteristic function
    resu = chf_heston(u,t,r,q,kappa,theta,eta,rho,v0);
    else
        throw(ME1)
    end

elseif strcmp(model,'Bates')
    if length(param) == 8
    % List of parameters
    kappa = param(1);
    theta = param(2);
    eta = param(3);
    rho = param(4);
    v0 = param(5);
    muJ = param(6);
    sigmaJ = param(7);
    lambda = param(8);
    % Give parameters to characteristic function
    resu = chf_bates(u,t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda);
    else
        throw(ME1)
    end

% elseif strcmp(model,'WMSV') 
    % if length(param) == 1
    % List of parameters
    % TODO
    % Give parameters to characteristic function
    % TODO
    % else
    %     throw(ME1)
    % end
else
    warning('Method not found')
end


end


%-------------------------------------------------------------------------%
% Here we effectively compute the characteristic function of various model
%  
% Below we can add more functions
%-------------------------------------------------------------------------%

% Black Scholes
function cf = chf_bs(u,t,r,q,sigma)
    i = complex(0,1);
    % Drift correction term
    mu = (r - q) - 0.5*sigma^2;
    % ChF for the Black Scholes model
    cf = exp(i*u*mu*t - 0.5*sigma^2*u.^2*t);  
end

% Merton Jump diffusion
function cf = chf_merton(u,t,r,q,sigma,muJ,sigmaJ,lambda)
    i = complex(0,1);
    % Term for E[exp(J) - 1]
    EX1 = exp(muJ + 0.5 * sigmaJ * sigmaJ) - 1.0;
    % Term for E[exp(iuJ) - 1]
    EX2 = exp(i*u*muJ - 0.5 * sigmaJ * sigmaJ * u.^2) - 1.0;
    % Drift correction term
    mu = (r - q) - 0.5 * sigma^2 - lambda * EX1;
    % ChF for the Merton model
    cf = exp(i*u*mu*t - 0.5*sigma*sigma*u.^2 + lambda*t*EX2);
end

% Kou 
function cf = chf_kou(u,t,r,q,sigma,p1,eta1,eta2,lambda)
    i = complex(0,1);
    p2 = 1 - p1;
    % Term for E[J]
    EX1 = (p1 * eta1)/(eta1 - 1.0) + (p2 * eta2)/(eta2 + 1.0);
    % Term for E[iuJ]
    EX2 = (p1 * eta1) ./ (eta1 - i*u) + (p2 * eta2) ./ (i*u + eta2);
    % Drift correction term
    mu = (r - q) - 0.5 * sigma^2 - lambda * EX1; 
    % Chf for the Kou model
    cf = exp(i*u*mu*t - 0.5*sigma*sigma*u.^2 + lambda*t*EX2);
end

% Variance Gamma [VG]
function cf = chf_vg(u,t,r,q,sigmaVG,beta,theta)
    i = complex(0,1);
    % Correction for martingality
    omega = (1/beta) * log(1 - theta*beta - 0.5*sigmaVG^2*beta);
    % Drift correction term
    mu = (r - q) + omega;
    % Levy process part
    psi = (1 - i*u*theta*beta + 0.5*sigmaVG^2*beta*u.^2) .^ (-t/beta); 
    % Chf for Variance Gamma model
    cf = exp(i*u*mu*t) .* psi;
    % cf = exp(i*u*mu*t + i*u*log(S0)) .* psi;
end

% Heston Stochastic Volatility (Little Trap, see Albrecher et. al)
function cf = chf_heston(u,t,r,q,kappa,theta,eta,rho,v0)
    i = complex(0,1);
    % Functions D and g
    D = sqrt((-kappa + i*rho*eta.*u).^2 + (u.^2 + i*u)*eta^2 );
    g = (kappa - i*rho*eta*u - D) ./ (kappa - i*rho*eta*u + D);
    % Complex-valued functions A and C
    C = (1/eta^2) * (1-exp(-D*t))./(1-g.*exp(-D*t)) .* (kappa - eta*rho*i*u - D);
    A = i*u*(r-q)*t + kappa*theta/eta^2 * (kappa - eta*rho*i*u - D)*t... 
        - kappa*theta/eta^2 * 2*log((1-g.*exp(-D*t))./(1-g));
    % ChF for the Heston model
    cf = exp(A + C * v0);
end

% Bates Jump diffusion stochastic volatility
function cf = chf_bates(u,t,r,q,kappa,theta,eta,rho,v0,muJ,sigmaJ,lambda)
    i = complex(0,1);
    % Functions D and g
    D = sqrt((-kappa + i*rho*eta.*u).^2 + (u.^2 + i*u)*eta^2 );
    g = (kappa - i*rho*eta*u - D) ./ (kappa - i*rho*eta*u + D);
    % Complex-valued functions A and C
    C = (1/eta^2) * (1-exp(-D*t))./(1-g.*exp(-D*t)) .* (kappa - eta*rho*i*u - D);
    A = i*u*(r-q)*t + kappa*theta/eta^2 * (kappa - eta*rho*i*u - D)*t... 
        - kappa*theta/eta^2 * 2*log((1-g.*exp(-D*t))./(1-g));
    % Adjustment for the Bates model
    A = A - lambda*i*u*t*(exp(muJ+1/2*sigmaJ^2)-1) + lambda*t*(exp(i*u*muJ-1/2*sigmaJ^2*u.^2)-1);
    % ChF for the Bates model
    cf = exp(A + C * v0);
end

% WMSV
function cf = chf_wmsv(u,t,r,q,M,R,Q,Sigma0,beta)
    i = complex(0,1);
    resu = zeros(1,length(u));
    rot = [0, 0]; 
    mki = zeros(2,length(u));     
    for k = 1:length(u) 
        w = u(1,k);
        MATEXP = expm(t*[M, -2.0*(Q'*Q); ...
                 0.5*i*w*(i*w - 1)*eye(2), -1.0*(M' + 2.0*i*w*R*Q)]); % step(3)
        % A11 = [MATEXP(1:2,1:2)];
        % A12 = [MATEXP(1:2,3:4)];
        A21 = [MATEXP(3:4,1:2)];
        A22 = [MATEXP(3:4,3:4)];
        [Pk,Dk] = eig(A22); 
        % Vector that contain the eigenvalues 
        lambda_eig = Dk(sub2ind(size(Dk),1:size(Dk,1),1:size(Dk,2))); 
        % DkLog = logm(Dk);
        DkLog = zeros(2,2); n = size(DkLog,1);
        for j = 1:n % step (5)
            val_eig = lambda_eig(1,j);
            % Step (6) - evaluate the complex logarithm
            dki = log(val_eig);
            % Step (7) - produce the sawtooth-like function
            mki(j,k) = mod(imag(dki),pi);
            % Step (8) - check whether rotation has occured
            if w == u(1,1) 
                DkLog = logm(diag(lambda_eig));
                break
            else
                if mki(j,k) - mki(j,k-1) < 0.5*pi
                    rot(j) = rot(j) + 1; % positive rotation
                elseif mki(j,k) - mki(j,k-1) < -0.5*pi
                    rot(j) = rot(j) - 1; % negative rotation
                end     
                % ---------- step (9) -----------% 
                % compute the correct branch of the imaginary part of complex log
                IMdki = mki(j,k) + pi * rot(j);      
                DkLog(1*j,1*j) = real(dki) + i*IMdki;
            end
        end
        % Step (10) - Compute the correct complex matrix log
        A22Log = Pk * DkLog * inv(Pk);
        % Computation of matrix function A(tau)
        A = A22\A21;
        % Computation of scalar function C(tau)
        % A22log = logm(A22);
        C = -0.5*beta * trace(A22Log + t*(M' + 2*R*Q)) + i*w*(r-q);
        % Characteristic Function for WMSV model
        resu(1,k) = exp(trace(A*Sigma0) + C);
    end
end
