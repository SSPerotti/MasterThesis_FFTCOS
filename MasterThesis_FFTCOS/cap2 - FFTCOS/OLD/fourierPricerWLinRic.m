function price=fourierPricerWLinRic(S0,K,Sigma0,A_i,A_j,M,Q,Gindikin,R,t,T,r_i,r_j)

%initial log-price
x=log(S0);

%dampening factor
alpha=1.2;

N=4096;                     %
eta=0.18;                   %
lambda=(2*pi)/(N*eta);      %equation 23 Carr and Madan
v=0 : eta : ((N-1)*eta);    %See Carr and Madan page 10
b=(N*lambda)/2;             %Equation 20 Carr and Madan

vv=v-(alpha+1)*1i;

tol=1e-20;
ncf=zeros(1,length(vv));
for j=1:length(vv)
ncf(j)=exp(-r_i*(T-t)).*computeLtW(x,Sigma0,1i.*vv(j),A_i,A_j,M,Q,Gindikin,R,t,T,r_i,r_j);

%if the cf is small then truncate
if abs(ncf(j))<tol
    ncf(j+1:end)=0;
    break
end

%if the cf is NaN then approximate with previous values since cf is
%continuous.
if (isnan(ncf(j))==1 && j == 1) 
        ncf(j)=0;
end

if (isnan(ncf(j))==1 && j ~= 1)
    ncf(j)=ncf(j-1);
end

end



dcf=(alpha^2+alpha-v.^2+1i*(2*alpha+1).*v);
cf=ncf./dcf;

%summation term in equation 22 Carr and Madan
tmp=cf.*exp(1i*v*b)*eta; 

%applying simpson's rule - formula 24 Carr and Madan
jvec=1:N;
tmp=(tmp/3).*(3+(-1).^jvec-((jvec-1)==0) ); 

ft=fft(tmp,N);
ft=ft';
ku=-b : lambda : ((N-1)*lambda-b); 
ku=ku';
cpvec=(exp(-alpha*ku(1:N)).*ft)/pi; %call price vector resulting in equation 22

%perform interpolation so to give the price
kvec=exp(ku);
price=zeros(length(K),1);
for j=1:length(K)
for i=1:N,
if kvec(i)>K(j),
break;
end
end
lowerstrike=kvec(i-1);
higherstrike=kvec(i);
lowercp=real(cpvec(i-1));
highercp=real(cpvec(i));
xp=[lowerstrike higherstrike];
yp=[lowercp highercp];
price(j)=max(real(interp1(xp,yp,K(j))),0);
end
