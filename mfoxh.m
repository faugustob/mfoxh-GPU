function out = mfoxh(z, Contour, an, Alphan, ap, Alphap, bq, Betaq, varargin)
% For dim >= 4, we recommend the use of GPU-enabled HPC servers
% an = [a1,...,an] and Alphan = [alpha,1,1 ... alpha,n,1;...; alpha,1,r...alpha,n,r]                      
% varargin form (i = 1 ..r): [ci,1 ...ci,n ; gammai,1 ... gammai,n], 
%[ci,n+1 ...ci,p ; gammai,n+1 ... gammai,p], 
%[di,1 ...di,m ; deltai,1 ... deltai,m], [di,m+1 ...di,q ; deltai,m+1 ... deltai,q]
% See notation in A. Mathai, The H-function, Theory and Applications, Annex A.1
%================================================================================================
dim    = size(Contour,1);

N      = 2^8 * 3^3 * 5^(2+dim); % For a better performance, N is made of powers of prime numbers
%========================== Multivariate Quasi Monte-Carlo Integration ==========================
p      = haltonset(dim,'Skip',1e3,'Leap',1e2);
p      = scramble(p,'RR2');
in     = gpuArray(net(p,N));
G      = gpuArray(ones(1,N));
C      = gpuArray(Contour);
points = kron(G,C(:,1)) + kron(G, C(:,2)-C(:,1)).* in';
mceval = Integrand(points);
mcsum  = sum(mceval,2);
v      = prod(C(:,2) - C(:,1)); % volume
out    = v * mcsum / N; % Integral
%========================================== Integrand ===========================================
function f = Integrand(s)
j       = sqrt(-1);
r       = length(z);
Nvar    = length(varargin);
for nvar = 1 : Nvar
 if(isempty(cell2mat(varargin(nvar))))
  varargin{nvar} = zeros(2,0);
 end
end
Phi = 1;
for i = 1 : r   
 cni = gpuArray(cell2mat(varargin(4*(i-1)+1)));
 cpi = gpuArray(cell2mat(varargin(4*(i-1)+2)));
 dmi = gpuArray(cell2mat(varargin(4*(i-1)+3)));
 dqi = gpuArray(cell2mat(varargin(4*(i-1)+4)));
 Phi = Phi .* ((GammaProd(1-cni(1,:),cni(2,:), s(i,:)).* GammaProd(dmi(1,:),-dmi(2,:), s(i,:)))...
  ./(GammaProd(cpi(1,:),-cpi(2,:), s(i,:)).* GammaProd(1-dqi(1,:),dqi(2,:), s(i,:)))).* z(i).^s(i,:);
end
Psi     = GammaProd(1-an,Alphan,s)./(GammaProd(ap,-Alphap,s).* GammaProd(1-bq,Betaq,s));
f       = (1/(2*j*pi)^r) * Phi .* Psi;
end
%========================================== GammaProd ===============================================
function output = GammaProd(p,m,s) 
if (isempty(p)|| isempty(m)) 
    output = ones(size(s(1,:)));
else
L1    = size(s,1);    
comb  = 0;
for i = 1 : L1
[pp ss]  = meshgrid(p,s(i,:));
 mm      = meshgrid(m(i,:),s(i,:));
 comb    =  comb + mm .* ss;
end
    output = reshape(prod(gammas(pp + comb),2),size(s(1,:)));
end
end
end
% gammas function here is the complex gamma, available in
% www.mathworks.com/matlabcentral/fileexchange/3572-gamma
end
