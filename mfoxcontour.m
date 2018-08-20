function c = mfoxcontour(W, dim, an, Alphan, varargin)
% an = [a1,...,an] and Alphan = [alpha,1,1 ... alpha,n,1;...; alpha,1,r...alpha,n,r]                      
% varargin form (i = 1 ..r): [ci,1 ...ci,n ; gammai,1 ... gammai,n], 
% [di,1 ...di,m ; deltai,1 ... deltai,m]
% See notation in A. Mathai, The H-function, Theory and Applications, Annex A.1
% W   : control the width of the integration interval in [-i\infty +i\infty]
% dim : stands for the dimension

Nvar    = length(varargin);
epsilon = 1/10;
f  = ones(1,dim);
Q  = -Alphan.';
b  = 1-an-epsilon;
lb = [];
ub = [];

for i = 1 : Nvar/2    
 cni = cell2mat(varargin(2*(i-1)+1)); % [c1,i ...cn,i;gamma_1,i...gamma_n,i]
 if(isempty(cni)) cni = [-1e10;1]; end
 dmi = cell2mat(varargin(2*(i-1)+2));  % [d1,i ...dm,i;delta_1,i...delta_m,i]
 if(isempty(dmi)) dmi = [1e10;1]; end
 lb  = [lb max((cni(1,:)-1)./cni(2,:))]+epsilon;
 ub  = [ub min(dmi(1,:)./dmi(2,:))]-epsilon;    
end

options = optimoptions('linprog','Algorithm','interior-point');
out = linprog(f, Q, b, [], [], lb, ub, options);
c = [out.' - 1i * W; out.' + 1i * W];  
