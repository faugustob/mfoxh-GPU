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
