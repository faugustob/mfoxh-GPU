function out = mfoxh(z, Contour, an, Alphan, ap, Alphap, bq, Betaq, varargin)
% For dim >= 4, we recommend the use of GPU-enabled HPC servers
% an = [a1,...,an] and Alphan = [alpha,1,1 ... alpha,n,1;...; alpha,1,r...alpha,n,r]                      
% varargin form (i = 1 ..r): [ci,1 ...ci,n ; gammai,1 ... gammai,n], 
%[ci,n+1 ...ci,p ; gammai,n+1 ... gammai,p], 
%[di,1 ...di,m ; deltai,1 ... deltai,m], [di,m+1 ...di,q ; deltai,m+1 ... deltai,q]
% See notation in A. Mathai, The H-function, Theory and Applications, Annex A.1
%================================================================================================

out    =  Integrate(@(s)mfox(s), Contour);
%========================================== Integrand ===========================================

function f = mfox(s)
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

end
