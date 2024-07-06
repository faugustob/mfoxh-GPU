function f = mfoxIntegrand(x,y,z)
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