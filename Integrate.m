function out = Integrate(fun, Contour)
dim    = size(Contour,1);

N      = 2^8 * 3^3 * 5^(2+dim); % For a better performance, N is made of powers of prime numbers
%========================== Multivariate Quasi Monte-Carlo Integration ==========================
p      = haltonset(dim,'Skip',1e3,'Leap',1e2);
p      = scramble(p,'RR2');
in     = gpuArray(net(p,N));
G      = gpuArray(ones(1,N));
C      = gpuArray(Contour);
points = kron(G,C(:,1)) + kron(G, C(:,2)-C(:,1)).* in';
mceval = fun(points);
mcsum  = sum(mceval,2);
v      = prod(C(:,2) - C(:,1)); % volume
out    = v * mcsum / N; % Integral

integrand = @(x,y,z) mfoxIntegrand(x,y,z);

% out = integral3(integrand,C(1,1),C(1,2),C(2,1),C(2,2),C(3,1),C(3,2))
% x=1;
end