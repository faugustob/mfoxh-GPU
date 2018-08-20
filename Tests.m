% Trivariate example
z = [3 2 0.5];
an = [1.5];
ap = [2];
Alphan = [1 ; 1 ; 1];
Alphap = [1 ; 1 ; 1];
bq = [2];
Betaq = [1 ; 1 ; 1];
Contour = mfoxcontour(10, 3, an, Alphan, [],[0;1],[],[3;1],[],[1;1]);
H1 = mfoxh(z, Contour, an, Alphan, ap, Alphap, bq, Betaq,[],[],[0;1],[],[],[],[3;1],[],[],[],[1;1],[])

% H1 = 0.4886 + 0.0035i

% Bivariate example
z = [3 2];
an = [1.5];
ap = [2];
Alphan = [1 ; 1];
Alphap = [1 ; 1];
bq = [2];
Betaq = [1 ; 1];
Contour = [-1.5-10i -1.5+10i ; 2.5-10i 2.5+10i]; % contour set manually
H2 = mfoxh(z, Contour, an, Alphan, ap, Alphap, bq, Betaq,[],[],[0;1],[],[],[],[3;1],[])
% H2 = -0.6014 + 0.0011i
