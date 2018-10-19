function F = Lw(kx,c)
global BDx BD2x BDy BD2y mu N ky mu0 tau

F = -(speye(N) + BD2x + ky^2*BD2y)^2 + mu0*speye(N) + c*(BDx + kx*BDy);