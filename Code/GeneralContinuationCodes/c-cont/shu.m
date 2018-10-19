%
function [F,J] = shu(u,kloc,dxu,Lapu,n,mu0,tau)
%global 
  Lu = -(kloc^2*dxu^(-2)*Lapu+speye(n))^2+mu0*speye(n);
%-(kloc^4*dxu^(-2)*Lapu2 + 2*kloc^2*dxu^(-2)*Lapu+speye(n))
% Right-hand side
%F =  Lu*u  + tau*u.^2 - u.^3;
%J = Lu + spdiags( 2*tau*u - 3*u.^2, 0, n, n);
F =  Lu*u  - u.^3;
J = Lu + spdiags( - 3*u.^2, 0, n, n);

end
