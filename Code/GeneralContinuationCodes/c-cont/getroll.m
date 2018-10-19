function U=getroll(kloc)
global dxu Lapu  xu mu0 n tau
  ustar = sqrt(abs(4/3*(mu0-(1-kloc^2)^2)))*cos(xu)';
  shu_rhs = @(u) shu(u,kloc,dxu,Lapu,n,mu0,tau);
  options = optimset('Jacobian','on','Algorithm','levenberg-marquardt','Display','off');
  U = fsolve(shu_rhs,ustar,options);
end
