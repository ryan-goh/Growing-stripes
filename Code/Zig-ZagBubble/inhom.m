function [G,J]=inhom(b,w,k)
  global xw m Lw Lder dxw chi1 mu mu0
  

%      mu=1;
    h=dxw/k;
    % dxwk=h;
    J=0*speye(m,m);
    
    
    G = (b*mu0 + b*mu)*w;
    J = -b*mu0 + b*mu; %spdiags(mu0*(-b) + H'*b,0,m,m);
    
    
    
end