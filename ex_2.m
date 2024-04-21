clear
% Parameters
xa = -1; xb = 1; ya = -1; yb = 1;
N = 20; % N^2 = interior points
% ----
exact = @(x,y) sin(pi*x).*cos(pi/2*y);
f = @(x,y) (5/4*pi^2*sin(pi*x).*cos(pi/2*y));
g = @(x,y) 0;
[u,X,Y] = poiss_fd(-1,1,-1,1,N,N,f,g);
exactfn = exacsin(pi*X).*cos(pi/2*Y);
% ----
subplot(1,2,1)
surf(X,Y,u)
title('Approximated solution')
subplot(1,2,2)
surf(X,Y,exact)
title('Exact solution')
% Error
kmax = 10;
Nvalues = zeros(kmax,1);
E = zeros(kmax,1);
for k = 1:kmax
    N = 2^(k+1);
    h = 2/(N+1);
    Nvalues(k) = N;
    [u,X,Y] = poiss_fd(-1,1,-1,1,N,N,f,g);
    exact = sin(pi*X).*cos(pi/2*Y); 
    e = h*norm(u - exact);
    E(k) = e;
end

loglog(Nvalues,E,'.-')
axis([0 2^11 1e-7 1e0]), grid on

