function [u,X,Y] = poiss_fd(xa,xb,ya,yb,Nx,Ny,fFunc,gFunc)
% Parameters
hx = (xb - xa)/(Nx+1); 
hy = (yb - ya)/(Ny+1);
% Define Grid
x = (xa+hx):hx:(xb-hx);
y = (ya+hy):hy:(yb-hy);
[X, Y] = meshgrid(x, y);
% ----
F = fFunc(X,Y);
F = reshape(F', Nx*Ny,1);
for i = 1:Ny
gl = gFunc(xa,y(i));
F((i-1)*Nx+1) = F((i-1)*Nx+1) + gl/hx^2;
gr = gFunc(xb,y(i));
F(i*Nx) = F(i*Nx) + gr/hx^2;
end
for j = 1:Nx
gl = gFunc(x(i),ya);
F(j) = F(j) + gl/hy^2;
gu = gFunc(x(i),yb);
F(Nx*Ny - Nx + j) = F(j) + gu/hy^2;
end
% ----
INy = spdiags(ones(Ny,1), 0, Ny,Ny);
INx = spdiags(ones(Nx,1), 0, Nx,Nx);
diag = 2*(1/hx^2+1/hy^2)*ones(Nx,1);
offdiag = -1/hx^2*ones(Nx,1);
C = spdiags([offdiag diag offdiag], -1:1, Nx,Nx);
D = -1/hy^2.*INx;
E = spdiags([ones(Ny,1) zeros(Ny,1) ones(Ny,1)], -1:1, Ny,Ny);
A = kron(INy,C) + kron(E,D);
% ----
u = A\F;
u = reshape(u, Nx,Ny)';