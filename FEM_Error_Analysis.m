% Boundary of rectangular domain [a, b]x[c, d]
a = -1;
b = 1;
c = -1;
d = 1;

kmax = 7;
Nvalues = zeros(kmax);
E = zeros(kmax);

for k = 1:kmax
% Number of points in the x and y directions
nx = 2^(k);
ny = 2^(k);
Nvalues(k) = nx;
% List of x and y points
x = linspace(a, b, nx);
y = linspace(c, d, ny);
% Create grid of x and y values
[X, Y] = meshgrid(x, y);
% Turn grid into list we can use to represent node
X = X(:);
Y = Y(:);
node = [X, Y];
% Generate the triangulation from the nodes
elem = delaunay(node);

% Define the coefficient functions, boundary conditions and exact solution
DFunc = @(X)1;
qFunc = @(X)0;
fFunc = @(X)-2*exp(X(1)+X(2));
uExact = @(X)exp(X(:,1) + X(:,2));
gFunc = uExact;

% Solve the PDE using the finite element method
u = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, gFunc);

NT = size(elem,1);
for j = 1:NT
P1 = node(elem(j, 1), :);
P2 = node(elem(j, 2), :);
P3 = node(elem(j, 3), :);
J = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
E(k) = E(k) + det(J) *( (uExact(P1) - u(elem(j, 1)))^2 + (uExact(P2) - u(elem(j, 2)))^2 + (uExact(P3) - u(elem(j, 3)))^2 )/6;
end
E(k)  = sqrt(E(k));
end

loglog(Nvalues,E,'-')