x1 = 0.5;
y1 = 0.5;
x2 = 1;
y2 = 1;
x3 = 0.5;
y3 = 1;

P1 = [x1, y1];
P2 = [x2, y2];
P3 = [x3, y3];

D = @(X)5;
q = @(X)1;
f = @(X)3;

%%% Problem 2: Try these functions after you have it working with constant D, q and f
% D = @(X)1 + X(1).*X(2);
% q = @(X)exp(sqrt(X(1).^2 + X(2).^2));
% f = @(X)sin(X(1)) + cos(X(2));

Kej = elem_stiff(P1, P2, P3, D);
%Mej = elem_mass(P1, P2, P3, q);
Fej = elem_load(P1, P2, P3, f);

function Fej = elem_load(P1, P2, P3, f)
J = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
absdetJ = abs(det(J));

Fej = zeros(3,1);
for i = 1:3
    Fej(i) = absdetJ*(f(P1)*(i==1)+f(P2)*(i==2)+f(P3)*(i==3))/6;
end
end

function Kej = elem_stiff(P1, P2, P3, D)
gradN = [-1, 1, 0; -1, 0, 1];
J = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
absdetJ = abs(det(J));

Kej = zeros(3,3);
for i = 1:3
    for j = 1:3
        Kej(i,j) = absdetJ*(inv(J)'*gradN(:,i))'*(inv(J)'*gradN(:,j))*(D(P1) + D(P2) + D(P3))/6;
    end
end
end
