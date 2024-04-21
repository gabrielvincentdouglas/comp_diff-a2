function u = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, gFunc)
N = size(node,1); NT = size(elem,1);
A = sparse(N, N); b = zeros(N,1); 
for j = 1:NT
P1 = node(elem(j, 1), :);
P2 = node(elem(j, 2), :);
P3 = node(elem(j, 3), :);
KEj = elem_stiff (P1, P2, P3, DFunc); % compute local stiffness matrix 
MEj = elem_mass (P1, P2, P3, qFunc);
FEj = elem_load(P1,P2,P3,fFunc); %compute local load vector 
A(elem(j , :), elem(j , :)) = A(elem(j , :), elem(j , :)) + KEj + MEj;
b(elem(j, :), 1) = b(elem(j, :), 1) + FEj ; 
end

xa = min(node(:,1));
xb = max(node(:,1));
ya = min(node(:,2));
yb = max(node(:,2));
eps = 10^(-14);
isLeftBnd = abs(node(:, 1) - xa) < eps;
isRightBnd = abs(node(:, 1) - xb) < eps; 
isButtomBnd = abs(node(:, 2) - ya) < eps; 
isTopBnd = abs(node(:, 2) - yb) < eps;

isBndNode = false(N, 1); 
isBndNode(isLeftBnd) = true; 
isBndNode(isRightBnd) = true; 
isBndNode(isButtomBnd) = true; 
isBndNode(isTopBnd) = true; 
bndNode = find(isBndNode); 
freeNode = find(~isBndNode);
u = zeros(N, 1);
u(bndNode) = gFunc(node(bndNode, :));
b = b - A*u;
u(freeNode) = A(freeNode, freeNode)\b(freeNode);


    function Fej = elem_load(P1, P2, P3, fFunc)
    J = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
    absdetJ = abs(det(J));

    Fej = zeros(3,1);
    for i = 1:3
    Fej(i) = absdetJ*(fFunc(P1)*(i==1)+fFunc(P2)*(i==2)+fFunc(P3)*(i==3))/6;
    end
    end

    function Kej = elem_stiff(P1, P2, P3, DFunc)
    gradN = [-1, 1, 0; -1, 0, 1];
    J = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
    absdetJ = abs(det(J));

    Kej = zeros(3,3);
    for i = 1:3
    for k = 1:3
    Kej(i,k) = absdetJ*(inv(J)'*gradN(:,i))'*(inv(J)'*gradN(:,k))*(DFunc(P1) + DFunc(P2) + DFunc(P3))/6;
    end
    end
    end

    function Mej = elem_mass(P1, P2, P3, qFunc)
    J = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
    absdetJ = abs(det(J));

    Mej = zeros(3,3);
    for i = 1:3
    for k = 1:3
    Mej(i,k) = absdetJ*(qFunc(P1)*(i==1)*(k==1)+qFunc(P2)*(i==2)*(k==2)+qFunc(P3)*(i==3)*(k==3))/6;
    end
    end
    end
end