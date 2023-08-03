% Compute the cotangent weight Laplacian.
% W is the symmetric cot Laplacian, and A are the area weights 
function [W, A] = cotLaplacian(mesh, L23, L13, L12)
X = mesh.vertices;
T = mesh.faces;
nv = size(X,1);

inputL = exist('L23', 'var') + exist('L13', 'var') + exist('L12', 'var');
if inputL < 3
    % Find orig edge lengths and angles
    L1 = normv(X(T(:,2),:)-X(T(:,3),:));
    L2 = normv(X(T(:,1),:)-X(T(:,3),:));
    L3 = normv(X(T(:,1),:)-X(T(:,2),:));
else
    L1 = L23;
    L2 = L13;
    L3 = L12;
end
EL = [L1,L2,L3];
A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
A = [A1,A2,A3];
A = acos(A);

% The Cot Laplacian 
I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
S = 0.5*cot([A(:,3);A(:,1);A(:,2)]);
In = [I;J;I;J];
Jn = [J;I;I;J];
Sn = [-S;-S;S;S];
W = sparse(double(In),double(Jn),Sn,nv,nv);

if inputL < 3
    % Use the Barycentric areas
    M = mass_matrix_barycentric(mesh);
    A = sum(M,2);
else
    M = mass_matrix_barycentric(mesh, L1, L2, L3);
    A = sum(M,2);
end


function nn = normv(V)
nn = sqrt(sum(V.^2,2));

function M = mass_matrix_barycentric(mesh, L1, L2, L3)

T = mesh.faces; 
inputL = exist('L1', 'var') + exist('L2', 'var') + exist('L3', 'var');
if inputL < 3
    Ar = mesh.ta;
else
    s=(L1+L2+L3)/2;
    Ar = sqrt(s.*(s-L1).*(s-L2).*(s-L3)); 
end
nv = mesh.nv;

I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
Mij = 1/12*[Ar; Ar; Ar];
Mji = Mij;
Mii = 1/6*[Ar; Ar; Ar];
In = [I;J;I];
Jn = [J;I;I];
Mn = [Mij;Mji;Mii];
M = sparse(double(In),double(Jn),Mn,nv,nv);
