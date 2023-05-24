%display("Second-order problems with P1 finite elements")
%% Simple FEA program for ME335-Finite Element Analysis
%% Authors: Gustavo Buscaglia and Adrian Lew
%% Names of variables (more or less in order of appearance)
% nel:      number of elements
% nod:      number of vertices in the mesh (nodes in this case)
% X:        coordinates of vertices
% npe:      number of vertices per element 
% dd:       number of spatial dimensions
% nunk:     number of unknowns/dimension of Wh ('m' in the notes)
% nbe:      number of boundary elements
% kk:       number of shape functions in an element (all elements the same)
% LG:       local-to-global map
% Etag:     constrained index set
% GG:       values of constrained components of uh
% hh:       values of natural boundary condition (one per element on the boundary, non-zero only when needed)
% ff:       values of the right hand side in an element (constant)
% diffcoef: diffusion coefficient of the equation in an element (constant)
% K:        stiffness matrix
% F:        load vector
% lge:      local-to-global map for an element 
% xe:       coordinate of the vertices of an element
% fe:       value of ff in an element
% ke:       value of lambda in an element
% Ke:       element stiffness matrix
% Fe:       element load vectorcalpkg load msh

%% Build a mesh
close all
[X, LV, BE, BN]=CP2Mesh(.1);
% Define local-to-global map and associated quantities
% we take LG=LV
LG=LV;
nel=size(LG,2);
npe=size(LG,1);
nod=size(X,2);
dd=size(X,1);
nunk=nod;
nbe=size(BE,2); %

%%
%% finite element solver begins
%% 

%% Parameters of the problem
% boundary values
h1 = 0.90; % Value of Neumann boundary condition for each edge in BE [Wft^-2]
h2 = 0.0;
h3 = -0.10;
h4 = 0.0;

% <<-------------------------------------------------------------------->>
% Complete with the value of EtaG, GG, hh, and ff, and diffcoeff
indices = (BN(2,:) == 5) | (BN(2,:) == 6);
EtaG= BN(1, indices); % Constrained indices

GG = zeros(size(EtaG));
EtaGPos = BN(2, BN(2,:) > 4);
GG(EtaGPos(:) == 5) = 2500;
GG(EtaGPos(:) == 6) = 100;  % Value of u for each constrained index [C]

hh = zeros(nbe, 1);   % Value of Neumann boundary condition for each edge in BE
hh(BE(3, :) == 1) = h1;
hh(BE(3, :) == 2) = h2;
hh(BE(3, :) == 3) = h3;
hh(BE(3, :) == 4) = h4;
ff = ones(nel,1) * 0.2;   % Value of f for each element [0.2W/ft3]

% material parameters
difcoeff=ones(nel)* 2.0* 10^-4; % Value of k for each element [W ft^-1 C^-1]
% <<-------------------------------------------------------------------->>

%% assembly, from local to global
K=zeros(nunk,nunk);
F=zeros(nunk,1);
for iel=1:nel
% setting the local data
  lge=LG(:,iel);
  xe(1:dd,1:npe)=X(1:dd,lge(1:npe));
  ke=difcoeff(iel);
  fe=ff(iel);
% computing element K and F
  [Ke Fe]=elementKandF(xe,ke,fe);
% assembly, from local to global
  K(lge,lge)=K(lge,lge)+Ke;
  F(lge)=F(lge)+Fe;
end
% nodes with specified value
ng=length(EtaG); 
II=eye(nunk);
for ig=1:ng
 K(EtaG(ig),:)=II(EtaG(ig),:);
 F(EtaG(ig))=GG(ig);
end
% Neumann boundaries
% <<-------------------------------------------------------------------->>
% Assemble Neumann boundary conditions
for ied = 1:nbe
    lged = BE(1:2, ied);
    xed(1:dd, 1:2) = X(1:dd, lged);
    Led = norm(xed(:, 1) - xed(:, 2));
    Hed = hh(ied);
    F(lged) = F(lged) + Hed * Led / 2;
end

% <<-------------------------------------------------------------------->>

%% solve algebraic system
U=K\F;
%% plot
figure 
trisurf(LG',X(1,:),X(2,:),U)
colorbar
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 18)

%% Evaluate the solution

xpArray = [0 , 0.9; 0.25, 0.25];

[u, du]=uValue(xpArray(1,:), X, LV, U);
fprintf('%f %f\n', xpArray(1,1), xpArray(1,2));
disp(u);
disp(du);

[u, du]=uValue(xpArray(2,:), X, LV, U);
fprintf('%f %f\n', xpArray(2,1), xpArray(2,2));
disp(u);
disp(du);
%% element matrix and load
function [Ke, Fe]=elementKandF(xe,ke,fe)
% <<-------------------------------------------------------------------->>
% Complete by computing Ke and Fe
    dN=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
    xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
    Ae2=dN(2,3)*dN(1,2)-dN(1,3)*dN(2,2);
    dN = dN / Ae2;
    Ke = Ae2 / 2 * ke * dN' * dN;
    Fe = Ae2 * fe * ones(3, 1)/6;
% <<-------------------------------------------------------------------->>
end

% Compute shape functions and derivatives to reconstruct function inside
function [NN, pdN]=P1Functions(xe,x)
% <<-------------------------------------------------------------------->>
% Complete by computing the shape functions NN and their gradient dN at x
% NN(a) contains the value of N_a^e(x), and dN is defined in the notes
    NN = zeros(1,3);
    pdN = zeros(2,3); % dN/dxi
    dN=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
    xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
    Ae2=dN(2,3)*dN(1,2)-dN(1,3)*dN(2,2);

    NN(1) = 1./Ae2 * (dN(1, 1)*(x(1) - xe(1, 2)) + dN(2, 1)*(x(2) - xe(2, 2)));
    NN(2) = 1./Ae2 * (dN(1, 2)*(x(1) - xe(1, 3)) + dN(2, 2)*(x(2) - xe(2, 3)));
    NN(3) = 1./Ae2 * (dN(1, 3)*(x(1) - xe(1, 1)) + dN(2, 3)*(x(2) - xe(2, 1)));
    
    if min(NN) < 0
        NN = [0,0,0];
    else
        pdN(1,1) = 1./Ae2 * dN(1, 1);
        pdN(2,1) = 1./Ae2 * dN(2, 1);

        pdN(1,2) = 1./Ae2 * dN(1, 2);
        pdN(2,2) = 1./Ae2 * dN(2, 2);

        pdN(1,3) = 1./Ae2 * dN(1, 3);
        pdN(2,3) = 1./Ae2 * dN(2, 3);
    end
    
% <<-------------------------------------------------------------------->>
end

function [u, du]=uValue(xp, X, LV, U)
% <<-------------------------------------------------------------------->>
% Complete by computing the value of u(xp) and du=[dudx(xp), dudy(xp)]
% Return 0 in u and du if xp is outside the domain
    trias = [192, 220, 237, 39, 101, 115, 29]; 
    i=1;
    %pel = brutalForceSearch(xp, LV, X);
    pel=0;
    while pel==0
        pel=quadtreeRec(xp, LV, X, trias(i));
        i=i+1;
    end
    
    els=LV(:,pel);
    uel=U(els,1);
    xe=X(:,els);
    
    %disp(xe);
    [NN,pdN]=P1Functions(xe,xp);
    u=NN*uel;
    du=pdN*uel;

% <<-------------------------------------------------------------------->>
end

function [pel]=brutalForceSearch(xp, LV, X)
    % Inputs:
    % 
    % xp : A 2D point in the format [x y] for which the containing triangle is to be found.
    % LV : A 3-by-N matrix where each column represents a triangle in the triangulation.
    %         The entries are indices to the points in X that form the vertices of the triangle.
    % X : A 2-by-M matrix where each column represents a point in the 2D space.
    %         The points are the vertices of the triangles in the triangulation. 
    % Output:
    % pel : The column index in LV of the triangle that contains the point xp. 
    %         If no such triangle is found, pel is 0.
    for i = 1:size(LV,2)
        x1 = LV(1,i);x2 = LV(2,i);x3 = LV(3,i); 
        A2 = ((X(2,x2)-X(2,x3))*(X(1,x1)-X(1,x2))+(X(1,x3)-X(1,x2))*(X(2,x1)-X(2,x2)));
        th1 = ((X(2,x2)-X(2,x3))*(xp(1)-X(1,x2))+(X(1,x3)-X(1,x2))*(xp(2)-X(2,x2))...
            )/A2;
        th2 = ((X(2,x3)-X(2,x1))*(xp(1)-X(1,x3))+(X(1,x1)-X(1,x3))*(xp(2)-X(2,x3))...
            )/A2;

        th3 = 1-th1-th2;

        if min([th1, th2, th3])<0
            continue;
        else
            pel = i;
        end
    end
end

function [pel]=quadtreeRec(xp, LV, X, i)
    % Inputs:
    % 
    % xp : A 2D point in the format [x y] for which the containing triangle is to be found.
    % LV : A 3-by-N matrix where each column represents a triangle in the triangulation.
    %         The entries are indices to the points in X that form the vertices of the triangle.
    % X : A 2-by-M matrix where each column represents a point in the 2D space.
    %         The points are the vertices of the triangles in the triangulation.
    % i : The column index in LV of the current triangle being processed.
    % Output:
    % pel : The column index in LV of the triangle that contains the point xp. 
    %         If no such triangle is found, pel is 0.

    x1 = LV(1,i);x2 = LV(2,i);x3 = LV(3,i); 
    A2 = ((X(2,x2)-X(2,x3))*(X(1,x1)-X(1,x2))+(X(1,x3)-X(1,x2))*(X(2,x1)-X(2,x2)));
    th1 = ((X(2,x2)-X(2,x3))*(xp(1)-X(1,x2))+(X(1,x3)-X(1,x2))*(xp(2)-X(2,x2))...
        )/A2;

    th2 = ((X(2,x3)-X(2,x1))*(xp(1)-X(1,x3))+(X(1,x1)-X(1,x3))*(xp(2)-X(2,x3))...
        )/A2;

    th3 = 1-th1-th2;
    if min([th1, th2, th3])<0
        if th1<0
            [row1,col1] = find(LV==x2);
            [row2,col2] = find(LV==x3); 
            [row3,col3] = find(LV==x1);
        elseif th2<0
            [row1,col1] = find(LV==x3);
            [row2,col2] = find(LV==x1); 
            [row3,col3] = find(LV==x2);
        else
            [row1,col1] = find(LV==x1);
            [row2,col2] = find(LV==x2); 
            [row3,col3] = find(LV==x3);
        end
        ic = setdiff(intersect(col1, col2), i);
        if isempty(ic)==1
            pel=0;
        else
            pel = quadtreeRec(xp, LV, X, ic);
        end
    else
        pel = i;
    end
end
