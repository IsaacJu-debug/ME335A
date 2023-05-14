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
[X, LV, BE, BN]=CP2Mesh(.1);
% Define local-to-global map and associated quantities
% we take LG=LV
LG=LV;
nel=size(LG,2);
npe=size(LG,1);
nod=size(X,2);
dd=size(X,1);
nunk=nod;
nbe=size(BE,2);

%%
%% finite element solver begins
%% 

%% Parameters of the problem
% boundary values

% <<-------------------------------------------------------------------->>
% Complete with the value of EtaG, GG, hh, and ff, and diffcoeff
EtaG=; % Constrained indices
GG=;   % Value of u for each constrained index
hh=;   % Value of Neumann boundary condition for each edge in BE
ff=;   % Value of f for each element

% material parameters
difcoeff=; % Value of k for each element
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
% <<-------------------------------------------------------------------->>

%% solve algebraic system
U=K\F;
%% plot
figure 
trisurf(LG',X(1,:),X(2,:),U)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 18)

%% Evaluate the solution



%% element matrix and load
function [Ke, Fe]=elementKandF(xe,ke,fe)
% <<-------------------------------------------------------------------->>
% Complete by computing Ke and Fe
% <<-------------------------------------------------------------------->>
end

% Compute shape functions and derivatives to reconstruct function inside
function [NN, dN]=P1Functions(xe,x)
% <<-------------------------------------------------------------------->>
% Complete by computing the shape functions NN and their gradient dN at x
% NN(a) contains the value of N_a^e(x), and dN is defined in the notes
% <<-------------------------------------------------------------------->>
end

function [u, du]=uValue(xp, X, LV, U)
% <<-------------------------------------------------------------------->>
% Complete by computing the value of u(xp) and du=[dudx(xp), dudy(xp)]
% Return 0 in u and du if xp is outside the domain
% <<-------------------------------------------------------------------->>
end
