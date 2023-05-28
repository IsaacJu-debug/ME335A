display("Elasticity problems with P2 finite elements")
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
% lve:      local connectivity array for an element
% xe:       coordinate of the vertices of an element
% fe:       value of ff in an element
% ke:       value of lambda in an element
% Ke:       element stiffness matrix
% Fe:       element load vectorcalpkg load msh

%% Build a mesh
[X, LV, BE, BN]=CP3Mesh(0.08,1);
fname='u4';
nel=size(LV,2);
npe=size(LV,1);
nod=size(X,2);
dd=size(X,1);
nunk=2*nod;
nbe=size(BE,2);
% Define local-to-global map and associated quantities
LG=[LV; LV+nod];

%%
%% finite element solver begins
%% 

%% Parameters of the problem
% boundary values
% <<-------------------------------------------------------------------->>
% Complete with the value of EtaG, GG.
EtaG=; % Constrained indices
GG=;   % Value of u for each constrained index
% <<-------------------------------------------------------------------->>

% material parameters
EE  = 1e6*ones(1,nel);
nu = 0.45*ones(1,nel);
rho = 10*ones(1,nel);
grav = [0;0];
bb = grav*rho;

%% assembly, from local to global
K=sparse(nunk,nunk);
F=zeros(nunk,1);
for iel=1:nel
% setting the local data
  lge=LG(:,iel);
  lve=LV(1:npe,iel);
  xe(1:dd,1:npe)=X(1:dd,lve(1:npe));
  Ee=EE(iel);
  nue=nu(iel);
  be=bb(:,iel);
% computing element K and F
  [Ke, Fe]=elementKandF(xe,Ee,nue,be);
% assembly, from local to global
  K(lge,lge)=K(lge,lge)+Ke;
  F(lge)=F(lge)+Fe;
end
% nodes with specified value
ng=length(EtaG); 
II=sparse([1:nunk],[1:nunk],1,nunk,nunk);
for ig=1:ng
 K(EtaG(ig),:)=II(EtaG(ig),:);
 F(EtaG(ig))=GG(ig);
end

%% solve algebraic system
U=K\F;

%% plot
% Displacements
Ux=U(1:nod);
Uy=U(nod+1:2*nod);
% Split the triangles into 4 triangles to plot the deformed mesh.
LVS=[LV([1 4 6],:),LV([2 5 4],:),LV([3 6 5],:),LV([4 5 6],:)];
alpha=1; % You can use this to enlarge the displacements
figure(3)
clf
triplot (LVS(1:3,:)', X(1,:), X(2,:),"linewidth",1);
grid on
%axis([0 1.1 -0.8 0.1])
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 18)
hold on
triplot (LVS(1:3,:)', X(1,:)+alpha*Ux', X(2,:)+alpha*Uy',"linewidth",1,'color','red');
title('Reference and deformed configurations');

% Displacement contours
figure(2)
clf
trimesh(LVS(1:3,:)',X(1,:),X(2,:),sqrt(Ux.^2+Uy.^2),"facecolor", "interp","linewidth",1)
view(2)
grid on
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 18)
title('Contours of the norm of the displacement field');

% Stress contours - The stress is discontinuous across element boundaries,
% just as the strains
% We split each element so that we can plot the discontinuous contours
Xn=X(:,LV(1:3,:));
Un=U(LG(:,:));
LVn=reshape([1:size(Xn,2)],[3,nel]);
Strains=zeros(2,2,size(Xn,2));
Stresses=Strains;
for iel=1:nel
  Ee=EE(iel);
  nue=nu(iel);
  lix=[iel*3-2:iel*3];
  [Strains(:,:,lix), Stresses(:,:,lix)]=computeStrainAndStress(Xn(:,lix),Un(:,iel),Ee,nue);
end
figure(4)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Strains(2,2,:),"facecolor", "interp","linewidth",1)
title('epsilon(2,2) component of the strain')
figure(5)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(2,2,:),"facecolor", "interp","linewidth",1)
title('Sigma(2,2) component of the stress')
%view([180, 0]);
figure(6)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(1,1,:),"facecolor", "interp","linewidth",1)
title('Sigma(1,1) component of the stress')
figure(7)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(1,2,:),"facecolor", "interp","linewidth",1)
title('Sigma(1,2) component of the stress')




%% Element matrix and load
function [Ke,Fe]=elementKandF(xe,Ee,nue,be);
% <<-------------------------------------------------------------------->>
% Complete by computing Ke and Fe with a 3-point quadrature
% <<-------------------------------------------------------------------->>
end


%% Stresses and Strains
function [Strains, Stresses]=computeStrainAndStress(xe,ue,Ee,nue)
% Returns the strains and stresses in the element at the corners
% A lot of repeated code, only done for simplicity of explanation.
% <<-------------------------------------------------------------------->>
% Complete by computing Strains and Stresses 
% <<-------------------------------------------------------------------->>
end