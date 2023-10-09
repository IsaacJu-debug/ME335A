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
[X, LV, BE, BN]=CP3Mesh(0.1,1.0);
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
gy2=0.00;

EtaG=[BN(1,BN(2,:)==2),BN(1,BN(2,:)==4)]; % Constrained indices
GG=zeros(size(EtaG,2),2);   % Value of u for each constrained index

Etag2=[BN(2,BN(2,:)==2),BN(2,BN(2,:)==4)];
% Only Uy displacement is non zeros on
% Dirichlet boundary corresponding to 2
GG(Etag2==2,2)=gy2; 
% <<-------------------------------------------------------------------->>

% material parameters
EE  = 0.001*1e6*ones(1,nel);
nu = 0.45*ones(1,nel);
rho = 10*ones(1,nel);
grav = [0;-9.81];
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
    F(EtaG(ig))=GG(ig,1);
    K(EtaG(ig)+nod,:)=II(EtaG(ig)+nod,:);
    F(EtaG(ig)+nod)=GG(ig,2);
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
legend('Reference','Deformed')
xlabel('X1 axis') 
ylabel('X2 axis') 
title('Reference and deformed configurations');

% Displacement contours
figure(2)
clf
trimesh(LVS(1:3,:)',X(1,:),X(2,:),sqrt(Ux.^2+Uy.^2),"facecolor", "interp","linewidth",1)
colorbar
view(2)
grid on
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 18)
xlabel('X1 axis') 
ylabel('X2 axis') 
title('Contours of the norm of the displacement field');

% Stress contours - The stress is discontinuous across element boundaries,
% just as the strains
% We split each element so that we can plot the discontinuous contours
Xn=X(:,LV(1:3,:));
Un=U(LG(:,:));
LVn=reshape([1:size(Xn,2)],[3,nel]);
Strains=zeros(2,2,size(Xn,2));
Stresses=zeros(2,2,size(Xn,2));
for iel=1:nel
  Ee=EE(iel);
  nue=nu(iel);
  lix=[iel*3-2:iel*3];
  [Strains(:,:,lix), Stresses(:,:,lix)]=computeStrainAndStress(Xn(:,lix),Un(:,iel),Ee,nue);
end
figure(4)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Strains(2,2,:),"facecolor", "interp","linewidth",1)
colorbar
% zlim([0.024 0.026])
title('epsilon(2,2) component of the strain')
xlabel('X1 axis') 
ylabel('X2 axis') 
figure(5)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(2,2,:),"facecolor", "interp","linewidth",1)
colorbar
% zlim([0.024 0.026])
title('Sigma(2,2) component of the stress')
xlabel('X1 axis') 
ylabel('X2 axis') 
%view([180, 0]);
figure(6)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(1,1,:),"facecolor", "interp","linewidth",1)
colorbar
% zlim([0.024 0.026])
title('Sigma(1,1) component of the stress')
xlabel('X1 axis') 
ylabel('X2 axis') 
figure(7)
clf
trimesh(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(1,2,:),"facecolor", "interp","linewidth",1)
colorbar
% zlim([0.024 0.026])
title('Sigma(1,2) component of the stress')
xlabel('X1 axis') 
ylabel('X2 axis') 




%% Element matrix and load
function [Ke,Fe]=elementKandF(xe,Ee,nue,be)
% <<-------------------------------------------------------------------->>
% Complete by computing Ke and Fe with a 3-point quadrature
mue=Ee/(1+nue);
lambdae=mue*nue/(1-2*nue);
dNdL=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
      xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
A2 = ((xe(2,2)-xe(2,3))*(xe(1,1)-xe(1,2))+(xe(1,3)-xe(1,2))*(xe(2,1)-xe(2,2)));
[dNx,Nx]=P2elementsX(xe,dNdL,A2);
[dNy,Ny]=P2elementsY(xe,dNdL,A2);
dNa=zeros(2,24,3); dEa=zeros(2,24,3); BB = zeros(2, 24, 3);
Na=zeros(2,12,3);

 % Define the 3 Gauss quadrature points and weights.
 gp = [1/6, 1/6, 2/3; 1/6, 2/3, 1/6; 2/3, 1/6, 1/6];
    
% for jj=1:3
%     dNa(:,:,jj)=[dNx(xe(:,jj)), dNy(xe(:,jj))];
%     Na(:,:,jj)=[Nx(xe(:,jj)), Ny(xe(:,jj))];
% end

for jj=1:3
    xq = xe(:, 1: 3)*gp(jj, :)'; % get the coordinates of quadrature points
    dNa(:,:,jj)=[dNx(xq), dNy(xq)];
    Na(:,:,jj)=[Nx(xq), Ny(xq)];
end

Ke=zeros(12,12); Fe=zeros(12,1);
for jj=1:3
    for kk=1:12
        Ui=dNa(:,2*kk-1:2*kk,jj);
        Du=sum(sum(eye(2).*Ui));
        dEa(:,2*kk-1:2*kk,jj)=mue*(Ui+Ui')/2+lambdae*Du*eye(2);
        BB(:,2*kk-1:2*kk,jj) = (Ui+Ui')/2;
    end
end
% for i1=1:12
%     for jj=1:3
%         for i2=1:12
%             Ke(i2, i1)=Ke(i2, i1)+sum(sum(dEa(:,2*i1-1:2*i1,jj).*...
%                 dNa(:,2*i2-1:2*i2,jj)))*A2/6;
%         end
%         Fe(i1)=Fe(i1)+dot(be,Na(:,i1,jj))*A2/6;
%     end
% end

for i1=1:12
    for jj=1:3
        for i2=1:12
            %Ui=dNa(:,2*i2-1:2*i2,jj);
%             if (i2 == 5)
%                 disp(sum(sum(dEa(:,2*i1-1:2*i1,jj).*...
%                 BB(:,2*i2-1:2*i2,jj))))
%             end
            Ke(i2, i1)=Ke(i2, i1)+sum(sum(dEa(:,2*i1-1:2*i1,jj).*...
                BB(:,2*i2-1:2*i2,jj)))*A2/6;
        end
        Fe(i1)=Fe(i1)+dot(be,Na(:,i1,jj))*A2/6;
    end
end

% Neumann Boundary conditions for all our problems are zero
% As traction free surfaces are considered.
% <<-------------------------------------------------------------------->>
end


%% Stresses and Strains
function [Strains, Stresses]=computeStrainAndStress(xe,ue,Ee,nue)
% Returns the strains and stresses in the element at the corners
% A lot of repeated code, only done for simplicity of explanation.
% <<-------------------------------------------------------------------->>
% Complete by computing Strains and Stresses 
mue=Ee/(1+nue);
lambdae=mue*nue/(1-2*nue);
dNdL=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
      xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
A2 = ((xe(2,2)-xe(2,3))*(xe(1,1)-xe(1,2))+(xe(1,3)-xe(1,2))*(xe(2,1)-xe(2,2)));
[dNx,Nx]=P2elementsX(xe,dNdL,A2);
[dNy,Ny]=P2elementsY(xe,dNdL,A2);
dNa=zeros(2,24,3); dEa=zeros(2,24,3);
Na=zeros(2,12,3);
for jj=1:3
    dNa(:,:,jj)=[dNx(xe(:,jj)), dNy(xe(:,jj))];
    Na(:,:,jj)=[Nx(xe(:,jj)), Ny(xe(:,jj))];
end
Strains=zeros(2,2,3); Stresses=zeros(2,2,3);
for jj=1:3
    for kk=1:12
        Ui=dNa(:,2*kk-1:2*kk,jj);
        Du=sum(sum(eye(2).*Ui));
        dEa(:,2*kk-1:2*kk,jj)=mue*(Ui+Ui')/2+lambdae*Du*eye(2);
        Strains(:,:,jj)=Strains(:,:,jj)+ue(kk)*(Ui+Ui')/2;
        Stresses(:,:,jj)=Stresses(:,:,jj)+ue(kk)*dEa(:,2*kk-1:2*kk,jj)/(10^6);
    end
end
% <<-------------------------------------------------------------------->>
end

function [dN,N]=P2elementsX(xe,dNdL,A2)
    th1 = @(x)((xe(2,2)-xe(2,3))*(x(1)-xe(1,2))+(xe(1,3)-xe(1,2))*(x(2)-xe(2,2))...
          )/A2;
    th2 = @(x)((xe(2,3)-xe(2,1))*(x(1)-xe(1,3))+(xe(1,1)-xe(1,3))*(x(2)-xe(2,3))...
          )/A2;
    th3 = @(x)((xe(2,1)-xe(2,2))*(x(1)-xe(1,1))+(xe(1,2)-xe(1,1))*(x(2)-xe(2,1))...
          )/A2;
    dN1=@(x) (1/A2).*[4.*th1(x)-1,0,0;0,0,0]*dNdL';
    dN2=@(x) (1/A2).*[0,4.*th2(x)-1,0;0,0,0]*dNdL';
    dN3=@(x) (1/A2).*[0,0,4.*th3(x)-1;0,0,0]*dNdL';

    dN4=@(x) (1/A2).*[4.*th2(x),4.*th1(x),0;0,0,0]*dNdL';
    dN5=@(x) (1/A2).*[0,4.*th3(x),4.*th2(x);0,0,0]*dNdL';
    dN6=@(x) (1/A2).*[4.*th3(x),0,4.*th1(x);0,0,0]*dNdL';

    dN=@(x) [dN1(x), dN2(x), dN3(x), dN4(x), dN5(x), dN6(x)];
    N=@(x) [2*th1(x).*(th1(x)-0.5), 2*th2(x).*(th2(x)-0.5), 2*th3(x).*(th3(x)-0.5),...
            4*th1(x).*th2(x), 4*th3(x).*th2(x), 4*th1(x).*th3(x);
            0,0,0,0,0,0];
end

function [dN,N]=P2elementsY(xe,dNdL,A2)
    th1 = @(x)((xe(2,2)-xe(2,3))*(x(1)-xe(1,2))+(xe(1,3)-xe(1,2))*(x(2)-xe(2,2))...
          )/A2;
    th2 = @(x)((xe(2,3)-xe(2,1))*(x(1)-xe(1,3))+(xe(1,1)-xe(1,3))*(x(2)-xe(2,3))...
          )/A2;
    th3 = @(x)((xe(2,1)-xe(2,2))*(x(1)-xe(1,1))+(xe(1,2)-xe(1,1))*(x(2)-xe(2,1))...
          )/A2;
    dN1=@(x) (1/A2).*[0,0,0;4.*th1(x)-1,0,0]*dNdL';
    dN2=@(x) (1/A2).*[0,0,0;0,4.*th2(x)-1,0]*dNdL';
    dN3=@(x) (1/A2).*[0,0,0;0,0,4.*th3(x)-1]*dNdL';

    dN4=@(x) (1/A2).*[0,0,0;4.*th2(x),4.*th1(x),0]*dNdL';
    dN5=@(x) (1/A2).*[0,0,0;0,4.*th3(x),4.*th2(x)]*dNdL';
    dN6=@(x) (1/A2).*[0,0,0;4.*th3(x),0,4.*th1(x)]*dNdL';
    
    dN=@(x) [dN1(x), dN2(x), dN3(x), dN4(x), dN5(x), dN6(x)];
    N=@(x) [0,0,0,0,0,0;
            2*th1(x).*(th1(x)-0.5), 2*th2(x).*(th2(x)-0.5), 2*th3(x).*(th3(x)-0.5),...
            4*th1(x).*th2(x), 4*th3(x).*th2(x), 4*th1(x).*th3(x)];
end


% matref=[3,-1,-1; 
%     3,-1,-1;
%     3,-1,-1;
%     4,4,0;
%     0,4,4;
%     4,0,4];
% map46=[0,0;0,0;0,0;
%     2,1;3,2;3,1];
