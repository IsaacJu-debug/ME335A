display("Second-order problems with P1 finite elements")
%% Simple FEA program for ME335-Finite Element Analysis
%% Authors: Gustavo Buscaglia and Adrian Lew
%% Names of variables (in order of appearance)
% nel:      number of elements
% nod:      number of vertices in the mesh (nodes in this case)
% X:        coordinates of vertices
% npe:      number of degrees of freedom per element 
% dd:       number of spatial dimensions
% L:        length of the interval
% nunk:     number of unknowns/dimension of Wh ('m' in the notes)
% kk:       number of shape functions in an element (all elements the same)
% LG:       local-to-global map
% Etag:     constrained index set
% GG:       values of constrained components of uh
% hh:       values of natural boundary condition (one per node, non-zero only when needed)
% ff:       values of the right hand side in an element (constant)
% lambda:   reaction coefficient of the equation in an element (constant)
% K:        stiffness matrix
% F:        load vector
% lge:      local-to-global map for an element 
% xe:       coordinate of the vertices of an element
% fe:       value of ff in an element
% lambdae:  value of lambda in an element
% hhe:      value of the natural boundary condition for each node of the element
% Ke:       element stiffness matrix
% Fe:       element load vector

%% Build a mesh  -- Elements are implicitly defined by two consecutive vertices
nel=100;
nod=nel+1;
X=[0:1:nel]/nel;

%% Useful Constants 
npe=2; 
dd=1; 
L=1;

%% Build local-to-global map
nunk=nod;
kk=2;
LG=[];
% <<-------------------------->>
% Complete with the value of LG
% <<-------------------------->>  
  
%% finite element solver begins
%%

%% parameters of the problem
EtaG=[nunk];
GG=[1];
hh=zeros(1,nod);
hh(1)=-1;
ff=-2*ones(1,nel);
lambda=10;

%% assembly, from local to global
K=zeros(nunk,nunk);
F=zeros(nunk,1);

%% Assembly of active indices
for iel=1:nel
  % setting the local data
  lge=LG(:,iel);
  xe(1:dd,1:npe)=X(1:dd,iel:iel+1);
  fe=ff(iel);
  hhe = hh(iel:iel+1);
  
  % computing element K and F
  [Ke, Fe]=elementKandF(xe,fe,lambda,hhe);
  % <<-------------------------->>
  % Complete with the assembly of Ke and Fe into K and F 
  % <<-------------------------->>  
end
  
%% Assembly of constrained indices
ng=length(EtaG);
for ig=1:ng
  % <<-------------------------->>
  % Complete with the assembly of K and F 
  % <<-------------------------->>  
end

%% Solve 
U=K\F;


%% plotting (mesh specific)
% plotting points (uniform spacing dx)
fig1=figure(1);
fig1.Name='Function';
clf(fig1);
dx=L/1000; 
xp=0:dx:L;
nxp=length(xp); 
up=zeros(1,nxp);
dup=zeros(1,nxp);
for kp=1:nxp
 xx=xp(kp); 
 nle=max(find(X<=xx));
 if (xx==L) 
     nle=nod-1; 
 end
 nri=nle+1;
 xl=(xx-X(nle))/(X(nri)-X(nle));
 [N1, N2, dN1, dN2]=P1(xl);
 % Please pay attention to how the values of up and dup are computed
 up(kp)=U(nle)*N1+U(nle+1)*N2;
 dup(kp)=(U(nle)*dN1+U(nle+1)*dN2)/(X(nri)-X(nle));
end
plot(X,U(1:nunk),"or","linewidth",2);
set(gca,"FontName","Arial")
set(gca,"FontSize",16)
hold on
plot(xp,up,"r-","linewidth",2);
hold off
fig2=figure(2);
fig2.Name='Derivative';
clf(fig2);
plot(xp,dup,"r-","linewidth",2);

% Exact solution
uep=zeros(1,nxp);
duep=zeros(1,nxp);
for kp=1:nxp
  % <<-------------------------->>
  % Complete uep and duep with the exact solution at each sampling point
  % <<-------------------------->>  
end
figure(1);
hold on
plot(xp,uep,"b-","linewidth",2);
hold off
figure(2)
hold on
plot(xp,duep,"b-","linewidth",2);
hold off

%% element matrix and load
function [Ke, Fe]=elementKandF(xe,fe,lambdae,hhe)
  h=xe(2)-xe(1);  % Element size
  % <<-------------------------->>
  % Complete with the values of Ke and Fe
  % <<-------------------------->>  
end


%% basis functions (only for plotting)
function [N1, N2, dN1, dN2]=P1(x)
  % <<-------------------------->>
  % Complete with the values of N1(x), N1'(x), N2(x), and N2'(x) over an interval [0,1]
  % <<-------------------------->> 
end



