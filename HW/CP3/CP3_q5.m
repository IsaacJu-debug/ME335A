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
eccArray = [1.0 , 0.5, 0.1];
hmaxArray = [0.1 , 0.05];

output_folder = 'plots';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

for ecc_i = 1:3
    for h_i = 1:2
        ecc = eccArray(ecc_i);
        hmax = hmaxArray(h_i);
        [X, LV, BE, BN]=CP3Mesh(hmax, ecc);
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

        g2y = 0.05; % q2

        indices = (BN(2,:) == 2) | (BN(2,:) == 4); % Dirichlet boundary
        EtaG= BN(1, indices); % Constrained indices

        GG = zeros(size(EtaG, 2), 2);
        EtaGPos = BN(2, indices);

        GG(EtaGPos(:) == 2, 2) = g2y;
        GG(EtaGPos(:) == 4, 1) = 0.0;  % Value of u for each constrained index [C]
        GG(EtaGPos(:) == 4, 2) = 0.0;  % Value of u for each constrained index [C]

        % as there is no traction applied, "do nothing boundary is applied"

        % <<-------------------------------------------------------------------->>

        % material parameters
        EE  = 1e6*ones(1,nel); % Young's modulus
        %nu = 0.0*ones(1,nel); % q2 poisson's ratio 
        nu = 0.45*ones(1,nel); % q3 poisson's ratio 
        rho = 10*ones(1,nel); % kg/m3
        grav = [0;-9.81]; % m/s^2
        bb = grav*rho;

        %% debugging local stiffness matrix
        % xe = [ 0.2170, 0.3072, 0.3052, 0.2621, 0.3062, 0.2611;
        %         -0.0516, -0.0978, -2.4275e-4, -0.0747, -0.0490, -0.0259];
        % iel = 1;
        % be=bb(:,iel);
        % Ee=EE(iel);
        % nue=nu(iel);   
        % [Ke,Fe]=elementKandF(xe,Ee,nue,be);

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
        colorbar

        figure(5)
        clf
        trisurf(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(2,2,:),"facecolor", "interp","linewidth",1)
        title(sprintf('Sigma(2,2) component of the stress\nEccentricity: %.2f, Max Mesh Size: %.2f', ecc, hmax));
        view([180, 0]);
        colorbar
        
        % Save the figure with both subplots as a PNG file
        output_filename = sprintf('ecc_%.2f_hmax_%.2f.png', ecc, hmax);
        full_output_path = fullfile(output_folder, output_filename);
        saveas(gcf, full_output_path);

    end
end
% 
% figure(6)
% clf
% trisurf(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(1,1,:),"facecolor", "interp","linewidth",1)
% title('Sigma(1,1) component of the stress')
% colorbar
% 
% figure(7)
% clf
% trisurf(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(1,2,:),"facecolor", "interp","linewidth",1)
% title('Sigma(1,2) component of the stress')
% colorbar
%%
% figure(8)
% clf
% trisurf(LVn(:,:)',Xn(1,:),Xn(2,:),Strains(2,2,:),"facecolor", "interp","linewidth",1)
% title('epsilon(2,2) component of the strain')
% colorbar
% 
% figure(9)
% clf
% trisurf(LVn(:,:)',Xn(1,:),Xn(2,:),Stresses(2,2,:),"facecolor", "interp","linewidth",1)
% title('Sigma(2,2) component of the stress')
% colorbar

%% Element matrix and load
function [Ke,Fe]=elementKandF(xe,Ee,nue,be)
% <<-------------------------------------------------------------------->>
% Complete by computing Ke and Fe with a 3-point quadrature
% input: xe coordintes [2,6]
%        Ee Young's modulus
%        nue Poisson's ratio
%        be Body force [2,1]
    nodeNum = size(xe, 2);
    dofNum = size(xe, 2)*2;
    Ke = zeros(dofNum, dofNum);
    Fe = zeros(dofNum, 1);
    
    mue=Ee/(1+nue)/2; %
    lambdae=mue*nue/(1-2*nue); %
    
    dLdx=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
          xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
    A2 = ((xe(2,2)-xe(2,3))*(xe(1,1)-xe(1,2))+(xe(1,3)-xe(1,2))*(xe(2,1)-xe(2,2)));
    
    nq = 3;
    % Define the 3 Gauss quadrature points and weights.
    gp = [1/6, 1/6, 2/3; 1/6, 2/3, 1/6; 2/3, 1/6, 1/6];
    w = [1/3, 1/3, 1/3]*A2/2;

    [dNdx,N1x]=P2elementsX1(xe,dLdx,A2); % dNdx [2, 12]; N1x [2, 6]
    [dNdy,N1y]=P2elementsX2(xe,dLdx,A2);
    
    for i = 1:nq
        xq = xe(:, 1: 3)*gp(i, :)'; % get the coordinates of quadrature points
        % Calculate the shape functions and their derivatives at the quadrature point.
        BB = zeros(2,2,dofNum);
        DD = zeros(dofNum,1);
        SS = zeros(2,2,dofNum);

        dNdx_q = dNdx(xq);
        dNdy_q = dNdy(xq);
        
        N1x_q = N1x(xq);
        N1y_q = N1y(xq);
        for iNode = 1:nodeNum
%             if (iNode == 5)
%                 disp(dNdx_q)
%             end
%             
            BB(1:2, 1:2, iNode) = [dNdx_q(1, iNode*2 - 1), 0.5*dNdx_q(1, iNode*2); 0.5*dNdx_q(1, iNode*2), 0.0];
            BB(1:2, 1:2, iNode + 6) =[0.0, 0.5*dNdy_q(2,iNode*2 -1); 0.5*dNdy_q(2, iNode*2-1), dNdy_q(2,iNode*2)];
            
            Fe(iNode,1) = Fe(iNode, 1) + w(i) * be'*N1x_q(:, iNode);
            Fe(iNode+6,1) = Fe(iNode+6, 1) + w(i) * be'* N1y_q(:, iNode);
        end
        
        for iDof = 1:dofNum
            DD(iDof) = sum(diag(BB(:, :, iDof)));
            SS(1:2, 1:2, iDof) = 2*mue*BB(:,:,iDof) + lambdae * DD(iDof)* eye(2);
        end
        
        % Build Ke and Fe
        for a = 1:dofNum
            for b = 1:dofNum
%                 if (a == 5)
%                     disp(sum(sum(SS(:, :, b).*BB(:,:, a))))
%                 end
                Ke(a,b) = Ke(a,b) + w(i)*sum(sum(SS(:, :, b).*BB(:,:, a)));
            end
        end
    end
% <<-------------------------------------------------------------------->>
end

function [dN,N]=P2elementsX1(xe,dLdx,A2)
    % dN is gradN_tilda [1,6]
    % N is P2 function [2,6]
    lambda1 = @(x)((xe(2,2)-xe(2,3))*(x(1)-xe(1,2))+(xe(1,3)-xe(1,2))*(x(2)-xe(2,2))...
          )/A2;
    lambda2 = @(x)((xe(2,3)-xe(2,1))*(x(1)-xe(1,3))+(xe(1,1)-xe(1,3))*(x(2)-xe(2,3))...
          )/A2;
    lambda3 = @(x)((xe(2,1)-xe(2,2))*(x(1)-xe(1,1))+(xe(1,2)-xe(1,1))*(x(2)-xe(2,1))...
          )/A2;
    
    dN1=@(x) (1/A2).*[4.*lambda1(x)-1,0,0;0,0,0]*dLdx';
    dN2=@(x) (1/A2).*[0,4.*lambda2(x)-1,0;0,0,0]*dLdx';
    dN3=@(x) (1/A2).*[0,0,4.*lambda3(x)-1;0,0,0]*dLdx';

    dN4=@(x) (1/A2).*[4.*lambda2(x),4.*lambda1(x),0;0,0,0]*dLdx';
    dN5=@(x) (1/A2).*[0,4.*lambda3(x),4.*lambda2(x);0,0,0]*dLdx';
    dN6=@(x) (1/A2).*[4.*lambda3(x),0,4.*lambda1(x);0,0,0]*dLdx';

    dN=@(x) [dN1(x), dN2(x), dN3(x), dN4(x), dN5(x), dN6(x)];
    N=@(x) [2*lambda1(x).*(lambda1(x)-0.5), 2*lambda2(x).*(lambda2(x)-0.5), 2*lambda3(x).*(lambda3(x)-0.5),...
            4*lambda1(x).*lambda2(x), 4*lambda3(x).*lambda2(x), 4*lambda1(x).*lambda3(x);
            0,0,0,0,0,0];
end

function [dN,N]=P2elementsX2(xe,dNdL,A2)
    % dN is gradN_tilda [1,6]
    % N is P2 function [2,6]
    lambda1 = @(x)((xe(2,2)-xe(2,3))*(x(1)-xe(1,2))+(xe(1,3)-xe(1,2))*(x(2)-xe(2,2))...
          )/A2;
    lambda2 = @(x)((xe(2,3)-xe(2,1))*(x(1)-xe(1,3))+(xe(1,1)-xe(1,3))*(x(2)-xe(2,3))...
          )/A2;
    lambda3 = @(x)((xe(2,1)-xe(2,2))*(x(1)-xe(1,1))+(xe(1,2)-xe(1,1))*(x(2)-xe(2,1))...
          )/A2;
    
    dN1=@(x) (1/A2).*[0,0,0;4.*lambda1(x)-1,0,0]*dNdL';
    dN2=@(x) (1/A2).*[0,0,0;0,4.*lambda2(x)-1,0]*dNdL';
    dN3=@(x) (1/A2).*[0,0,0;0,0,4.*lambda3(x)-1]*dNdL';

    dN4=@(x) (1/A2).*[0,0,0;4.*lambda2(x),4.*lambda1(x),0]*dNdL';
    dN5=@(x) (1/A2).*[0,0,0;0,4.*lambda3(x),4.*lambda2(x)]*dNdL';
    dN6=@(x) (1/A2).*[0,0,0;4.*lambda3(x),0,4.*lambda1(x)]*dNdL';

    dN=@(x) [dN1(x), dN2(x), dN3(x), dN4(x), dN5(x), dN6(x)];
    N=@(x) [0,0,0,0,0,0;
            2*lambda1(x).*(lambda1(x)-0.5), 2*lambda2(x).*(lambda2(x)-0.5), 2*lambda3(x).*(lambda3(x)-0.5),...
            4*lambda1(x).*lambda2(x), 4*lambda3(x).*lambda2(x), 4*lambda1(x).*lambda3(x)];
end

%% Stresses and Strains
function [Strains, Stresses]=computeStrainAndStress(xe,ue,Ee,nue)
% Returns the strains and stresses in the element at the corners
% A lot of repeated code, only done for simplicity of explanation.
% <<-------------------------------------------------------------------->>
% Complete by computing Strains and Stresses 
% input: ue [12, 1]
    nodeNum = 6;
    dofNum = 12;
    Ke = zeros(dofNum, dofNum);
    Fe = zeros(dofNum, 1);
    
    mue=Ee/(1+nue)/2; %
    lambdae=mue*nue/(1-2*nue); %
    
    dLdx=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
          xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
    A2 = ((xe(2,2)-xe(2,3))*(xe(1,1)-xe(1,2))+(xe(1,3)-xe(1,2))*(xe(2,1)-xe(2,2)));
    
    nq = 3;
    [dNdx,N1x]=P2elementsX1(xe,dLdx,A2); % dNdx [2, 12]; N1x [2, 6]
    [dNdy,N1y]=P2elementsX2(xe,dLdx,A2);
    
    BB = zeros(2,2,dofNum);
    DD = zeros(dofNum);
    SS = zeros(2,2,dofNum);
    Strains=zeros(2,2,3);
    Stresses=zeros(2,2,3);
    
    for i = 1:nq
        xq = xe(:, i); % get the coordinates of quadrature points
        % Calculate the shape functions and their derivatives at the quadrature point.

        dNdx_q = dNdx(xq);
        dNdy_q = dNdy(xq);
        
        for iNode = 1:nodeNum
            BB(1:2, 1:2, iNode) = [dNdx_q(1, iNode*2 - 1), 0.5*dNdx_q(1, iNode*2); 0.5*dNdx_q(1, iNode*2), 0.0];
            BB(1:2, 1:2, iNode + 6) =[0.0, 0.5*dNdy_q(2,iNode*2 -1); 0.5*dNdy_q(2, iNode*2-1), dNdy_q(2,iNode*2)]; 
            
        end
        
        for iDof = 1:dofNum
            DD(iDof) = sum(diag(BB(:, :, iDof)));
            SS(1:2, 1:2, iDof) = 2*mue*BB(:,:,iDof) + lambdae * DD(iDof)* eye(2);
            Strains(:, :, i) = Strains(:, :, i) + BB(:,:,iDof)*ue(iDof);
            Stresses(:, :, i) = Stresses(:, :, i) + SS(:, :, iDof)*ue(iDof);
        end
    end
    Stresses = Stresses ./ 1e6; % convert to MPa
% <<-------------------------------------------------------------------->>
end