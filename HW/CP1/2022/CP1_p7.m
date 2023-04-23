display("Second-order problems with P1 finite elements")
%% Simple FEA program for ME335-Finite Element Analysis
%% Authors: Gustavo Buscaglia and Adrian Lew
%% Names of variables (in order of appearance)
% nel:      number of elements
% nod:      number of vertices in the mesh (nodes in this case)
% X:        coordinates of vertices
% npe:      number of nodes per element 
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
nel_array = [10, 100, 500];
lam_array = [-10, 0, 10];
err_matrx = zeros(3, 3);
derr_matrx = zeros(3, 3);
i_fig1 = 1; i_fig2 = 2;
for n_1=1:length(nel_array)
    for k_1=1:length(lam_array)
        %% Build a mesh  -- Elements are implicitly defined by two consecutive vertices
        nel=nel_array(n_1);
        nod=nel+1;
        X=[0:1:nel]/nel;
        
        %% Useful Constants 
        npe=2; 
        dd=1; 
        L=1;
        
        %% Build local-to-global map
        nunk=nod;
        kk=2;
        LG=zeros(npe, nel);
        % <<-------------------------->>
        % Complete with the value of LG
        for ig=1:nel
            LG(:,ig)=[ig; ig+1];
        end
        % <<-------------------------->>  
          
        %% finite element solver begins
        %%
        
        %% parameters of the problem
        EtaG=[nunk];
        GG=[1];
        hh=zeros(1,nod);
        hh(1)=-1;
        ff=-2*ones(1,nel);
        lambda=lam_array(k_1);
        
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
          K(lge(1), lge(1)) = K(lge(1), lge(1)) + Ke(1,:);
          K(lge(1), lge(2)) = K(lge(1), lge(2)) + Ke(2,:);
          if iel<nel
              K(lge(2), lge(2)) = K(lge(2), lge(2)) + Ke(1,:);
              K(lge(2), lge(1)) = K(lge(2), lge(1)) + Ke(2,:);
          end
        
          for ig=1:npe
              if lge(ig)<nel+1
                  F(lge(ig),:) = F(lge(ig),:) + Fe(ig,:);
              end
          end
          
          % <<-------------------------->>  
        end
          
        %% Assembly of constrained indices
        ng=length(EtaG);
        for ig=1:ng
          % <<-------------------------->>
          % Complete with the assembly of K and F 
          K(EtaG(ig), EtaG(ig)) = 1;
          F(EtaG(ig),:) = 1;
          % <<-------------------------->>  
        end
        %% Solve 
        U=K\F;
        
        
        %% plotting (mesh specific)
        % plotting points (uniform spacing dx)
        fig1=figure(i_fig1);
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
         %disp(nle)%, disp(U(nle))
         [N1, N2, dN1, dN2]=P1(xl); %, X(nle), X(nle+1)
         % Please pay attention to how the values of up and dup are computed
         up(kp)=U(nle)*N1+U(nle+1)*N2;
         dup(kp)=(U(nle)*dN1+U(nle+1)*dN2)/(X(nri)-X(nle));
        end
        plot(X,U(1:nunk),"or","linewidth",2);
        title(['Function, Elements: ',num2str(nel),', Lambda: ',num2str(lambda)])
        set(gca,"FontName","Arial")
        set(gca,"FontSize",16)
        hold on
        plot(xp,up,"r-","linewidth",2);
        hold off
        fig2=figure(i_fig2);
        fig2.Name='Derivative';
        clf(fig2);
        plot(xp,dup,"r-","linewidth",2);
        title(['Derivative, Elements: ',num2str(nel),', Lambda: ',num2str(lambda)])
        % Exact solution
        uep=zeros(1,nxp);
        duep=zeros(1,nxp);
        for kp=1:nxp
          % <<-------------------------->>
          % Complete uep and duep with the exact solution at each sampling point
          if lambda==0
              uep(kp) = xp(kp)^2 + xp(kp) - 1;
              duep(kp) = 2*xp(kp) + 1;
          elseif lambda<0
              sq_lambda = sqrt(-lambda);
              A = 1/sq_lambda; C = -2/lambda;
              B = (lambda + sq_lambda*sin(sq_lambda) + 2)/(lambda*cos(sq_lambda));
              uep(kp) = A*sin(sq_lambda*xp(kp)) + B*cos(sq_lambda*xp(kp)) + C;
              duep(kp) = sq_lambda*(A*cos(sq_lambda*xp(kp)) - B*sin(sq_lambda*xp(kp)));
          else
              sq_lambda = sqrt(lambda);
              C = -2/lambda;
              A = (lambda + 2 + sq_lambda*exp(-sq_lambda))/(lambda*(exp(-sq_lambda) + exp(sq_lambda)));
              B = (lambda + 2 - sq_lambda*exp(sq_lambda))/(lambda*(exp(-sq_lambda) + exp(sq_lambda)));
              uep(kp) = A*exp(sq_lambda*xp(kp)) + B*exp(-sq_lambda*xp(kp)) + C;
              duep(kp) = sq_lambda*(A*exp(sq_lambda*xp(kp)) - B*exp(-sq_lambda*xp(kp)));
          end
          % <<-------------------------->>  
        end
        figure(i_fig1);
        hold on
        plot(xp,uep,"b-","linewidth",2);
        hold off
        figure(i_fig2)
        hold on
        plot(xp,duep,"b-","linewidth",2);
        hold off
        i_fig1 = i_fig1 + 2;
        i_fig2 = i_fig2 + 2;

        err_matrx(n_1, k_1) = norm(uep - up);
        derr_matrx(n_1, k_1) = norm(duep - dup);
    end
end
colors = ["b", "r", "g"];
for i=1:3
    x = log(nel_array);
    y = log(err_matrx(:, i));
    c = polyfit(x,y,1);
    figure(i_fig2+1);
    hold on
    plot(x,y,colors(i),'linestyle','none','marker','o', 'DisplayName',strcat('rate=',num2str(abs(c(1))),'  lambda=',num2str(lam_array(i))));
    y_est = polyval(c,x);
    % Add trend line to plot
    hold on
    plot(x,y_est,colors(i),'LineWidth',2,'DisplayName',strcat('Fitting - ', num2str(i)))
    title(['Function Error Norm'])
    hold off
end
legend('show')

for i=1:3
    y = log(derr_matrx(:, i));
    c = polyfit(x,y,1);
    figure(i_fig2+2)
    hold on
    plot(x,y,colors(i),'linestyle','none','marker','o', 'DisplayName',strcat('rate=',num2str(abs(c(1))),'  lambda=',num2str(lam_array(i))));
    y_est = polyval(c,x);
    % Add trend line to plot
    hold on
    plot(x,y_est,colors(i),'LineWidth',2,'DisplayName',strcat('Fitting - ', num2str(i)))
    title(['Derivative Error Norm'])
    hold off
end
legend('show')


%% element matrix and load
function [Ke, Fe]=elementKandF(xe,fe,lambdae,hhe)
  h=xe(2)-xe(1);  % Element size
  % <<-------------------------->>
  % Complete with the values of Ke and Fe
  Ke = zeros(2, 1);
  Fe = zeros(2, 1);
  Ke(1,:) = (1/h)+lambdae*h/3;
  Ke(2,:) = -(1/h)+lambdae*h/6;
  % hhe(1,:), fe*h/2 , Fe(1,:)
  Fe(1,:) = fe*h/2 + hhe(1); %
  Fe(2,:) = fe*h/2 + hhe(2);
  % <<-------------------------->>  
end


%% basis functions (only for plotting)
function [N1, N2, dN1, dN2]=P1(x) %, i, j
  % <<-------------------------->>
  % Complete with the values of N1(x), N1'(x), N2(x), and N2'(x) over an interval [0,1]
  N1 = 1-x; %(x-j)/(i-j);
  N2 = x; %(x-i)/(j-i);
  dN1 = -1; dN2 = 1;
  %dN1 = 1/(i-j); dN2 = 1/(j-i);
  % <<-------------------------->> 
end



