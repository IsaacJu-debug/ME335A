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
nel_array = [10, 100, 500]; % input array for element number
lambda_array = [-10, 0, 10]; % input array for lambda
fig_count1 = 1; % counter for function and its derivative
output_folder = 'plots';
u_error = zeros(9,1);
u_prime_error = zeros(9,1);
error_count = 1; % counter for recording error

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

for k=1:length(lambda_array)
    subplot_num = 1;
    fig = figure;
    fig.Name = 'Function and Derivative';

    for n=1:length(nel_array)    
        %nel=100;
        nel = nel_array(n);
        lambda = lambda_array(k);
        
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
        % Create a vector of element indices (1 to n_elements)
        element_indices = 1:nel;
        % Create the local-to-global mapping using vectorized operations
        LG = element_indices + (0:(kk - 1))';
        % Reshape the LG array to the desired shape [kk, nel]
        LG = reshape(LG, [kk, numel(element_indices)]);

        % <<-------------------------->>  

        %% finite element solver begins
        %%

        %% parameters of the problem
        EtaG=[nunk];
        GG=[1];
        hh=zeros(1,nod);
        hh(1)=-1;
        ff=-2*ones(1,nel);
        

        %% assembly, from local to global
        K=zeros(nunk,nunk);
        F=zeros(nunk,1);

        %% Assembly of active indices
        for iel=1:nel
          % setting the local data
          % loop through elements
          lge=LG(:,iel);
          xe(1:dd,1:npe)=X(1:dd,iel:iel+1);
          fe=ff(iel);
          hhe = hh(iel:iel+1);

          % computing element K and F
          [Ke, Fe]=elementKandF(xe,fe,lambda,hhe);
          % <<-------------------------->>
          % Complete with the assembly of Ke and Fe into K and F 

          K( lge(1):lge(2), lge(1):lge(2)) =  K( lge(1):lge(2), lge(1):lge(2)) + Ke(:, :);

          F(lge(1):lge(2)) = F(lge(1):lge(2)) + Fe(:);

          % <<-------------------------->>  
        end

        %% Assembly of constrained indices
        ng=length(EtaG);
        for ig=1:ng
          % <<-------------------------->>
          % Complete with the assembly of K and F 

          K(EtaG(ig), :) = 0;
          K(EtaG(ig), EtaG(ig)) = 1;
          F(EtaG(ig)) = 1; % g = 1
          % <<-------------------------->>  
        end

        %% Solve 
        U=K\F;

    
        %% plotting (mesh specific)
        % plotting points (uniform spacing dx)
        %clf(fig1);
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
        
        % Exact solution
        uep=zeros(1,nxp);
        duep=zeros(1,nxp);
        for kp=1:nxp
          % <<-------------------------->>
          % Complete uep and duep with the exact solution at each sampling point
          if lambda<0
              sq_lambda = sqrt(-lambda);
              A = 1/sq_lambda; C = -2/lambda;
              B = (lambda + sq_lambda*sin(sq_lambda) + 2)/(lambda*cos(sq_lambda));
              uep(kp) = A*sin(sq_lambda*xp(kp)) + C +  B*cos(sq_lambda*xp(kp)) ;
              duep(kp) = sq_lambda*(A*cos(sq_lambda*xp(kp)) - B*sin(sq_lambda*xp(kp)));
          elseif lambda==0
              uep(kp) = xp(kp)^2 + xp(kp) - 1;
              duep(kp) = 2*xp(kp) + 1;
          else
              sq_lambda = sqrt(lambda);
              C = -2/lambda;
              A = (lambda + 2 + sq_lambda*exp(-sq_lambda))/(lambda*(exp(-sq_lambda) + exp(sq_lambda)));
              B = (lambda + 2 - sq_lambda*exp(sq_lambda))/(lambda*(exp(-sq_lambda) + exp(sq_lambda)));
              uep(kp) = A*exp(sq_lambda*xp(kp)) +  C + B*exp(-sq_lambda*xp(kp));
              duep(kp) = sq_lambda*(A*exp(sq_lambda*xp(kp)) - B*exp(-sq_lambda*xp(kp)));
          end
          % <<-------------------------->>  
        end
        
        % First plot
        subplot(2, 3, n);
        plot(X, U(1:nunk), "or", "linewidth", 2);
        title(['$u(x)$, nel: ', num2str(nel), ', $\lambda$: ', num2str(lambda)], 'Interpreter', 'latex');
        set(gca, "FontName", "Arial");
        set(gca, "FontSize", 16);
        hold on;
        plot(xp, up, "r-", "linewidth", 2);
        plot(xp, uep, "b-", "linewidth", 2);
        hold off;

        % Second plot
        subplot(2, 3, n + 3);
        plot(xp, dup, "r-", "linewidth", 2);
        set(gca, "FontName", "Arial");
        set(gca, "FontSize", 16);
        title(['$\frac{du(x)}{dx}$, nel: ', num2str(nel), ', $\lambda$: ', num2str(lambda)], 'Interpreter', 'latex');
        hold on;
        plot(xp, duep, "b-", "linewidth", 2);
        hold off;
        
        % save error
        u_error(error_count) = norm(up - uep);
        u_prime_error(error_count) = norm(dup - duep);
        error_count = error_count + 1;
    end
    
    % Save the figure with both subplots as a PNG file
    output_filename = sprintf('Function_and_Derivative_%d.png', fig_count1);
    full_output_path = fullfile(output_folder, output_filename);
    saveas(gcf, full_output_path);

    % Update fig_count1
    fig_count1 = fig_count1 + 1;
    
    % display error table
    disptable(lambda, nel_array, u_error((k - 1)*3 + 1: (k)*3),...
                u_prime_error((k - 1)*3 + 1: (k)*3));
    
end

%% error plots
close all
nel_2D_array = repmat(nel_array, 3, 1);

full_output_path = fullfile(output_folder, 'error_u.png');
loglog_error_plot(reshape(u_error, [3, 3]), nel_2D_array, 1, ...
                 full_output_path,lambda_array );

full_output_path = fullfile(output_folder, 'error_u_prime.png');
loglog_error_plot(reshape(u_prime_error, [3, 3]), nel_2D_array, 2,...
                full_output_path, lambda_array);

%% element matrix and load
function [Ke, Fe]=elementKandF(xe,fe,lambdae,hhe)
  h=xe(2)-xe(1);  % Element size
  % <<-------------------------->>
  % Complete with the values of Ke and Fe
  % local stiffness matrix
  Ke = zeros(2,2);
  Fe = zeros(2,1);
  Ke(1,1) = (-1/h ) * (-1/h) * h + lambdae * h/3; % 1st order shape functions has constant derivative
  Ke(1,2) = (-1/h ) * (1/h) * h +lambdae*h/6;
  Ke(2,1) = Ke(1,2); % symmetry
  Ke(2,2) = (1/h ) * (1/h) * h + lambdae * h/3;
  
  Fe(1) = fe * h /2 + hhe(1);
  Fe(2) = fe * h /2 + hhe(2);
  % <<-------------------------->>  
end


%% basis functions (only for plotting)
function [N1, N2, dN1, dN2]=P1(x)
  % <<-------------------------->>
  % Complete with the values of N1(x), N1'(x), N2(x), and N2'(x) over an interval [0,1]
  N1 = 1-x; 
  N2 = x; 
  dN1 = -1; 
  dN2 = 1;
  % <<-------------------------->> 
end

%% display table
function disptable( lambda, nel_array, error_u, error_u_prime)
    % Compute the 'h' array
    disp('Table: Error')
    
    h = 1./ nel_array;
    
    % Combine the arrays into a single matrix
    lambda = ones(size(nel_array)) * lambda;
    data = [lambda(:), nel_array(:), h(:), error_u(:), error_u_prime(:)];

    % Define the column names
    colname = ["Lambda", "nel", "h", "||u-u_h||", "||u'-u_h'||"];

    % Display the table
    disp(array2table(data, 'VariableNames', colname))
end

%% log-log error plot
function loglog_error_plot(error_array, h_array, opt, full_output_path, lambda_array)
    figure;
    hold on;

    % Define colors for lines and markers
    colors = {'r', 'g', 'b'};

    for i = 1:3
        % Extract error and h data from each column
        error_data = error_array(:, i);
        h_data = h_array(i, :);

        % Plot the scatter points for each row using loglog
        loglog(h_data, error_data, 'o', 'MarkerEdgeColor', colors{i}, 'MarkerFaceColor', colors{i});

        % Fit a linear line to the log-log data
        p = polyfit(log(h_data), log(error_data), 1);

        % Calculate the fitted line
        fitted_line = h_data .^ p(1) * exp(p(2));

        % Plot the fitted line using loglog
        loglog(h_data, fitted_line, colors{i}, 'LineWidth', 1.5);
    end

    grid on;

    % Set labels and title
    xlabel('log(1/h)');
    ylabel('Error');

    if (opt == 1)
        title(('$u(x)$, Nel vs. Error '), 'Interpreter', 'latex');
    else
        title(('$\frac{du(x)}{dx}$, Nel vs. Error'), 'Interpreter', 'latex');
    end

    legend_text = strings(1, 6);
    for i = 1:3
        legend_text(2*i-1) = ['Error ', num2str(i), '. $\lambda = ', num2str(lambda_array(i)), '$'];
        legend_text(2*i) = ['Fit ', num2str(i)];
    end
    lgd = legend(legend_text, 'Location', 'Best', 'Interpreter', 'latex');
    lgd.FontSize = 14;

    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    saveas(gcf, full_output_path);
end
