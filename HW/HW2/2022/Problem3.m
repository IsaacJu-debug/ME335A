syms x
n_pts = 50;
x_vec = linspace(0, 1, n_pts);
an_sol = - (sin(2*pi*x_vec)/(2*pi)).^2 / 2;
figure
plot(x_vec, an_sol, '--', 'DisplayName','Analytical')
hold on
for ii = 1:3
    k = 2*ii + 1;
    k
    K_mat = zeros([k-1, k-1]);
    F_mat = zeros([k-1, 0]);
    for jj = 1:k-1
        K_mat(jj, jj) = int((Legendre(jj+1, x) * Legendre(jj+1, x)), 0, 1);
        F_mat(jj, 1) = integral_F(jj+1, x);
    end
    disp(F_mat)
    disp(K_mat)
    func_mat = zeros([n_pts, 1]);
    for jj = 1:n_pts
        for kk = 1:k-1
            
            func_mat(jj, 1) = func_mat(jj, 1) + F_mat(kk, 1)*iLegendre(kk+1, x_vec(jj));
        end
    end
    plot(x_vec,func_mat, 'DisplayName',strcat('n=',num2str(k)))
    hold on
end
hold off
legend
