%% HW3(a) 
syms u(x)
syms lambda
Du = diff(u);

ode = - diff(u,x,2) + lambda * diff(u,x,1) == x;
cond1 = u(0) == 0.1;
cond2 = Du(0) == 1;

conds = [cond1 cond2];
uSol(x) = dsolve(ode,conds);
uSol
%uSol = simplify(uSol)

%% HW3(c)
syms u(x)
syms lambda
syms v(x)
syms u_prime
syms v_prime
syms l_inte
syms c1
syms c2
%c1 = 0.1;
%c2 = (lambda^2 - lambda - 1)/(lambda^3* exp(lambda));
uSol = c1 + c2 * exp(lambda * x) + x^2/(2*lambda) + x/(lambda^2);
v(x) = 2*x^3 - 2*x^2;
u_prime = diff(uSol, x);
v_prime = diff(v, x);
l_inte_part1 = u_prime * v_prime;
l_inte_part2 = lambda * u_prime * v;
lhs_part1 = int(l_inte_part1, 0, 1);
lhs_part1
lhs_part2 = int(l_inte_part2, 0, 1);
lhs_part2
lhs = lhs_part1 + lhs_part2;
lhs
rhs = int(x * v(x) , 0, 1)+ v(1);
rhs

%% HW3(d)
syms u(x)
syms lambda
syms v(x)
syms u_prime
syms v_prime
syms l_inte
syms c1
syms c2
%c1 = 0.1;
%c2 = (lambda^2 - lambda - 1)/(lambda^3* exp(lambda));
uSol = c1 + c2 * exp(lambda * x) + x^2/(2*lambda) + x/(lambda^2);
v(x) = 3*x^3 - x^2;
u_prime = diff(uSol, x);
v_prime = diff(v, x);
l_inte_part1 = u_prime * v_prime;
l_inte_part2 = lambda * u_prime * v;
lhs_part1 = int(l_inte_part1, 0, 1);
lhs_part1
lhs_part2 = int(l_inte_part2, 0, 1);
lhs_part2
lhs = lhs_part1 + lhs_part2;
lhs

rhs = int(x * v(x) , 0, 1) + v(1)*1.0;
rhs

%% HW3(d.2)
% back solve c2
syms c2_d
eqn = lhs == rhs;
sol_c2 = solve(eqn, c2);
c2 = sol_c2;
c2_d = c2;
new_uSol = c1 + c2 * exp(lambda * x) + x^2/(2*lambda) + x/(lambda^2);
new_uSol
new_d_uSol = diff(new_uSol, x);
simplify(subs(new_d_uSol, x, 1))

%% HW3(e)
syms u(x)
syms lambda
syms v(x)
syms u_prime
syms v_prime
syms l_inte
syms c1
syms c2
%c1 = 0.1;
%c2 = (lambda^2 - lambda - 1)/(lambda^3* exp(lambda));
c2 = c2_d;
uSol = c1 + c2 * exp(lambda * x) + x^2/(2*lambda) + x/(lambda^2);
v(x) = x;
u_prime = diff(uSol, x);
v_prime = diff(v, x);
l_inte_part1 = u_prime * v_prime;
l_inte_part2 = lambda * u_prime * v;
lhs_part1 = int(l_inte_part1, 0, 1);
lhs_part1
lhs_part2 = int(l_inte_part2, 0, 1);
lhs_part2
lhs = lhs_part1 + lhs_part2;
lhs
rhs = int(x * v(x) , 0, 1) + v(1) - 0.1 + subs(uSol, x, 0);
rhs

eqn = lhs == rhs;
sol_c1 = solve(eqn, c1);
simplify(sol_c1)

%% hw3(e)_2
c1 = sol_c1
new_uSol = c1 + c2 * exp(lambda * x) + x^2/(2*lambda) + x/(lambda^2);
simplify(subs(new_uSol, x, 0))


