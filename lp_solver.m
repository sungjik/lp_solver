function x = lp_solver(A, b, c)
[m,n] = size(A);
if m >= n || rank(A) ~= m
    error('input matrix A must be fat and full column rank');
end
x0 = A\b;
t0 = 2 - min(x0);
v = barrier([A -A*ones(n,1)], [zeros(n,1);1], [x0+(t0-1)*ones(n,1); t0], true);
z = v(1:n);
t = v(end);
x_feasible = z + (1-t)*ones(n,1);
x = barrier(A,c,x_feasible, false);
end