close all;
clear all;

%generate problem data that results in well-conditioned matrices
n = 500;
m = 100;
rng(1,'twister');
A = randn(m,n);
A = A + [ones(1,n); zeros(m-1,n)];
x0 = abs(randn(n,1));
b = A*x0;
c = randn(n,1)+m*ones(n,1);

%solve using cvx library
cvx_begin
cvx_quiet(true)
    variable x(n)
    minimize( c'*x )
    subject to
        A*x == b
        x >= 0
cvx_end
cvx_optval

%solve using lp_solver function
x_lp = lp_solver(A, b, c);
c'*x_lp