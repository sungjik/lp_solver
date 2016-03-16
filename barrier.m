function x = barrier(A, c, x0, phase_one)
[m,n] = size(A);
x = x0;
u = 50;
newton_steps = [];
duality_gaps = [];
num_centering_itrs=0;
k=1;
while n/k >= 1e-3
    [x,itrs] = newton_eq(A,c,k,x,phase_one);
    num_centering_itrs = num_centering_itrs + 1;
    newton_steps(num_centering_itrs) = itrs;
    duality_gaps(num_centering_itrs) = n/k;
    k=k*u;
    %if phase_one == true, find a strictly feasible point and exit
    %immediately for phase 2
    if phase_one && x(end) < 1
        break
    end
end
close all;
figure;
[xx, yy] = stairs(cumsum(newton_steps),duality_gaps);
semilogy(xx,yy);
end