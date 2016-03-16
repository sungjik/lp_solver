function [x, num_descent_itrs] = newton_eq(A, c, k, x0, phase_one)
alpha = 0.1;
beta = 0.9;
h = 1e-6;
x = x0;
lambdas = [];
objectives = [];
num_descent_itrs = 0;

while true
    %compute newton step with gradient and Hessian
    grad_fx = k*c - 1./x;
    H = diag(1./x.^2);
    H_inv = diag(x.^2);
    w = inv(-A*H_inv*A')*A*H_inv*grad_fx;
    x_ns = -H_inv*(grad_fx + A'*w);
    
    %stopping criterion
    lambda_sq = -x_ns'*grad_fx;
    if lambda_sq/2 <= h
        lambdas(num_descent_itrs+1) = lambda_sq/2;
        break;
    end
    %early stopping for phase one (can stop if a strictly feasible point is
    %found)
    if phase_one && x(end) < 1
        break;
    end
    
    %backtracking line search
    t = 1;
    num_backtracking_itrs = 0;
    num_domain_itrs = 0;
    
    %find x_step in domain of f
    while ~all(pos(x + t*x_ns))
        t = t*beta;
        num_domain_itrs = num_domain_itrs + 1;
    end
        
    while k*c'*(x + t*x_ns) - sum(reallog(x + t*x_ns)) > k*c'*x - sum(reallog(x)) + alpha*t*grad_fx'*x_ns
        t = beta*t;
        num_backtracking_itrs = num_backtracking_itrs + 1;
    end
    
    %update x
    x = x + t*x_ns;
    num_descent_itrs = num_descent_itrs + 1;
    lambdas(num_descent_itrs) = lambda_sq/2;
    objectives(num_descent_itrs) = k*c'*x - sum(reallog(x));
end

num_descent_itrs;
% figure;
% hold on;
% plot([1:num_descent_itrs+1],lambdas, 'r');
% plot([1:num_descent_itrs], objectives, 'b');
