function [x,flag,iter,diff_v,first_temp_gradient,second_temp_gradient] = pcg4(A,b,tol,maxit,M1,M2, inject_error,bitflip_pos,bitflip_iter)

flag = 0;
N = length(A);
x = zeros(N, 1);
x1 = x;
x2 = x;

r = b - A*x;
z = M2\(M1\r);
p = z;

iter = 0;
normr = norm(r);
normb = norm(b);
relres = normr/normb;

while (iter < maxit) && (relres > tol)
    %inject error in vector p
    if (iter == bitflip_iter) && (inject_error == 1)
        diff_v = max_p;
        p(bitflip_pos) = p(bitflip_pos) + diff_v;
    end
    
    q = A*p;
    v = z'*r;
    alpha = (r'*z)/(p'*q);
    x2 = x1;
    x1 = x;
    x = x + alpha * p;         
    
    % look at properties of vector x
    if (iter == bitflip_iter - 1) && (inject_error == 1)
        first_temp_gradient = abs(x(bitflip_pos)- x1(bitflip_pos));
        second_temp_gradient = abs(x(bitflip_pos)- 2*x1(bitflip_pos) + x2(bitflip_pos));
    end
    
    r = r - alpha * q;
    z = M2\(M1\r);
    beta = z'*r/v;
    p = z + beta*p;
    
    normr = norm(r);
    relres = normr/normb;
    
    iter = iter + 1;
end

if relres <= tol
    flag = 0; % converged
elseif iter == maxit
    flag = 1; % no converge and maxit reached
end

end