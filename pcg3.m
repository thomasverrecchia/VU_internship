function [x, flag, iter, diff_v] = pcg3(A, b, tol, maxit, M1, M2, inject_error, bitflip_pos, bitflip_iter)

diff_v = 0;
flag = 0;           % status of the final result 
N = length(A);
x = zeros(N, 1);    % initial guess of solution - all-0 vector
r = b - A*x;
z = M2\(M1\r);      % M = M1M2 is the preconditioner matrix 
p = z;

iter = 0;
normr = norm(r);
normb = norm(b);
relres = normr/normb;

while (iter < maxit) && (relres > tol)
    %inject error in vector p
    if (iter == bitflip_iter) && (inject_error == 1)
        diff_v = max(p); 
        p(bitflip_pos) = p(bitflip_pos) + diff_v;
    end
    
    q = A*p;
    v = z'*r;
    alpha = v/(p'*q);
    x = x + alpha * p;         
    r = r - alpha * q;
    z = M2\(M1\r);
    beta = (z'*r)/v;
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

end % end of function 