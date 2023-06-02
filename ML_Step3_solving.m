function ML_Step3_solving(matrixname,nb_experiment)

% This function injects errors randomly to a specific matrix.
% matrixnname: name of the matrix 
% bitflip_iter: iteration number to inject error 
% bitflip_element: element to inject error


M = nb_experiment; % number of experiments, each injects error at a random vector element

% load matrix
matrixfile = ['./matrices/', matrixname, '.mat'];
load(matrixfile);
A = Problem.A;
disp('Done loading matrix');
drawnow('update');
[N, ~] = size(A);
result_filename = ['./data/ML/Step3_', 'test', '.xls'];

% load preconditioner of matrix
precond_filename = ['./matrices/', matrixname, '_precond.mat'];
load(precond_filename);
disp('Done loading incomplete Cholesky factorization');
drawnow('update');

% load row 2-norm of matrix
norms_filename = ['./matrices/', matrixname, '_norms.mat'];
load(norms_filename);
disp('Done loading row 2-norms of matrix');
drawnow('update');

% setup
xx = ones(N, 1);    % all-1 vector
b = A*xx;           % b is set as A times the all-1 vector 
tol = 1e-6;
max_iter = 10000;

%% start pcg 
for m = 0:M
    if m == 0
        % Experiment 0: make an error-free run
        inject_error = 0;
        [~,flag,iter,~] = pcg4(A, b, tol, max_iter, L, L', inject_error, 0, 0);
        
        if flag == 1
           disp('error-free execution does not converge');
           return;
        end
        noerror_converge = iter;   % number of iterations in error-free run
        error_max_iter = noerror_converge*100;   % set max number of iterations to run when injecting errors (100x)
        disp(['Matrix = ', matrixname, ', Experiment=', num2str(m), ', converge=', num2str(noerror_converge)]);
        
        % create error file 
        error_filename = ['./matrices/ML/', matrixname, '_errorML.mat'];
        if exist(error_filename, 'file') < 1 % file does not exist 
            if noerror_converge <= M % if number of iteration in error-free run is <= number of experiments, then inject error in all elements
                E1 = [1:noerror_converge]'; % error in every iteration
                M = noerror_converge;
                if N > M
                E2 = datasample([1:N], M);
                else
                E2 = [1:N]'; 
                M = N;
                E1 = datasample([1:noerror_converge],M);
                end
            else
                E1 = datasample([1:noerror_converge], M);
                E2 = datasample([1:N], M);
            end
            E= cat(2,E1.',E2.');
            save(error_filename, 'E');
        else % file already exists, then just load the file 
            load(error_filename, 'E');
            M = length(E);
        end
    else
        % Inject errors from Experiment 1 to M, each at a random iteration 
        inject_error = 1;
        bitflip_iter = E(m,1);
        bitflip_element = E(m,2);
        
        [~,flag,iter,diff_v,first_abs_gradient,first_rel_gradient,xval] = pcg4(A, b, tol, error_max_iter, L, L', inject_error, bitflip_element, bitflip_iter);
        converge = iter;   % number of iterations in error-injecting run  

        grad_abs(:, m) = first_abs_gradient;
        grad_rel(:, m) = first_rel_gradient;

        result = [N,flag,bitflip_iter,bitflip_element,diff_v,A_row_2norm(bitflip_element),grad_abs(bitflip_element),grad_rel(bitflip_element),xval(bitflip_element),converge/noerror_converge];
        writematrix(result,result_filename,'WriteMode','append');
        
        disp(['Matrix = ', matrixname, ', Experiment=', num2str(m), ', converge=', num2str(converge)]);
        if flag == 1
            disp('did not converged');
        end
        
    end
end

end