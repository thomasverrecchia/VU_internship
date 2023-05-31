function Step3Bis_matrices

% This function calls Step3Bis_solving, which injects error to a specific
% matrix and specific element. 

%matrices = {'cvxbqp1', 'thermal1', 'nd6k', ...
%    'bcsstk18', 'bodyy5', 'cbuckle', 'Pres_Poisson', 'bcsstk36', 'ct20stif', 'gyro_m', 't2dah_e', 'm_t1', 'msc23052', '2cubes_sphere', 'pwtk', 'G2_circuit', 'raefsky4', ...
%    'Trefethen_20000', 'vanbody','wathen100'};
matrices = {'bcsstk18'}; 
  
num_matrices = length(matrices);
% bitflip_element = 1;


for m = 1:num_matrices
    matrixname = matrices{m};
    matrixfile = ['./matrices/', matrixname, '.mat'];
    load(matrixfile);
    A = Problem.A;
    disp('Done loading matrix');
    drawnow('update');
    [N, ~] = size(A);
    K = datasample([1:N], 15);
    for k = 1:length(K)  
    Step3Bis_solving(matrixname, K(k));
    end
end 

end