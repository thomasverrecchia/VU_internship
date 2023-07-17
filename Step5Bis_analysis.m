
function Step5Bis_analysis

% This function analyzes the slowdown associated with our selective 
% protection algorithm and the random protection algorithm under 
% different fractions of protection (from 1% to 100%). 

comments = 'Step5';

%matrices = {'cvxbqp1', 'thermal1', 'nd6k', ...
%    'bcsstk18', 'bodyy5', 'cbuckle', 'Pres_Poisson', 'bcsstk36', 'ct20stif', 'gyro_m', 't2dah_e', 'm_t1', 'msc23052', '2cubes_sphere', 'pwtk', 'G2_circuit', 'raefsky4', ...
%    'Trefethen_20000', 'vanbody','wathen100'};
matrices = {'ct20stif'};
num_matrices = length(matrices);

for m = 1:num_matrices
               
    %% Store the predicted slowdown in a matrix

    % load matrix
    matrixname = matrices{m};
    matrixfile = ['./matrices/', matrixname, '.mat'];
    load(matrixfile);
    A = Problem.A;
    disp('Done loading matrix');
    drawnow('update');
    [N, ~] = size(A);
    
    % load preconditioner of matrix
    precond_filename = ['./matrices/', matrixname, '_precond.mat'];
    load(precond_filename);
    disp('Done loading incomplete Cholesky factorization');
    drawnow('update');
    
    % load row 2-norm of matrix
    norms_filename = ['./matrices/', matrixname, '_norms.mat'];
    load(norms_filename);
    N = length(A_row_2norm);
    [~, sorted_positions] = sort(A_row_2norm, 'descend');
    rand_positions = randperm(N);
    
    % setup
    xx = ones(N, 1);    % all-1 vector
    b = A*xx;           % b is set as A times the all-1 vector 
    tol = 1e-6;
    max_iter = 10000;

    inject_error = 0;
    [~,~,iter,~] = pcg4(A, b, tol, max_iter, L, L', inject_error, 0, 0);
    noerror_converge = iter;   % number of iterations in error-free run
    disp(matrixname);
    
    modelname = ['data/ML/TrainedModel/', matrixname, '_RandomForest_V2.mat'];
    trained_model = load(modelname);
 

    predictname = ['data/ML/PredictedSlowdown/prediction_',matrixname,'_V2.mat'];
    if exist(predictname, 'file') < 1
        iter = zeros(N*noerror_converge,1);
        norm = zeros(N*noerror_converge,1); 
        predicted_slowdown = zeros(N*noerror_converge,1);
        compteur = 1;
        for i = 1:N
            for j = 1:noerror_converge
                iter(compteur) = j;
                norm(compteur) = A_row_2norm(i);
                compteur = compteur + 1;
            end
        end
        parfor i = 1:size(iter,1)
            predicted_slowdown(i) = trained_model.predictFcn([iter(i),norm(i)]);
            disp([num2str(i)," out of ", num2str(size(iter,1))])
        end    
        save(predictname,'predicted_slowdown')
    else
        load(predictname);
    end

    %% Choosing the elements to protect
    
    protectname = ['data/ML/PredictedSlowdown/protected_',matrixname,'_V2.dat'];
    if exist(protectname, 'file') < 1
        protected = zeros(size(predicted_slowdown,1),2);
        compteur = 1;
        for i = 1:size(predicted_slowdown,1)
            disp(size(predicted_slowdown,1)\i)
            if predicted_slowdown(i) == 1
                element = fix(i/noerror_converge);
                iter = i-element*noerror_converge;
                protected(compteur,1) = iter;
                protected(compteur,2) = element;
                compteur = compteur + 1;
            end
        end

        protected_iter = protected(:,1);
        protected_element = protected(:,2);
        protected_iter = protected_iter(protected_iter ~=0);
        protected_element = protected_element(protected_element ~= 0);
        protected_element = protected_element(1:size(protected_iter,1),:);
        protected_V2 = cat(2,protected_iter,protected_element);
        writematrix(protected_V2,protectname)
        
    else
        protected = readmatrix(protectname);
    end
    num_protects = size(protected,1);

    %% load experimental data
    result_filename = ['./data/ML/Step3ML_', matrixname, '.xls'];
    result = xlsread(result_filename);

    error_iters = result(:, 3);
    error_elements = result(:,4);
    converge_ratios = result(:, 10);
    num_exps = length(error_iters);

    %% analyze data by ML
    protect_method = 'ML_Model'; 
    analysis_filename = ['./data/ML/Analysis/', comments, '_', matrixname, '_', protect_method, '_V2.dat'];
    slowdowns = zeros(num_exps,1);
    for p = 1:num_exps
        error_iter = error_iters(p);
        error_element = error_elements(p);
        converge_ratio = converge_ratios(p);
        if ismember([error_iter,error_element], protected)
            slowdowns(p) = 1;
        else
            slowdowns(p) = converge_ratio; 
        end
    end
    writematrix(slowdowns,analysis_filename);

    protects = [0:0.01:1];   % percentage of protection
    num_protects = length(protects);

    %% Analyse data by Row-2 Norm
    protect_method = 'Arow2norm'; 
    analysis_filename = ['./data/ML/Analysis/', comments, '_', matrixname, '_', protect_method, '.dat'];
    slowdowns = zeros(num_protects, num_exps);
    for p = 1:num_protects
        protect_percent = protects(p);
        protect_number = ceil(protect_percent * N);
        protect_positions = sorted_positions(1:protect_number);
        for e = 1:num_exps
            error_position = error_elements(e);
            converge_ratio = converge_ratios(e);
            if ismember(error_position, protect_positions)
                slowdowns(p, e) = 1;
            else
                slowdowns(p, e) = converge_ratio; 
            end
        end
    end
    writematrix(slowdowns,analysis_filename);
    
    %% analyze data by random
    protect_method = 'random'; 
    analysis_filename = ['./data/ML/Analysis/', comments, '_', matrixname, '_', protect_method, '.dat'];
    slowdowns = zeros(num_protects, num_exps);
    for p = 1:num_protects
        protect_percent = protects(p);
        protect_number = ceil(protect_percent * N);
        protect_positions = rand_positions(1:protect_number);
        for e = 1:num_exps
            error_position = error_elements(e);
            converge_ratio = converge_ratios(e);
            if ismember(error_position, protect_positions)
                slowdowns(p, e) = 1;
            else
                slowdowns(p, e) = converge_ratio; 
            end
        end
    end
    writematrix(slowdowns,analysis_filename);
end 


end

