function Step6Bis_overhead

% This function plots Figure 6, 8 (combined) and Figure 7 in the paper

comments = 'Step5';
close all;

%matrices = {'bcsstk18', 'bodyy5', 'cbuckle', 'G2_circuit'};
matrices = {'bodyy5','bcsstk18','cbuckle','G2_circuit'};
num_matrices = length(matrices);

bitflip_iter = 1;
protects = [0:0.01:1];
num_protects = length(protects);

for m = 1:num_matrices

    matrixname = matrices{m};
    disp(matrixname);

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
    
    % setup
    xx = ones(N, 1);    % all-1 vector
    b = A*xx;           % b is set as A times the all-1 vector 
    tol = 1e-6;
    max_iter = 10000;

    inject_error = 0;
    [~,~,iter,~] = pcg4(A, b, tol, max_iter, L, L', inject_error, 0, 0);
    noerror_converge = iter;   % number of iterations in error-free run
    disp(matrixname);
    toload = ['data\ML\PredictedSlowdown\protected_',matrixname,'_V2.dat'];
    protected = load(toload);

    
    %% load matrix file
    matrixfile = ['./matrices/', matrixname, '.mat'];
    load(matrixfile);
    A = Problem.A;
    
    %% load experimental data
    result_filename = ['./data/Step3_', matrixname, '_iter=', num2str(bitflip_iter), '.dat'];
    result = dlmread(result_filename);
    noerror_converges = result(:, 7);
    converges = result(:, 8);
    converge_ratios = converges./noerror_converges;
    num_exps = length(converge_ratios);
    
    %% load analysis data by A_row 2-norm
    
    protect_method = 'Arow2norm'; 
    analysis_filename = ['./data/ML/Analysis/', comments, '_', matrixname, '_', protect_method, '.dat'];
    slowdowns_Arow2norm = readmatrix(analysis_filename);
    slowdowns_Arow2norm = slowdowns_Arow2norm';
    mean_slowdowns_Arow2norm = mean(slowdowns_Arow2norm);
    std_slowdowns_Arow2norm = std(slowdowns_Arow2norm);
    sem_slowdowns_Arow2norm = std_slowdowns_Arow2norm/sqrt(num_exps);
    overheads1_Arow2norm = slowdowns_Arow2norm;
    for p = 1:num_protects
        protect = protects(p);
        overheads1_Arow2norm(:, p) = 100*((protect+1)*slowdowns_Arow2norm(:, p)-1);
    end
    mean_overheads1_Arow2norm = mean(overheads1_Arow2norm);
    disp(min(mean_overheads1_Arow2norm))
    std_overheads1_Arow2norm = std(overheads1_Arow2norm);
    sem_overheads1_Arow2norm = std_overheads1_Arow2norm/sqrt(num_exps);
    
    %% load analysis data by random
    
    protect_method = 'random'; 
    analysis_filename = ['./data/ML/Analysis/', comments, '_', matrixname, '_', protect_method, '.dat'];
    slowdowns_random = readmatrix(analysis_filename);
    slowdowns_random = slowdowns_random';
    mean_slowdowns_random = mean(slowdowns_random);
    std_slowdowns_random = std(slowdowns_random);
    sem_slowdowns_random = std_slowdowns_random/sqrt(num_exps);
    overheads1_random = slowdowns_random;
    for p = 1:num_protects
        protect = protects(p);
        overheads1_random(:, p) = 100*((protect+1)*slowdowns_random(:, p)-1);
    end
    mean_overheads1_random = mean(overheads1_random);
    std_overheads1_random = std(overheads1_random);
    sem_overheads1_random = std_overheads1_random/sqrt(num_exps);

    %% load analysis data for ML    
    protect_method = 'ML_Model'; 
    analysis_filename = ['./data/ML/Analysis/', comments, '_', matrixname, '_', protect_method, '_V2.dat'];
    protectname = ['data/ML/PredictedSlowdown/protected_',matrixname,'_V2.dat'];
    protected = readmatrix(protectname);
    slowdowns_ML = readmatrix(analysis_filename);
    mean_slowdowns_ML = mean(slowdowns_ML);
    overhead_ML = mean(100*(((N*noerror_converge\size(protected,1))+1)*slowdowns_ML-1));
    disp(overhead_ML)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %% plot figure
    figure;
    set(0,'defaultAxesTickLabelInterpreter','none');  
    hold on;
    plot(protects, mean_overheads1_Arow2norm, 'b', 'LineWidth',2);
    plot(protects+0.005, mean_overheads1_random, 'r', 'LineWidth',2);
    plot(protects,overhead_ML*ones(size(protects)),'LineWidth',2);
    %errorbar(protects+0.005, mean_overheads1_Arow2norm, CI95_overheads1_Arow2norm, 'CapSize',1, 'LineWidth',0.7, 'Color','blue');
    %errorbar(protects, mean_overheads1_random, CI95_overheads1_random, 'CapSize',1, 'LineWidth',0.7, 'Color','red');
    
    xlim([0, 1]);
    ylim([0, inf]);
    xlabel('Fraction of protection');
    ylabel('Average overhead (%)');
    title(matrixname, 'interpreter', 'none');
    hold off;
    set(gca,'FontSize',15);
    print(['./figures/overhead_compare_', matrixname, '_bitflip_iter=', num2str(bitflip_iter)], '-dpng');
    
    figure;
    hold on;
    plot(protects, mean_slowdowns_Arow2norm, 'b', 'LineWidth',2);
    plot(protects+0.005, mean_slowdowns_random, 'r', 'LineWidth',2);
    plot(protects,mean_slowdowns_ML*ones(size(protects)),'LineWidth',2)
    %errorbar(protects+0.005, mean_slowdowns_random, CI95_slowdowns_random, 'CapSize',1, 'LineWidth',0.7, 'Color','red');
    %errorbar(protects, mean_slowdowns_Arow2norm, CI95_slowdowns_Arow2norm, 'CapSize',1, 'LineWidth',0.7, 'Color', 'blue');
    legend('Row-2 Norm', 'Random','ML', 'Location', 'NorthEast');
    xlabel('Fraction of protection');
    ylabel('Average slowdown (x times)');
    title(matrixname, 'interpreter', 'none');
    xlim([0, 1]);
    hold off;
    set(gca,'FontSize',15);
    print(['./figures/slowdown_compare_', matrixname, '_bitflip_iter=', num2str(bitflip_iter)], '-dpng');
end 
    
end