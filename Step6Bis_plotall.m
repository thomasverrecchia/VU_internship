function Step6_plotall

% This function plots Figure 9 in the paper

close all;

% matrices = {'cvxbqp1','thermal1','nd6k', ...
%     't2dah_e','bcsstk18','cbuckle','Pres_Poisson','gyro_m','bodyy5','raefsky4','msc23052','bcsstk36','ct20stif','m_t1','2cubes_sphere','G2_circuit','pwtk',...
%     'Trefethen_20000','vanbody','wathen100'};
matrices = {'bcsstk18'};
num_matrices = length(matrices);

bitflip_iter = 1;
protects = [0:0.01:1];
num_protects = length(protects);

overheads = ones(4, num_matrices);

for m = 1:num_matrices
    matrixname = matrices{m};
    disp(matrixname);
    
    %% load matrix file
    matrixfile = ['./matrices/', matrixname, '.mat'];
    load(matrixfile);
    A = Problem.A;
    
    %% load matrix norms
    norms_filename = ['./matrices/', matrixname, '_norms.mat'];
    load(norms_filename);
    all_A_row_2norms = A_norms(:, 1);
    all_sorted_A_row_2norms = sort(all_A_row_2norms);
    
    %% load experimental data
    result_filename = ['./data/Step3_', matrixname, '_iter=', num2str(bitflip_iter), '.dat'];
    result = dlmread(result_filename);
    N = length(A);
    A_row_2norms = result(:, 6);
    noerror_converges = result(:, 7);
    converges = result(:, 8);
    converge_ratios = converges./noerror_converges;
    num_exps = length(converge_ratios);
    
    %% load analysis data by A_row 2-norm
    protect_method = 'Arow2norm'; 
    analysis_filename = ['./data/Step5_', matrixname, '_iter=', num2str(bitflip_iter), '_', protect_method, '.dat'];
    slowdowns_Arow2norm = dlmread(analysis_filename);
    slowdowns_Arow2norm = slowdowns_Arow2norm';
    overheads1_Arow2norm = slowdowns_Arow2norm;
    for p = 1:num_protects
        protect = protects(p);
        overheads1_Arow2norm(:, p) = 100*((protect+1)*slowdowns_Arow2norm(:, p)-1);
    end
    mean_overheads1_Arow2norm = mean(overheads1_Arow2norm);
    
    %% curve fitting
    [p2,~,mu2] = polyfit(A_row_2norms, converge_ratios, 2);
    [p3,~,mu3] = polyfit(A_row_2norms, converge_ratios, 3);
    if p2(1) > 0
        p = p2;
        mu = mu2;
        n = 2;
    elseif p3(1) > 0 
        p = p3;
        mu = mu3;
        n = 3;
    else
        p = p2;
        mu = mu2;
        n = 2;
    end 
    %p = p2;
    %mu = mu2;
    %n = 2;
    
    min_A_row_2norms = min(A_row_2norms);
    max_A_row_2norms = max(A_row_2norms);
    x = [min_A_row_2norms:(max_A_row_2norms-min_A_row_2norms)/1000:max_A_row_2norms];
    y_fit = polyval(p,x,[],mu);
    
    predicted_converge_ratios = polyval(p,A_row_2norms,[],mu);
    temp = corrcoef(predicted_converge_ratios, converge_ratios);
    R2 = temp(1,2)^2;
    
    %% predict overhead
    all_sorted_converge_ratios = polyval(p,all_sorted_A_row_2norms,[],mu);
    expected_slowdowns = zeros(1, 101);
    expected_overheads = zeros(1, 101);
    index = 1;
    for i = 0:0.01:1 % fraction
        num_protected_elements = ceil(N*i);
        expected_slowdowns(index) = (1/N)*(sum(all_sorted_converge_ratios(1:N-num_protected_elements))+num_protected_elements);
        expected_overheads(index) = 100*(expected_slowdowns(index)*(1+i)-1);
        index = index + 1;
    end
    [min_expected_overhead, minind_expected_overhead] = min(expected_overheads);
    min_alg_overhead = mean_overheads1_Arow2norm(minind_expected_overhead);
    
    %% our algorithm's actual overhead
    %[minind_expected_overhead, min_alg_overhead]
    overheads(1, m) = min_alg_overhead; %min_mean_overheads1_Arow2norm;
    
    %% load analysis data by random
    %CI95 = tinv([0.975], num_exps-1);  % 95% confidence interval
    
    protect_method = 'random'; 
    analysis_filename = ['./data/Step5_', matrixname, '_iter=', num2str(bitflip_iter), '_', solver, '_', protect_method, '.dat'];
    slowdowns_random = dlmread(analysis_filename);
    slowdowns_random = slowdowns_random';
    overheads1_random = slowdowns_random;
    for p = 1:num_protects
        protect = protects(p);
        overheads1_random(:, p) = 100*((protect+1)*slowdowns_random(:, p)-1);
    end
    mean_overheads1_random = mean(overheads1_random);
    [min_mean_overheads1_random, minind_mean_overheads1_random] = min(mean_overheads1_random);
    
    overheads(2, m) = min_mean_overheads1_random;
    
    %% load zero-protect and all-protect data
    overheads(3, m) = mean_overheads1_Arow2norm(1); % no protect
    overheads(4, m) = mean_overheads1_Arow2norm(length(mean_overheads1_Arow2norm)); % all protect
end 

%% plot figure
figure_filename = ['./figures/plotall'];
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2, 1, 1);
x = [1:num_matrices];
b = bar(x, overheads); 
b(1).FaceColor = [0.4660 0.6740 0.1880];
b(2).FaceColor = [0 0.4470 0.7410];
b(3).FaceColor = [0.8500 0.3250 0.0980];
b(4).FaceColor = [0.4940 0.1840 0.5560];
ylabel('Average overheads (%)');
ylim([0, 200]);
xticks([1:num_matrices]);
xticklabels(matrices);
yticks([0:50:200]);
yticklabels({'0', '50', '100', '150', '200'});
xtickangle(45);
legend('Our selective-protection','Optimal random-protection', 'Zero-protection', 'Full-protection', 'Location', 'northoutside', 'Orientation', 'horizontal');
set(gca,'TickLabelInterpreter','none');

print(figure_filename, '-dpng');
%close all;
   
end