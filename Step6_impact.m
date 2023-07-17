function Step6_impact

% This function plots Figure 1 of the paper. 

close all;
comments = 'impact';
color = 'b';

%matrices = {'cvxbqp1', 'thermal1', 'nd6k', ...
%    'bcsstk18', 'bodyy5', 'cbuckle', 'Pres_Poisson', 'bcsstk36', 'ct20stif', 'gyro_m', 't2dah_e', 'm_t1', 'msc23052', '2cubes_sphere', 'pwtk', 'G2_circuit', 'raefsky4', ...
%    'Trefethen_20000', 'vanbody','wathen100'};
matrices = {'bcsstk18'};
num_matrices = length(matrices);

bitflip_iter = 1;

for m = 1:num_matrices
    matrixname = matrices{m};
    disp(matrixname);
    
    %% load experimental data
    result_filename = ['./data/Step3_', matrixname, '_iter=', num2str(bitflip_iter), '.dat'];
    result = dlmread(result_filename);
    noerror_converges = result(:, 7);
    converges = result(:, 8);
    converge_ratios = converges./noerror_converges;
    
    %% plot figure
    figure;
    bar(sort(converge_ratios), 'FaceColor', color);
    title(matrixname, 'interpreter', 'none');
    %set(gca,'xscale','log');
    xlabel('Sample runs');
    ylabel('Slowdown (x times)');
    set(gca,'FontSize',15);
    
    figure_filename = ['./figures/', comments, '_', matrixname, '_bitflip_iter=', num2str(bitflip_iter)];
    print(figure_filename, '-dpng');
    %close all;
end 
    
end