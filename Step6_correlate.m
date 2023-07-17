function Step6_correlate

% This function plots Figure 2 of the paper

close all;
comments = 'correlate_slowdown';
mrk = 'o';
mrk_size = 15;
color = 'b';

%matrices = {'cvxbqp1', 'thermal1', 'nd6k', ...
%    'bcsstk18', 'bodyy5', 'cbuckle', 'Pres_Poisson', 'bcsstk36', 'ct20stif', 'gyro_m', 't2dah_e', 'm_t1', 'msc23052', '2cubes_sphere', 'pwtk', 'G2_circuit', 'raefsky4', ...
%    'Trefethen_20000', 'vanbody','wathen100'};
matrices = {'bcsstk18'};
num_matrices = length(matrices);

bitflip_iter = 400;

for m = 1:num_matrices
    matrixname = matrices{m};
    disp(matrixname);
    
    %% load experimental data
    result_filename = ['./data/Step3_', matrixname, '_iter=', num2str(bitflip_iter), '.dat'];
    result = dlmread(result_filename);
    A_row_2norms = result(:, 6);
    noerror_converges = result(:, 7);
    converges = result(:, 8);
    converge_ratios = converges./noerror_converges;
   
    %% slowdown figure
    figure;
    scatter(A_row_2norms, converge_ratios, mrk_size, mrk, 'filled', color);
    set(gca,'xscale','log');
    xlabel('Row 2-norm of matrix A');
    ylabel('Slowdown (x times)');
    titlename = [matrixname, '_iter : ', num2str(bitflip_iter)];
    title(titlename, 'interpreter', 'none');
    set(gca,'FontSize',15);
    hold off;
    figure_filename = ['./figures/', comments, '_', matrixname, '_bitflip_iter=', num2str(bitflip_iter)];
    print(figure_filename, '-dpng');
end 
    
end