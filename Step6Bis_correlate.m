function Step6Bis_correlate

% This function plots Figure 2 of the paper

close all;
comments = 'correlate_slowdown_Bis';
mrk = 'o';
mrk_size = 15;
color = 'b';

%matrices = {'cvxbqp1', 'thermal1', 'nd6k', ...
%    'bcsstk18', 'bodyy5', 'cbuckle', 'Pres_Poisson', 'bcsstk36', 'ct20stif', 'gyro_m', 't2dah_e', 'm_t1', 'msc23052', '2cubes_sphere', 'pwtk', 'G2_circuit', 'raefsky4', ...
%    'Trefethen_20000', 'vanbody','wathen100'};
matrices = {'bcsstk18'};
num_matrices = length(matrices);

bitflip_element = 8235;

for m = 1:num_matrices
    matrixname = matrices{m};
    disp(matrixname);
    
    %% load experimental data
    result_filename = ['./data/Step3Bis_', matrixname, '_element=', num2str(bitflip_element), '.dat'];
    result = dlmread(result_filename);
    bitflip_iter = result(:, 3);
    noerror_converges = result(:, 7);
    converges = result(:, 8);
    converge_ratios = converges./noerror_converges;
   
    %% slowdown figure
    figure;
    scatter(bitflip_iter, converge_ratios, mrk_size, mrk, 'filled', color);
    set(gca,'xscale','log');
    xlabel('Injection iteration');
    ylabel('Slowdown (x times)');
    norm_2_A = result(:, 6);
    titlename = ['A_2_Norm =',num2str(norm_2_A(1))];
    title(titlename, 'interpreter', 'none');
    set(gca,'FontSize',15);
    hold off;
    figure_filename = ['./figures/', comments, '_', matrixname, '_bitflip_element=', num2str(bitflip_element)];
    print(figure_filename, '-dpng');
end 
    
end