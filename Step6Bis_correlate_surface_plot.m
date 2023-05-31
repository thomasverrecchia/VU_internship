function Step6Bis_correlate_surface_plot

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

elements = [185,515,715,1577,2806,4220,6873,9812,11256,11424];
num_elements = length(elements);
display_result = [];

for m = 1:num_matrices
    figure;
    for k = 1:num_elements
        matrixname = matrices{m};
    
        %% load experimental data
        result_filename = ['./data/Step3Bis_', matrixname, '_element=', num2str(elements(k)), '.dat'];
        result = dlmread(result_filename);
        bitflip_iter = result(:, 3);
        norm_2_A = result(:, 6);
        noerror_converges = result(:, 7);
        converges = result(:, 8);
        converge_ratios = converges./noerror_converges;

        %% slowdown figure
        
        scatter(bitflip_iter, norm_2_A,1:30,converge_ratios,'filled');
        set(gca,'xscale','log');
        xlabel('Injection iteration');
        ylabel('Row 2-norm');
        h = colorbar;
        ylabel(h, 'Slowdown (x times)')
        set(gca,'FontSize',15);
        hold on;            
    end
    hold off;
    figure_filename = ['./figures/', comments, '_', matrixname, '_surface_plot'];
    print(figure_filename, '-dpng');
end 
    
end