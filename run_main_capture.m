clearvars;
close all;
clc;

out_dir = fullfile(pwd, 'generated_results', datestr(now, 'yyyymmdd_HHMMSS_main'));
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

diary(fullfile(out_dir, 'run_main_log.txt'));
main;

result_file = fullfile(out_dir, 'result_main.mat');
save(result_file, 'alg', 'alg_list', 'alg_name', 'data_obj', 'data_lambda', ...
    'dll', 'x', 'V0', 'x_min', 'nx', 'ny', 'nm', 'nd', 'eeee', 'rho', ...
    'dmm0', 'epsilon', 'scale0', 'coord_x', 'irr', 't_iter');

figs = findall(0, 'Type', 'figure');
for k = 1:numel(figs)
    fig = figs(k);
    fig_file = fullfile(out_dir, sprintf('main_figure_%02d.png', fig.Number));
    try
        exportgraphics(fig, fig_file, 'Resolution', 200);
    catch
        saveas(fig, fig_file);
    end
end

for i_alg = alg
    last_idx = find(data_obj(:, i_alg) ~= 0, 1, 'last');
    if ~isempty(last_idx)
        fprintf('SUMMARY Alg %s: Eig1 %.12f Eig2 %.12f Eig3 %.12f volume %.12f\n', ...
            alg_list(i_alg), data_lambda(last_idx, i_alg, 1), ...
            data_lambda(last_idx, i_alg, 2), data_lambda(last_idx, i_alg, 3), dll' * x);
    end
end

fprintf('Saved main.m captured results to: %s\n', out_dir);
diary off;
