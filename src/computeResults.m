function [tb] = computeResults(result_table)

    avg_err = [];
    std_err = [];
    eps = [];
    time_build_rb = [];
    time_compute_fine = [];
    time_compute_rb = [];
    time_ratio = [];
    skeletons = [];
    rb_dim = [];
    time_rb_total = [];

    fine_time = result_table{1, 2}.debug_timing_data{12, 1};

    for i = 1:size(result_table, 1)
        
        rb = result_table{i, 2};
        if size(rb.debug_timing_data, 1) > 12
            indx_offset = 0;
        else
            indx_offset = 1;
        end
        
        eps = [eps, result_table{i, 1}];
        avg_err = [avg_err, rb.debug_diagnostic_error_data.mean(2)];
        std_err = [std_err, rb.debug_diagnostic_error_data.stdev(2)];
        time_build_rb = [time_build_rb, rb.debug_timing_data{11, 1}];
        time_compute_fine = [time_compute_fine, fine_time];
        time_compute_rb = [time_compute_rb, rb.debug_timing_data{13 - indx_offset, 1}];
        time_rb_total = [time_rb_total, time_build_rb(i) + time_compute_rb(i)];
        time_ratio = [time_ratio, time_compute_fine(i) / time_rb_total(i)];
        skeletons = [skeletons, rb.debug_diagnostic_data{7, 1}];
        rb_dim = [rb_dim, rb.debug_diagnostic_data{8, 1}];
    end

    skeletons_dbl = double(skeletons);
    rb_dim_dbl = double(rb_dim);

    tb = table(eps', skeletons_dbl', rb_dim_dbl', time_build_rb', time_compute_rb', ...,
        time_compute_fine', time_ratio', avg_err');
    tb.Properties.VariableNames = {'eps', 'SKELETON_COUNT', 'RB_DIM', 'TIME_BUILD_RB', 'TIME_COMPUTE_RB', ...
       'TIME_FINE', 'TIME_RATIO', 'AVG_ERR'};

    figure(1);

    rat = 1./eps;
    loglog(rat, avg_err, '-o');

    input = struct();
    input.data = tb;
    input.dataFormat = {'%.4f',1, '%i', 2, '%.2f', 1, '%.2f', 1, '%.2f', 1, '%.1f', 1, '%.4f', 1};

    ltx = latexTable(input);

end
