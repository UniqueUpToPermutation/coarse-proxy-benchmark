if exist('flag_load_flam', 'var')
    cd ../FLAM/
    startup
    cd ../RB-Framework/
end

[oracle, ss] = GuassianRTEHarderNonSym();

debug_all_fine_solutions = [];
debug_sample_perm = [];

eps = [0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002];
rbs = [];
for ep = eps
    rb = RBObject(oracle, ss);
    rb.m_selection_epsilon = ep; % Set selection epsilon
    rb.computeReducedBasis();
    if isempty(debug_all_fine_solutions)
        rb.computeDebugData();
        debug_all_fine_solutions = rb.debug_all_fine_solutions;
        debug_sample_perm = rb.debug_sample_perm;
    end
    rb.outputDiagnosticsReuseDebugData(debug_all_fine_solutions, ...
        debug_sample_perm);
    rb.clean_unnecessary(); % Save storage space
    rbs = [rbs, rb];
end

result_table = table(eps', rbs');
save result_table_BIE.mat result_table;

