[oracle, ss] = MoreSamplesSingularBIE();

oracle.setCoarseResolution(128);
oracle.setFineResolution(2048);

debug_all_fine_solutions = [];
debug_sample_perm = [];

eps = [0.0001, 0.00005, 0.00002, 0.00001, 0.000005, 0.000002, 0.000001, 0.0000005];
rbs = [];
for ep = eps
    rb = RBObject(oracle, ss);
    rb.m_selection_epsilon = ep; % Set selection epsilon
    rb.b_enable_gramm_schmidt_for_operator_selection = true;
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