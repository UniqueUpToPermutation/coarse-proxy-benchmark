debug_all_fine_solutions = [];
debug_sample_perm = [];
b_fine_solutions_loaded = false;

if isfile('BIE_fine_solutions.mat')
    [oracle, ss] = BenchmarkBIE('final');
    oracle.setCoarseResolution(128);
    oracle.setFineResolution(2048);

    strct = load('BIE_fine_solutions.mat');
    debug_all_fine_solutions = strct.debug_all_fine_solutions;
    debug_sample_perm = 1:size(ss, 2);
    b_fine_solutions_loaded = true;
end

eps = [0.0001, 0.00005, 0.00002, 0.00001, 0.000005, 0.000002, 0.000001, 0.0000005];
rbs = [];
for ep = eps
    % Each RBObject must have its own oracle and ss
    [oracle, ss] = BenchmarkBIE('final');
    oracle.setCoarseResolution(128);
    oracle.setFineResolution(2048);

    rb = RBObject(oracle, ss);
    rb.m_selection_epsilon = ep; % Set selection epsilon
    
    rb.b_enable_gramm_schmidt_for_operator_selection = true;
    rb.b_enable_run_parallel = true;
    
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
if ~b_fine_solutions_loaded
    save BIE_fine_solutions.mat debug_all_fine_solutions;
end

