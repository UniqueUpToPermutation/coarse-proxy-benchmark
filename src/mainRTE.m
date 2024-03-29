debug_all_fine_solutions = [];
debug_sample_perm = [];
b_fine_solutions_loaded = false;

if isfile('RTE_fine_solutions.mat')
    [oracle, ss] = BenchmarkRTE('final');
    strct = load('RTE_fine_solutions.mat');
    debug_all_fine_solutions = strct.debug_all_fine_solutions;
    debug_sample_perm = 1:size(ss, 2);
    b_fine_solutions_loaded = true;
end

eps = [0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002, 0.0001];
rbs = [];
for ep = eps
    % Each RBObject must have its own oracle and ss
    [oracle, ss] = BenchmarkRTE('final');

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
save result_table_RTE.mat result_table;
if ~b_fine_solutions_loaded
    save RTE_fine_solutions.mat debug_all_fine_solutions;
end

