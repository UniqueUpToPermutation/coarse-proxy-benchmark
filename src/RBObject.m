classdef RBObject < handle

    properties (Constant = true)      
        % Defaults
        g_selection_epsilon_default = 1E-2;
        g_operator_selection_multiplier_default = 1.5;
        g_n_operator_samples_default = 10;
        g_b_additional_operator_skeletons = true;
        g_b_additional_skeleton_solutions = true;
        g_b_enable_normalized_qr = false;
        g_b_enable_gramm_schmidt = false;
        g_b_enable_normalized_qr_operator_samples = false;
        g_b_ignore_operator_sampling_time = false;
        g_b_enable_gramm_schmidt_for_operator_selection = false;
        g_b_enable_run_parallel = false;
        g_b_enable_compute_quantiles = false;
    end
   
    % Reduced basis properties
    properties
        % The reduced basis (n_f^2 x n_rb matrix)
        m_reduced_basis
        % The coarse solutions deemed "important" (n_c^2 x n_s matrix)
        m_coarse_solutions
        % The indices of the solutions deemed "important" (list)
        m_skeleton_indices
        % The fine solution analogues of the important coarse solutions
        % (n_f^2 x n_s matrix)
        m_fine_solutions
        % The fine left hand operators 
        m_fine_l_operators
        % Auxiliary data generated during fine solves
        m_fine_aux
        % The mixing matrix (n_s x n_sample_space_size matrix)
        m_mixing_matrix
        % The operators A projected onto the RB space (n_rb x n_rb x n_s
        % matrix)
        m_proj_l_operators
        % Important coarse operator problems
        m_operator_skeleton_indices
        % Samples we take from the fine operators, used for mixing matrix
        % construction
        m_fine_operator_samples
        % Additional basis samples computed with operators
        m_additional_skeleton_indices
    end
    
     % Problem data
    properties
        % The sample space for this problem
        m_sample_space
        % The oracle for this problem
        m_problem_oracle
        % The R-factor cutoff at which to select skeletons
        m_selection_epsilon = RBObject.g_selection_epsilon_default
        % This, multiplied by the selection epsilon gives the R-factor 
        % cutoff at which to select skeletons from operators
        m_operator_selection_multiplier = ...
            RBObject.g_operator_selection_multiplier_default
        % Number of rows to sample for mixing matrix
        m_n_operator_samples = RBObject.g_n_operator_samples_default
    end
    
    % Flags
    properties
        % Whether or not to compute additional skeletons using operator
        % samples
        b_enable_additional_operator_skeletons = ...
            RBObject.g_b_additional_operator_skeletons
        % Whether or not to add additional skeleton solutions on top of
        % additional operator skeletons.
        b_enable_additional_skeleton_solutions = ...
            RBObject.g_b_additional_skeleton_solutions
        % Whether or not to normalize input data to QR factorizations on
        % coarse data
        b_enable_normalized_qr = ...
            RBObject.g_b_enable_normalized_qr;
        % Whether or not to normalize input data to QR factorizations on
        % operator sample data
        b_enable_normalized_qr_operator_samples = ...
            RBObject.g_b_enable_normalized_qr_operator_samples;
        % Whether to use modified gramm-schmidt to compute the operator
        % samples with skeleton operator samples projected out (for
        % additional operator selection).
        b_enable_gramm_schmidt = ...
            RBObject.g_b_enable_gramm_schmidt;
        
        % Whether to ignore operator sampling time
        b_ignore_operator_sampling_time = ...
            RBObject.g_b_ignore_operator_sampling_time;
        
        % Whether to use modified gramm schmidt to actually select
        % additional skeleton operators
        b_enable_gramm_schmidt_for_operator_selection = ...
            RBObject.g_b_enable_gramm_schmidt_for_operator_selection
        
        % Whether or not to run operations in parallel
        b_enable_run_parallel = ...
            RBObject.g_b_enable_run_parallel;
        
        % Whether or not to compute error quantiles. (requires statistics or
        % ML toolkit)
        b_enable_compute_quantiles = ...
            RBObject.g_b_enable_compute_quantiles;
    end
 
    % DEBUG properties
    properties
       % An array of all fine solutions (n_f^2 x n_sample_space matrix)
       debug_all_fine_solutions
       % An array of all true operators A projected onto the RB space
       % (n_rb x n_rb x n_sample_space matrix)
       debug_proj_l_operators
       % Permutation of subsampled debug data
       debug_sample_perm
       % DEBUG: Drop the fine solutions selected through the coarse model
       b_debug_drop_solutions = false
       % DEBUG: Timing data
       debug_timing_data;
       % DEBUG: Debug diagnostics
       debug_diagnostic_data;
       % DEBUG: Error diagnostic data
       debug_diagnostic_error_data;
       % DEBUG: Coarse diagonal R values
       debug_coarse_r_diags;
       % DEBUG: Reduced basis singular values
       debug_rb_sing_values;
    end
    
    % Cleaning methods
    methods
        function dirtyRB(this)
            this.m_reduced_basis = [];
            this.m_mixing_matrix = [];
        end
        
        function clean(this)
           this.dirtyRB();
           this.m_fine_solutions = [];
           this.m_additional_skeleton_indices = [];
           this.m_coarse_solutions = [];
           this.m_fine_aux = {};
           this.m_fine_operator_samples = [];
           this.m_skeleton_indices = [];
           this.m_fine_l_operators = {};
           this.m_proj_l_operators = [];
           this.debug_timing_data = [];
           this.m_fine_operator_samples = [];
           this.clean_unnecessary();
        end
        
        function clean_debug(this)
           this.debug_all_fine_solutions = [];
           this.debug_proj_l_operators = [];
        end
        
        function clean_unnecessary(this)
           this.clean_debug();
           this.m_fine_solutions = [];
           this.m_coarse_solutions = [];
           this.m_fine_aux = {};
           this.m_fine_l_operators = {};
           this.m_fine_operator_samples = [];
        end
    end
    
    % Protected methods for internal routines
    methods (Access = protected)
        % Compute all coarse solutions using the problem oracle
        function [coarse_solutions] = compute_coarse(this, permutation)
            disp('Computing coarse solutions...');
            
            n = length(permutation);
            n_c2 = this.m_problem_oracle.getCoarseSolutionSize();
            coarse_solutions = zeros(n_c2, n);
            
            if ~this.b_enable_run_parallel
                for i=1:n
                    fprintf('Solving coarse system %i...\t(%i / %i)\r', ...
                        permutation(i), i, n);
                    omega = this.m_sample_space(:, permutation(i));
                    coarse_solutions(:, i) = ...
                        this.m_problem_oracle.solveCoarse(omega);
                end
            else
                sample_space = this.m_sample_space;
                problem_oracle = this.m_problem_oracle;
                parfor i=1:n
                    fprintf('Solving coarse system %i...\t(%i / %i)\r', ...
                        permutation(i), i, n);
                    omega = sample_space(:, permutation(i));
                    coarse_solutions(:, i) = ...
                        problem_oracle.solveCoarse(omega);
                end
            end
        end
        
        % Compute all fine solutions using the problem oracle
        function [fine_solutions, l_operators, auxs] = ...
                compute_fine(this, permutation)

            disp('Computing fine solutions...');

            n = length(permutation);
            n_f2 = this.m_problem_oracle.getFineSolutionSize();
            fine_solutions = zeros(n_f2, n);
            l_operators = cell(1, n);
            auxs = cell(1, n);

            if ~this.b_enable_run_parallel
                % Solve fine systems in serial
                for i=1:n
                    fprintf('Solving fine system %i...\t(%i / %i)\r', ...
                        permutation(i), i, n);
                    omega = this.m_sample_space(:, permutation(i));
                    [uC, op, aux] = this.m_problem_oracle.solveFine(omega);
                    fine_solutions(:, i) = uC;
                    l_operators{i} = op;
                    auxs{i} = aux;
                end
            else
                % Solve fine systems in parallel
                sample_space = this.m_sample_space;
                problem_oracle = this.m_problem_oracle;
                parfor i=1:n
                    fprintf('Solving fine system %i...\t(%i / %i)\r', ...
                        permutation(i), i, n);
                    omega = sample_space(:, permutation(i));
                    [uC, op, aux] = problem_oracle.solveFine(omega);
                    fine_solutions(:, i) = uC;
                    l_operators{i} = op;
                    auxs{i} = aux;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes additional skeletons by examining sampled rows and
        % the diagonal of our sample operators for each sample. 
        % First projects out the operator samples for the skeletons 
        % we've already selected.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [indices] = compute_additional_operator_skeletons(this)
            
            m_skeleton_op_samples = this.m_fine_operator_samples(:, ...
                this.m_skeleton_indices);
            
            eps = this.m_selection_epsilon * ...
                this.m_operator_selection_multiplier;
            
            fine_op_samples = this.m_fine_operator_samples;
            % For normalized qr, normalize the columns of the fine samples
            % to 1.
            if this.b_enable_normalized_qr_operator_samples
                fine_op_samples = fine_op_samples ./ vecnorm(fine_op_samples);
            end
            
            if ~this.b_enable_gramm_schmidt
                % Create a projector for the operators we've seen
                [Q, ~, ~] = qr(m_skeleton_op_samples, 0);

                % Project away stuff we've already captured in our skeletons
                m_op_samples_proj = fine_op_samples - ...
                    Q * (Q' * fine_op_samples);
            else
                % Use modified gramm-schmidt to project out skeletons
                m_op_samples_proj = fine_op_samples;
                for i_projector = 1:size(m_skeleton_op_samples, 2)
                    projector = m_skeleton_op_samples(:, i_projector);
                    projector = projector / norm(projector);
                    m_op_samples_proj = m_op_samples_proj - ...
                        projector * (projector' * m_op_samples_proj);
                end
            end
            
            if ~this.b_enable_normalized_qr_operator_samples
                base = max(vecnorm(this.m_fine_operator_samples));
                fprintf('Largest energy in column samples: %f\n', base);
            else
                base = 1.0;
            end
            
            % Select important columns
            if ~this.b_enable_gramm_schmidt_for_operator_selection
                [~, R, E] = qr(m_op_samples_proj, 0);
                aux = diag(abs(R));
                fprintf('Largest energy in projected column samples: %f\n', ...
                    aux(1));
                nSk = sum(aux > (base * eps));
                indices = sort(E(1:nSk));
            else
                % Use hacked implementation of Gramm-Schmidt for early
                % termination
                residual_vecs = m_op_samples_proj;
                
                n_vecs = size(residual_vecs, 2);
                n_gramm_schmidt_vec_count = 0;
                indices = zeros(1, n_vecs);
                
                threshold = base * eps;
                while true
                    
                    %l2norms = sum((residual_vecs .* residual_vecs), 1);
                    l2norms = vecnorm(residual_vecs); l2norms = l2norms .* l2norms;
                    [largest_norm, largest_index] = max(l2norms);
                    
                    if largest_norm <= threshold
                        break;
                    end
                    
                    largest_vec = residual_vecs(:, largest_index);
                    largest_vec = largest_vec / norm(largest_vec);
                    residual_vecs = residual_vecs - largest_vec * ...
                        (largest_vec' * residual_vecs);
                    
                    n_gramm_schmidt_vec_count = ...
                        n_gramm_schmidt_vec_count + 1;
                    indices(n_gramm_schmidt_vec_count) = largest_index;
                end
                
                indices = indices(1:n_gramm_schmidt_vec_count);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Merge additional skeletons into the skeletons we've already
        % selected.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function merge_skeletons(this, additional)
            new = setdiff(additional, this.m_skeleton_indices);
           
            if ~isempty(new)
                if (this.b_enable_additional_skeleton_solutions)
                    fprintf('Merging additional %i/%i basis elements...\n', ...
                        length(new), size(this.m_sample_space, 2));

                    disp('Solving additional skeleton fine systems...');
                    [uCs, l_ops, aux] = this.compute_fine(new);

                    disp('Merging fine solutions...');
                    this.m_fine_solutions = [this.m_fine_solutions, uCs];
                    this.m_fine_l_operators = [this.m_fine_l_operators, l_ops];
                    this.m_fine_aux = [this.m_fine_aux, aux];
                    this.m_skeleton_indices = [this.m_skeleton_indices, new];

                    dirtyRB(this);
                else
                    fprintf('Merging additional %i/%i operators...\n', ...
                        length(new), size(this.m_sample_space, 2));

                    [~, l_ops, ~] = this.compute_fine(new);

                    disp('Merging fine solutions...');
                    this.m_fine_l_operators = [this.m_fine_l_operators, l_ops];
                    this.m_skeleton_indices = [this.m_skeleton_indices, new];

                    dirtyRB(this);
                end
            else
                disp('Nothing to merge!');
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes mixing matrix for the projected operators
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [mixing] = compute_mixing_matrix(this, op_samples)
            % Solve to get the mixing (or interpolation) weight
            Tmp = op_samples(:, this.m_skeleton_indices);
            mixing = Tmp \ op_samples;
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes operators projected into the reduced basis space
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function compute_projected_operators(this)

            nBs = this.getReducedBasisSize();
            nSk = this.getSkeletonSolutionCount();
            
            this.m_proj_l_operators = zeros(nBs, nBs, nSk);
            for g=1:nSk
                this.m_proj_l_operators(:, :, g) = ...
                    this.m_problem_oracle.projectFineOperator(...
                    this.m_reduced_basis, this.m_fine_l_operators{g});
            end
        end
    end
    
    methods
        % Constructor
        % Note: the RBObject owns the problem oracle after this, multiple problem
        % oracles cannot be shared between RBObjects!
        function obj = RBObject(problem_oracle, sample_space)
            obj.m_problem_oracle = problem_oracle;
            obj.m_sample_space = sample_space;
            
            obj.debug_timing_data = cell2table(cell(0, 1));
            obj.debug_timing_data.Properties.VariableNames = {'Time'};
            obj.debug_diagnostic_data = cell2table(cell(0, 1));
            obj.debug_diagnostic_data.Properties.VariableNames = {'Value'};
        end
        
        % Helper functions
        function [n] = getSampleSpaceSize(this)
           n = size(this.m_sample_space, 2);
        end
        
        function [n] = getReducedBasisSize(this)
           n = size(this.m_reduced_basis, 2); 
        end
        
        function [n] = getSkeletonSolutionCount(this)
           n = size(this.m_fine_solutions, 2); 
        end
        
        function [n] = getSkeletonOperatorCount(this)
           n = length(this.m_skeleton_indices); 
        end
        
        % For debug use, add a profiling time to the record of times
        function setTime(this, task_name, time)
            if ismember(task_name, this.debug_timing_data.Properties.RowNames)
                this.debug_timing_data{task_name, 1} = time;
            else
                n = size(this.debug_timing_data, 1);
                this.debug_timing_data = [this.debug_timing_data; {time}];
                this.debug_timing_data.Properties.RowNames{n + 1} = task_name;
            end
        end
        
        % Get the time associated with a particular task
        function [time] = getTime(this, task_name)
            indx = find(ismember(this.debug_timing_data.Properties.RowNames, task_name));
            time = this.debug_timing_data{indx, 1};
        end
        
        % For debug use, add a diagnostic to the record of diagnostics
        function setDiagnostic(this, diagnostic_name, diagnostic_value)
             if ismember(diagnostic_name, this.debug_diagnostic_data.Properties.RowNames)
                this.debug_diagnostic_data{diagnostic_name, 1} = diagnostic_value;
            else
                n = size(this.debug_diagnostic_data, 1);
                this.debug_diagnostic_data = [this.debug_diagnostic_data; {diagnostic_value}];
                this.debug_diagnostic_data.Properties.RowNames{n + 1} = diagnostic_name;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assembles a reduced basis operator via interpolation.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [interp] = interpMatrix(this, index, mixing_mat, samples)

            tmp = mixing_mat(:, index);
            nSk = size(samples, 3);
            nR = this.getReducedBasisSize();

           	interp = zeros(nR, nR);
            for it = 1:nSk
                interp = interp + samples(:, :, it) * tmp(it);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assembles a reduced basis operator via interpolation.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [interp] = interpMatrixWithOffset(this, index, mixing_mat, ...
                samples)

            tmp = mixing_mat(:, index);
            nSk = size(samples, 3);

           	interp = this.m_problem_oracle.getAffineOffsetOperator(...
                this.m_reduced_basis);
            
            for it = 1:nSk
                interp = interp + samples(:, :, it) * tmp(it);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assembles a reduced basis vector via interpolation.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [interp] = interpVector(this, index, mixing_mat, samples)

            tmp = mixing_mat(:, index);
            nSk = size(samples, 2);
            nR = this.getReducedBasisSize();

           	interp = zeros(nR, 1);
            for it = 1:nSk
                interp = interp + samples(:, it) * tmp(it);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute a reduced basis for the sample space given
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeReducedBasis(this)
            
            globalClock = tic();
            samples_size = this.getSampleSpaceSize();
            
            this.setDiagnostic('Sample Space Size', this.getSampleSpaceSize());
            this.setDiagnostic('Selection Epsilon', this.m_selection_epsilon);
            
            % Compute all coarse solutions
            if isempty(this.m_coarse_solutions)
                clock = tic();
                this.m_coarse_solutions = compute_coarse(this, 1:samples_size);
                this.setTime('Compute Coarse Solutions', toc(clock));
            end
                    
            % Extract important basis elements
            if isempty(this.m_skeleton_indices)
                clock = tic();
                [this.m_skeleton_indices, this.debug_coarse_r_diags] = ...
                    extract_basis(this.m_coarse_solutions, ...
                    this.m_selection_epsilon, ...
                    this.b_enable_normalized_qr);
                 this.setTime('Compute Skeleton Indices', toc(clock));
            end        
            
            % Solve corresponding fine systems
            if isempty(this.m_fine_solutions)
                clock = tic();
                [uCs, ops, auxs] = compute_fine(this, this.m_skeleton_indices);
                this.m_fine_solutions = uCs;
                this.m_fine_l_operators = ops;
                this.m_fine_aux = auxs;
                this.setTime('Compute Fine Solutions', toc(clock));
                this.setDiagnostic('Initial Fine Solve Count', ...
                    size(this.m_skeleton_indices, 2)); 
            end
            
            % Sample the fine operators
            if isempty(this.m_fine_operator_samples)
                clock = tic();
                disp('Sampling fine operators...');
                this.m_fine_operator_samples = ...
                    this.m_problem_oracle.sampleFineOperators(...
                        this.m_sample_space, this.m_n_operator_samples);
                this.setTime('Sample Fine Operators', toc(clock));
                this.setDiagnostic('Fine Operator Sample Count', ...
                    this.m_n_operator_samples); 
            end
            
            % Add additional skeletons using samples from operators
            if this.b_enable_additional_operator_skeletons

                this.setDiagnostic('Operator Selection Multiplier', ...
                    this.m_operator_selection_multiplier);
                
                % DEBUG: Drop fine solutions computed from before
                if this.b_debug_drop_solutions
                    disp('DEBUG: Dropping solutions...');
                    this.clean_debug_drop_solutions();
                end
                
                clock = tic();
                disp('Computing operator-informed skeletons...');
                if isempty(this.m_additional_skeleton_indices)
                    this.m_additional_skeleton_indices = ...
                        this.compute_additional_operator_skeletons();
                end
                this.setTime('Compute Additional Skeleton Indices', toc(clock));
                this.setDiagnostic('Additional Skeleton Count', ...
                    size(this.m_additional_skeleton_indices, 2));
                
                clock = tic();
                this.merge_skeletons(this.m_additional_skeleton_indices);
                this.setTime('Compute Additional Skeleton Solutions', toc(clock));
                this.setDiagnostic('Total Skeleton Count', ...
                    size(this.m_skeleton_indices, 2));
            end
            
            % Compute mixing matrix
            if isempty(this.m_mixing_matrix)
                clock = tic();
                disp('Computing mixing matrix...');
                this.m_mixing_matrix = this.compute_mixing_matrix(...
                    this.m_fine_operator_samples);
                this.setTime('Compute Mixing Matrix', toc(clock));
            end
   
            % Orthogonalize important fine systems
            if isempty(this.m_reduced_basis)
                clock = tic();
                disp('Computing reduced basis...');
                fprintf('Skeleton size (solutions): %i/%i\n', ...
                    this.getSkeletonSolutionCount(), this.getSampleSpaceSize());
                fprintf('Skeleton size (operators): %i/%i\n', ...
                    this.getSkeletonOperatorCount(), this.getSampleSpaceSize());
                
                % Compute reduced basis from skeleton samples
                [this.m_reduced_basis, this.debug_rb_sing_values] = ...
                    compute_ortho_basis(this.m_fine_solutions, ...
                    this.m_selection_epsilon);
                
                fprintf('Reduced basis dimension: %i/%i\n', ...
                    this.getReducedBasisSize(), this.getSampleSpaceSize());
                this.setTime('Reduced Basis POD', toc(clock));
                this.setDiagnostic('Reduced Basis Dimension', ...
                    this.getReducedBasisSize());
            end

            % Project fine operators onto RB space
            if isempty(this.m_proj_l_operators)
                clock = tic();
                disp('Computing projected operators...');
                this.compute_projected_operators();
                this.setTime('Project Operators', toc(clock));
            end
            
            % Compute auxiliary data
            clock = tic();
            this.m_problem_oracle.computeAuxData(this);
            this.setTime('Compute Auxiliary Data', toc(clock));
            
            % Compute total time
            totalTime = toc(globalClock);
            if this.b_ignore_operator_sampling_time
                totalTime = totalTime - this.getTime('Sample Fine Operators');
            end
            this.setTime('Compute Reduced Basis (Total)', totalTime);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes a reduced basis approximation to the solution for
        % a given sample in the sample space.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [solution] = computeRBSolution(this, sample_index)
            
            if isempty(this.m_reduced_basis)
                throw(MException('RBObject:nobasis', ...
                    'Reduced basis not available! Run computeReducedBasis'));
            end

            omega = this.m_sample_space(:, sample_index);
            
            % Form the interpolated RB operator
            proj_l_operator = this.interpMatrixWithOffset(sample_index, ...
                this.m_mixing_matrix, this.m_proj_l_operators);
            
            % Use the oracle to form the right hand side
            proj_right_hand = ...
                this.m_problem_oracle.assembleRightHand(this, omega, ...
                sample_index);

            % Solve reduced problem
            rb_solution = proj_l_operator \ proj_right_hand;

            % Return solution
            solution = this.m_reduced_basis * rb_solution;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % View a random RB approximation result
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function viewRandomResult(this, figure_id)
            if ~exist('figure_id', 'var')
                figure_id = 1;
            end
            
            % Solve a random instance and output the solutions
            sample_index = randi(this.getSampleSpaceSize());
            viewResult(this, sample_index, figure_id);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % View the RB approximation result for the sample_index equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function viewResult(this, sample_index, figure_id)
            if ~exist('figure_id', 'var')
                figure_id = 1;
            end
            
            if isempty(this.m_reduced_basis)
                throw(MException('RBObject:nobasis', ...
                    'Reduced basis not available! Run computeReducedBasis'));
            end
            
            omega = this.m_sample_space(:, sample_index);
           
            [fine_solution, ~, ~] = this.m_problem_oracle.solveFine(omega);
            rb_solution = this.computeRBSolution(sample_index);
            
            % View solutions via problem oracle
            this.m_problem_oracle.viewResult(rb_solution, fine_solution, ...
                figure_id);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % View the reduced basis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function viewReducedBasis(this, figure_id)
            if ~exist('figure_id', 'var')
                figure_id = 1;
            end
            
            if ~isempty(this.m_reduced_basis)
                % Show the actual basis functions
                this.m_problem_oracle.viewBasis(this.m_reduced_basis, ...
                    figure_id);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute all fine solutions and true projected operators for
        % debug diagnostics.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeDebugData(this, subsample_factor)
            if ~exist('subsample_factor', 'var')
               subsample_factor = 1; 
            end
            
            disp('Computing DEBUG data...');
            
            nSk = this.getSampleSpaceSize();
            nR = this.getReducedBasisSize();
            nF2 = this.m_problem_oracle.getFineSolutionSize();

            subsample_factor = max(min(subsample_factor, 1), 0);
            nSk_sub = int32(nSk * subsample_factor);
            range = sort(randperm(nSk, nSk_sub));

            % Initialize all debug data to zero
            tmp_debug_all_fine_solutions = zeros(nF2, nSk_sub);
            tmp_debug_proj_l_operators = zeros(nR, nR, nSk_sub);
            proj_l_operator_offset = ...
                this.m_problem_oracle.getAffineOffsetOperator(...
                    this.m_reduced_basis);
            
            % Fill debug data by solving selected problems
            clock = tic();
            if ~this.b_enable_run_parallel
                for i=1:nSk_sub
                    fprintf('Solving fine system %i...\t(%i / %i)\r', range(i), ...
                        i, nSk_sub);

                    omega = this.m_sample_space(:, range(i));

                    [uC, MC, ~] = this.m_problem_oracle.solveFine(omega);

                    tmp_debug_all_fine_solutions(:, i) = uC;
                    tmp_debug_proj_l_operators(:, :, i) = ...
                        proj_l_operator_offset + ...
                        this.m_problem_oracle.projectFineOperator(...
                            this.m_reduced_basis, MC);
                end
            else
                sample_space = this.m_sample_space;
                problem_oracle = this.m_problem_oracle;
                reduced_basis = this.m_reduced_basis;
                parfor i=1:nSk_sub
                    fprintf('Solving fine system %i...\t(%i / %i)\r', range(i), ...
                        i, nSk_sub);

                    omega = sample_space(:, range(i));

                    [uC, MC, ~] = problem_oracle.solveFine(omega);

                    tmp_debug_all_fine_solutions(:, i) = uC;
                    tmp_debug_proj_l_operators(:, :, i) = ...
                        proj_l_operator_offset + ...
                        problem_oracle.projectFineOperator(...
                            reduced_basis, MC);
                end
            end

            this.debug_sample_perm = range;
            
            % Record timing data
            this.setTime('DEBUG: Total Fine Solution Time', toc(clock));
            this.setDiagnostic('DEBUG: Total Samples', nSk_sub);
            
            % Store data after computation is done
            this.debug_all_fine_solutions = tmp_debug_all_fine_solutions;
            this.debug_proj_l_operators = tmp_debug_proj_l_operators;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Output diagnostics on the error of projected true solutions.
        % (i.e. how well does the RB space approximate the solution space)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [stats] = outputRBSpaceDiagnostics(this, b_verbose)
            
            if ~exist('b_verbose', 'var')
               b_verbose = false; 
            end
            
            if isempty(this.debug_all_fine_solutions)
                throw(MException('RBObject:nodata', ...
                    'fine solutions not available!'));
            end

            uCs = this.debug_all_fine_solutions;
            rb = this.m_reduced_basis;
            proj = uCs - rb * (rb' * uCs);
            n = size(this.debug_all_fine_solutions, 2);
            errs = zeros(n, 1);
            for i = 1:n
                mag = norm(uCs(:,i));
                mag_proj = norm(proj(:,i));
                errs(i) = mag_proj / mag;
            end
            
            print_statistics(errs, b_verbose, 'projected fine solutions', ...
                this.b_enable_compute_quantiles);
            
            stats = make_statistics_struct(errs, this.b_enable_compute_quantiles);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Output diagnostics on interpolant operator error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [stats, stats_inv] = outputOperatorDiagnostics(this, b_verbose)
            
            if ~exist('b_verbose', 'var')
               b_verbose = false; 
            end
            
            if isempty(this.debug_all_fine_solutions)
                throw(MException('RBObject:nodata', ...
                    'fine solutions not available!'));
            end

            n = size(this.debug_proj_l_operators, 3);
            l_errs = zeros(n, 1);
            l_errs_inv = zeros(n, 1);
          
            % Ground truth operators
            true_l_ops = this.debug_proj_l_operators;
           
            % Assemble the interp operators and compare to ground truth
            for i = 1:n          
                op = this.interpMatrixWithOffset(this.debug_sample_perm(i), ...
                    this.m_mixing_matrix, this.m_proj_l_operators);
                true_op = true_l_ops(:, :, i);
                
                l_errs(i) = norm(op - true_op) / norm(op);
                
                op_inv = inv(op);
                true_op_inv = inv(true_op);
                
                l_errs_inv(i) = norm(op_inv - true_op_inv) / norm(op_inv);
            end

            % Print error statistics
            print_statistics(l_errs, b_verbose, 'projected LH operators', ...
                this.b_enable_compute_quantiles);
            print_statistics(l_errs_inv, b_verbose, ...
                'inverses of projected LH operators', ...
                this.b_enable_compute_quantiles);
            
            % Output error statistics
            stats = make_statistics_struct(l_errs, ...
                this.b_enable_compute_quantiles);
            stats_inv = make_statistics_struct(l_errs_inv, ...
                this.b_enable_compute_quantiles);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Output diagnostics on final RB approximation error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [stats] = outputRBErrorDiagnostics(this, b_verbose)
            
            clock = tic();
            
            if ~exist('b_verbose', 'var')
               b_verbose = false; 
            end
            
            n = size(this.debug_all_fine_solutions, 2);
            errs = zeros(n, 1);
            
            % Solve all problems via rb and compare to true solutions
            for i = 1:n
                rb_solution = this.computeRBSolution(this.debug_sample_perm(i));
                fine_solution = this.debug_all_fine_solutions(:, i);

                errs(i) = norm(rb_solution - fine_solution) / ...
                    norm(fine_solution);
            end
            
            print_statistics(errs, b_verbose, 'rb solutions', ...
                this.b_enable_compute_quantiles);
            
            % Output error statistics
            stats = make_statistics_struct(errs, ...
                this.b_enable_compute_quantiles);
            
            % Record times
            this.setTime('DEBUG: Total Reduced Basis Time', toc(clock));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Output diagnostics on the error of the computed RB space
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [combined_error_stats] = outputDiagnostics(this, b_verbose)
            if ~exist('b_verbose', 'var')
               b_verbose = false; 
            end
            
            stats_RBspace = outputRBSpaceDiagnostics(this, b_verbose);
            [stats_operator, stats_operator_inv] = ...
                outputOperatorDiagnostics(this, b_verbose);
            stats_RBerror = outputRBErrorDiagnostics(this, b_verbose);
            
            combined_error_stats = [stats_RBspace; ...
                stats_operator; ...
                stats_operator_inv; ...
                stats_RBerror];
            combined_error_stats.Properties.RowNames = {'RB Projection Error', ...
                'Operator Projection Error', ...
                'Projected Operator Inverse Error', ...
                'Final Error' };
            
            this.debug_diagnostic_error_data = combined_error_stats;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Output diagnostics but reuse fine solutions from old computation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [combined_error_stats] = outputDiagnosticsReuseDebugData(this, ...
                debug_all_fine_solutions, debug_sample_perm, b_verbose)
            if ~exist('b_verbose', 'var')
               b_verbose = false; 
            end
            
            this.debug_all_fine_solutions = debug_all_fine_solutions;
            this.debug_sample_perm = debug_sample_perm;
            
            stats_RBspace = outputRBSpaceDiagnostics(this, b_verbose);
            stats_RBerror = outputRBErrorDiagnostics(this, b_verbose);
            
            combined_error_stats = [stats_RBspace; ...
                stats_RBerror];
            combined_error_stats.Properties.RowNames = {'RB Projection Error', ...
                'Final Error' };
            
            this.debug_diagnostic_error_data = combined_error_stats;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to output error statistics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_statistics(errs, b_verbose, str_name, b_enable_compute_quantiles)
    % Print error statistics
    fprintf('\nDEBUG: Relative errors of ');
    fprintf(str_name);
    fprintf(':\n');
    if b_verbose
        for i = 1:n
            fprintf('Solution %i:\t%f (%.2f %%)\n', i, ...
                k_errs(i), 100 * k_errs(i));
        end
        fprintf('\n');
    end

    mx = max(errs);
    mn = mean(errs);
    mdn = median(errs);
    stdev = std(errs);
    
    fprintf('Maximum error:\t\t%f (%.2f %%)\n', mx, 100 * mx);
    fprintf('Average error:\t\t%f (%.2f %%)\n', mn, 100 * mn);
    fprintf('Median error:\t\t%f (%.2f %%)\n', mdn, 100 * mdn);
    fprintf('Standard Deviation:\t%f (%.2f %%)\n', stdev, 100 * stdev);
    
    if b_enable_compute_quantiles
        q1 = quantile(errs, 0.25);
        q3 = quantile(errs, 0.75);
        fprintf('Q1 error:\t\t%f (%.2f %%)\n', q1, 100 * q1); 
        fprintf('Q3 error:\t\t%f (%.2f %%)\n', q3, 100 * q3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to output statistical data struct about errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stat_struct] = make_statistics_struct(errs, b_enable_compute_quantiles)
    stat_max = max(errs);
    stat_mean = mean(errs);
    stat_median = median(errs);
    stat_stdev = std(errs);
    
    if b_enable_compute_quantiles
        stat_q1 = quantile(errs, 0.25);
        stat_q3 = quantile(errs, 0.75);
        stats = {stat_mean, stat_stdev, stat_q1, stat_median, stat_q3, stat_max};
        stat_names = {'mean', 'stdev', 'Q1', 'median', 'Q3', 'max'};
    else
        stats = {stat_mean, stat_stdev, stat_median, stat_max};
        stat_names = {'mean', 'stdev', 'median', 'max'};
    end
    
    stat_struct = cell2table(stats);
    stat_struct.Properties.VariableNames = stat_names;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract important samples from the coarse solutions
%
% Parameters:
% coarse_solutions: coarse solutions from which to extract
% epsilon: threshold for determining what is important
%
% Return values:
% basis_indices: indices of important basis elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [basis_indices, r_vals] = extract_basis(coarse_solutions, epsilon, b_normalize)
    fprintf('Selecting basis elements...\n');
    
    % Normalize columns first!
    if b_normalize
        coarse_solution_norms = vecnorm(coarse_solutions);
        coarse_solutions = coarse_solutions ./ coarse_solution_norms;
    end

    % Take pivoted QR decomposition of coarse solutions
    [~, R, E] = qr(coarse_solutions, 0);
    aux = diag(abs(R));
    
    % Select all columns with R factor greater than a factor of epsilon
    % times first R factor
    nSk = sum(aux > aux(1) * epsilon);
    
    fprintf('Selected %i/%i basis elements!\n', nSk, ...
        size(coarse_solutions, 2));
    
    % Return the corresponding indices
    basis_indices = sort(E(1:nSk));
    
    % Return the r values
    r_vals = aux;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes an orthogonal basis set for a given space spanned by the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ortho_basis, singular_vals] = compute_ortho_basis(basis, epsilon)
    [U, S, ~] = svd(basis,'econ');    
    s = diag(S);
    nBs = sum(s > s(1) * epsilon);
    ortho_basis = U(:, 1:nBs);
    
    % Return singular values
    singular_vals = s;
end