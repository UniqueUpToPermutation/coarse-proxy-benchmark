classdef OracleRTE < ProblemOracle

    properties (Constant = true)
        % HIF Parameters
        g_occ = 64;
        g_p = 64;
        g_rank_or_tol = 1e-6;
        g_skip = 1;
        g_symm = 'n'; % not symmetric in this case
        g_verbose = 0;
        
        % Number of quadrature points
        g_quadrature_points = 5; 
        g_theta = (1:OracleRTE.g_p) * 2 * pi / OracleRTE.g_p;
        g_proxy = 1.5 * [cos(OracleRTE.g_theta); sin(OracleRTE.g_theta)];
        
        % Default coarse and fine grid dimensions and selection epsilon 
        g_n_coarse_grid_default = 32;
        g_n_fine_grid_default = 128;
    end
    
    % Problem data
    properties
        % The underlying force function
        m_force_func
        % Transmition coefficient function
        m_mu_t_func
        % Scattering coefficient function
        m_mu_s_func
        % The vectors mu_s K f projected onto the RB space (n_rb x n_s matrix)
        m_proj_mus_k_f_vectors
        % The vector f projected onto the RB space
        m_proj_f_vector
    end
    
    methods
        function obj = OracleRTE(force_func, mu_t_func, mu_s_func)
            
            obj@ProblemOracle(OracleRTE.g_n_coarse_grid_default, ...
                OracleRTE.g_n_fine_grid_default);
            
            obj.m_force_func = force_func;
            obj.m_mu_t_func = mu_t_func;
            obj.m_mu_s_func = mu_s_func;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Implementation of super class abstract methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [n] = getCoarseSolutionSize(this)
            n = this.m_n_coarse * this.m_n_coarse;
        end
        
        % Get the vector size of fine solutions 
        function [n] = getFineSolutionSize(this)
            n = this.m_n_fine * this.m_n_fine; 
        end
        
        function [solution] = solveCoarse(this, omega)
            solution = compute_coarse(this.m_force_func, this.m_mu_t_func, ...
                this.m_mu_s_func, omega, this.m_n_coarse);
        end
        
        % Solve fine systems
        % We use aux to store the diagonal elements of mu_s^{-1}
        % These are later used to compute the interpolation samples Q^T K(w) f
        % in the computeAux method
        function [solution, loperator, aux] = solveFine(this, omega)
            [loperator, aux, solution] = compute_fine(this.m_force_func, ...
                this.m_mu_t_func, this.m_mu_s_func, omega, this.m_n_fine);
        end
        
        function [proj_operator] = projectFineOperator(~, ...
                reduced_basis, operator)
            proj_operator = reduced_basis' * hifie_mv(operator, reduced_basis);
        end
        
        function [samples] = sampleFineOperators(this, sample_space, n_samples)
            N = this.m_n_fine * this.m_n_fine;
            pck = randperm(N, n_samples);
            nRand = size(sample_space, 2);

            % each column of Bpck contains several columns of the M matix 
           	samples = zeros(n_samples * N, nRand); 
            for g=1:nRand
                omega = sample_space(:, g);
                musfuntmp = @(x1, x2) this.m_mu_s_func(x1, x2, omega);
                mutfuntmp = @(x1, x2) this.m_mu_t_func(x1, x2, omega);
                
                [N, ~, param, ~] = setupparam(...
                    OracleRTE.g_quadrature_points, this.m_n_fine, ...
                    OracleRTE.g_proxy, musfuntmp, mutfuntmp);

                % Sample the operator mu_s K
                Mtmp = Bfunnonsym(1:N, pck, param, false);
                
                samples(:, g) = Mtmp(:);
            end
        end
        
        function [samples] = sampleEntireFineOperatorsDebug(this, sample_space)
            N = this.m_n_fine * this.m_n_fine;
            pck = 1:N;
            nRand = size(sample_space, 2);

            % each column of Apck contains several columns of the M matix 
            % and the diagonal of the D matrix
           	samples = zeros(N * N, nRand); 
            for g=1:nRand
                omega = sample_space(:, g);
                musfuntmp = @(x1, x2) this.m_mu_s_func(x1, x2, omega);
                mutfuntmp = @(x1, x2) this.m_mu_t_func(x1, x2, omega);

               	[N, ~, param, ~] = setupparam(...
                    OracleRTE.g_quadrature_points, this.m_n_fine, ...
                    OracleRTE.g_proxy, musfuntmp, mutfuntmp);

                % Sample the operator mu_s K
                Mtmp = Bfunnonsym(1:N, pck, param, false);
                
                samples(:, g) = Mtmp(:);
            end
        end
        
        % Compute interpolation samples for projected vectors Q^T K f
        % Note that K = mu_s^{-1} - A
        % And mu_s^{-1} is stored in rb_obj.m_fine_aux
        function computeAuxData(this, rb_object)
            
            nF = this.m_n_fine;
            [x1, x2] = ndgrid((1/2:nF) / nF);
            xs = [x1(:) x2(:)]';
            fC = this.m_force_func(xs(1,:), xs(2,:));        
            fC = fC(:);
            
            nBs = size(rb_object.m_reduced_basis, 2);
            nSk = size(rb_object.m_fine_l_operators, 2);
            
            this.m_proj_mus_k_f_vectors = zeros(nBs, nSk);
            
            for g=1:nSk
                % The auxiliary data contains the mu_s^{-1} diagonals
                 this.m_proj_mus_k_f_vectors(:, g) = ...
                    rb_object.m_reduced_basis' * (hifie_mv(rb_object.m_fine_l_operators{g}, fC));
            end
            
            this.m_proj_f_vector = rb_object.m_reduced_basis' * fC;
        end
        
         % Assemble the right hand of a reduced basis solve
         function [proj_f] = assembleRightHand(this, rb_object, ~, sample_index)
             proj_f = - rb_object.interpVector(sample_index, ...
                 rb_object.m_mixing_matrix, this.m_proj_mus_k_f_vectors);
         end
         
         % We need to add I to the operator mu_s K which we interpolate
        function [offset_op] = getAffineOffsetOperator(~, reduced_basis)
            n = size(reduced_basis, 2);
            offset_op = eye(n);
        end
         
        % View the reduced basis
        function viewBasis(this, reduced_basis, figure_id)
            display_basis_set(reduced_basis, figure_id, this.m_n_fine);
        end
       
        % Display a solution
        function viewResult(this, result_rb, result_fine, figure_id)
            display_two_solutions(result_rb, result_fine, figure_id, ...
                this.m_n_fine);
        end
        
        % Display a fine solution
        function viewFineResult(this, result_fine, figure_id)
            figure(figure_id);
            imagesc(reshape(result_fine, this.m_n_fine, this.m_n_fine));
            colorbar;
        end
        
        % Display a coarse solution
        function viewCoarseResult(this, result_coarse, figure_id)
            figure(figure_id);
            imagesc(reshape(result_coarse, this.m_n_coarse, this.m_n_coarse));
            colorbar;
        end
    end
end

function display_basis_set(uCs, figure_id, nC)
    figure(figure_id);
    nSk = size(uCs, 2);
    nSkrt = ceil(sqrt(nSk));
    m = max(max(uCs));
    for ind = 1:nSk
        subplot(nSkrt,nSkrt,ind);
        basis_vec = uCs(:,ind);
        imagesc(reshape(basis_vec, nC, nC));
        caxis manual
        caxis([0.0 m]);
        set(gca, 'xticklabel', []);
        set(gca, 'yticklabel', []);
        % colorbar;
    end
end

function display_two_solutions(uCApprox, uCTrue, figure_id, nC)
    figure(figure_id);
    
    m = max(max(uCApprox), max(uCTrue));
    
    subplot(1, 2, 1);
    imagesc(reshape(uCApprox, nC, nC));
    caxis manual
    caxis([0.0 m]);
    colorbar;
    title('Reduced Basis Approx.');
    
    subplot(1, 2, 2);
    imagesc(reshape(uCTrue, nC, nC));
    caxis manual
    caxis([0.0 m]);
    colorbar;
    title('Ground Truth');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute one coarse solution for the radiative transport
% equation.
%
% Parameters:
% force_func(x1, x2): the source term, force function
% mu_t_func(x1, x2, param): the transmission coefficient
% mu_s_func(x1, x2, param): the scattering coefficient
% omega: the sample which to use in the solve
% nC: grid size
%
% Return values:
% result: the solution of the equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = compute_coarse(force_func, mu_t_func, ...
    mu_s_func, omega, nC)

    % Get parameters from sample space
    musfuntmp = @(x1, x2) mu_s_func(x1, x2, omega);
    mutfuntmp = @(x1, x2) mu_t_func(x1, x2, omega);
    
    [N,xs,param,~] = setupparam(OracleRTE.g_quadrature_points, ...
        nC, OracleRTE.g_proxy, musfuntmp, mutfuntmp);
    
    fC2 = force_func(xs(1, :), xs(2, :));        
    fC2 = fC2(:);
    MC2 = Bfunnonsym(1:N, 1:N, param, true);
    
    % Solve coarse radiative transport system and save solution
    result2 = MC2 \ fC2 - fC2;

    result = result2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute one fine solution for the radiative transport
% equation.
%
% Parameters:
% force_func(x1, x2): the source term, force function
% mu_t_func(x1, x2, param): the transmission coefficient
% mu_s_func(x1, x2, param): the scattering coefficient
% omega: the sample which to use in the solve
% nF: grid size
%
% Return values:
% result: the solution of the equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [operator, diag_elements, fine_solution] = compute_fine(...
    force_func, mu_t_func, mu_s_func, omega, nF)
    
    % Specialize mu_s and mu_t function using parameter sample
    musfuntmp = @(x1, x2) mu_s_func(x1, x2, omega);
    mutfuntmp = @(x1, x2) mu_t_func(x1, x2, omega);

    [N,xs,param,mus] = setupparam(OracleRTE.g_quadrature_points, ...
        nF, OracleRTE.g_proxy, musfuntmp, mutfuntmp);

    diag_elements = spdiags(1./mus, 0, N, N);
    fC = force_func(xs(1,:), xs(2,:));        
    fC = fC(:);

    % Factor operator using HIF
    opts = struct('skip', OracleRTE.g_skip, 'symm', OracleRTE.g_symm, ...
        'verb', OracleRTE.g_verbose);
    Bfuntmp = @(i,j) Bfunnonsym(i, j, param, true);
    pxyfuntmp = @(x, slf, nbr, l, ctr) pxyfun(x, slf, nbr, l, ctr, param);
    operator2 = hifie2(Bfuntmp, xs, OracleRTE.g_occ, ...
        OracleRTE.g_rank_or_tol, pxyfuntmp, opts);

    % Solve system for fine_solution
    fine_solution2 = hifie_sv(operator2, fC) - fC;
    
    % Factorize mu_s K
    musKfuntmp = @(i,j) Bfunnonsym(i, j, param, false);
    operator = hifie2(musKfuntmp, xs, OracleRTE.g_occ, ...
        OracleRTE.g_rank_or_tol, pxyfuntmp, opts);
        
    fine_solution = fine_solution2;
end