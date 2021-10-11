classdef (Abstract) ProblemOracle < handle

    properties (Access = protected)
        % The coarse resolution
        m_n_coarse
        % The fine resolution
        m_n_fine
    end
    
    methods
        % Class constructor, takes in coarse and fine problem resolutions
        function obj = ProblemOracle(n_coarse, n_fine)
            obj.setCoarseResolution(n_coarse);
            obj.setFineResolution(n_fine);
        end
        
        % Set coarse resolution
        function setCoarseResolution(this, n_coarse)
            this.m_n_coarse = n_coarse;
        end
        
        % Set fine resolution
        function setFineResolution(this, n_fine)
           this.m_n_fine = n_fine; 
        end
        
        % Get coarse resolution
        function [n] = getCoarseResolution(this)
            n = this.m_n_coarse;
        end
        
        % Set fine resolution
        function [n] = getFineResolution(this)
            n = this.m_n_fine;
        end
        
        % Get an offset operator for the reduced operators we form
        function [offset_op] = getAffineOffsetOperator(~, reduced_basis)
            n = size(reduced_basis, 2);
            offset_op = zeros(n);
        end
    end
      
    % Abstract methods to be implemented by subclasses
    methods (Abstract)
             
        % Get the vector size of coarse solutions 
        [n] = getCoarseSolutionSize(this)
        
        % Get the vector size of fine solutions 
        [n] = getFineSolutionSize(this)
        
        % Solve the coarse system for sample omega
        % Returns solution vector
        [solution] = solveCoarse(this, omega)
        
        % Solve the fine system for sample omega
        % Returns vector solution, left hand operator, and auxiliary information
        [solution, loperator, aux] = solveFine(this, omega)
        
        % Project an operator into the reduced basis
        [proj_operator] = projectFineOperator(this, reduced_basis, operator)
        
        % Sample important sections of all operators in the sample space
        [samples] = sampleFineOperators(this, sample_space, n_samples)
        
        % Compute auxiliary rb data needed to assemble the right hand side of
        % equations.
        computeAuxData(this, rb_object)
        
        % Assemble the right hand of a reduced basis solve
        [proj_f] = assembleRightHand(this, rb_object, omega, sample_index)
        
        % View the reduced basis
        viewBasis(this, reduced_basis, figure_id)
       
        % Display a fine solution
        viewResult(this, result_rb, result_fine, figure_id)
    end
end

