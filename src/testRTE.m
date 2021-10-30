% Construct an oracle and samples space for the Radiative Transport Equation
[oracle, sample_space] = BenchmarkRTE('test');

% oracle is an object which can be used to solve fine and coarse versions
% of a given problem. The sample_space is a matrix where the columns
% correspond to parameter values for specific instances of the problem we
% are looking to solve. See GaussianRTE for more information about how to
% set up an oracle and a sample space.

% We can set the coarse and fine resolution through the problem oracle
oracle.setCoarseResolution(16);
oracle.setFineResolution(32);

% We now create an RBOBject. This object is responsible for constructing
% the reduced basis using the oracle.
rb_obj = RBObject(oracle, sample_space);

% Change the R-factor cutoff for selecting skeleton samples
rb_obj.m_selection_epsilon = 1E-2;
% Number of rows to sample from operators
rb_obj.m_n_operator_samples = 10; 
% Whether or not to use the operators to select additional skeletons
rb_obj.b_enable_additional_operator_skeletons = true;
rb_obj.b_enable_additional_skeleton_solutions = true;

% Compute our reduced basis
rb_obj.computeReducedBasis();

% View the reduced basis, parameter is the figure id
rb_obj.viewReducedBasis(1);
% View a random result, parameter is the figure id
rb_obj.viewRandomResult(2);

% Use this to compute a specific RB solution, parameter is the sample 
% index in the sample space
[~] = rb_obj.computeRBSolution(1);

% To compute error diagonostics
subsample_factor = 0.1; % Only sample this fraction of the true solutions
rb_obj.computeDebugData(subsample_factor);
rb_obj.outputDiagnostics();