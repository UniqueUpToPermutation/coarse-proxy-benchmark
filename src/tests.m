% Runs a test to make sure that everything is working correctly

fprintf('Testing Radiative Transport Equation...\n\n');

[oracle, ss] = GaussianRTE();
oracle.setCoarseResolution(8);
oracle.setFineResolution(16);
rbobj = RBObject(oracle, ss);
rbobj.b_enable_gramm_schmidt_for_operator_selection = true;
rbobj.computeReducedBasis();
rbobj.computeDebugData();
rbobj.outputDiagnostics();

fprintf('\n\nTesting Boundary Integral Equation..\n\n');

[oracle, ss] = SingularBIE();
rbobj = RBObject(oracle, ss);
rbobj.b_enable_gramm_schmidt_for_operator_selection = true;
rbobj.computeReducedBasis();
rbobj.computeDebugData();
rbobj.outputDiagnostics();