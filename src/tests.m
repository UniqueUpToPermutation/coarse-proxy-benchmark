% Runs a test to make sure that everything is working correctly

disp('Testing Radiative Transport Eqation...');

[oracle, ss] = GaussianRTE();
oracle.setCoarseResolution(8);
oracle.setFineResolution(16);
rbobj = RBObject(oracle, ss);
rbobj.b_enable_gramm_schmidt_for_operator_selection = true;
rbobj.computeReducedBasis();
rbobj.computeDebugData();
rbobj.outputDiagnostics();

disp('Testing Boundary Integral Equation..');

[oracle, ss] = SingularBIE();
rbobj = RBObject(oracle, ss);
rbobj.b_enable_gramm_schmidt_for_operator_selection = true;
rbobj.computeReducedBasis();
rbobj.computeDebugData();
rbobj.outputDiagnostics();