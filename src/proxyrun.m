cd ../FLAM/
startup
cd ../RB-Framework/
[oracle, ss] = GaussianRTE();
oracle.setCoarseResolution(16);
oracle.setFineResolution(32);
rb = RBObject(oracle, ss);
rb.computeReducedBasis();
