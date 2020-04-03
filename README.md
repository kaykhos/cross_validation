# cross_validation
Implement cross fidelity validation, between qiskit device/simulator, simulator/simulator and parallel device/device

+ Can cross validate between circuits in a single job, or spesific circuits
  between mutiple jobs (multi-job default is all)

+ Makes Haar random single qubit rotations, and can spesify parallel jobs on a
  single circuit (not yet suported for cross fidelity comparison)

TODO: 
* Add batch processing: multiple CF in single job
* Add bootstrapping for error estimation
* Add sub-system fidelity (i.e. trace out some measurement results)
