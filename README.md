# cross_validation
Implement cross fidelity validation, between qiskit device/simulator, simulator/simulator and parallel device/device

+ Can cross validate between circuits in a single job, or spesific circuits
  between mutiple jobs (multi-job default is all)

+ Makes Haar random single qubit rotations, and can spesify parallel jobs on a
  single circuit (not yet suported for cross fidelity comparison)

TODO: 
* Add bootstrapping for error estimation
* Add sub-system fidelity (i.e. trace out some measurement results)


# Examples

# Define simple circuit
    backend = provider_free.get_backend('ibmq_qasm_simulator')
    instance = qk.aqua.QuantumInstance(backend, shots=2048, optimization_level=3)
    
    
    # Create two circuits that have cross F = 0
    circ1 = qk.QuantumCircuit(qk.QuantumRegister(2, 'regs_1'), name='circ1')
    circ1.rx(pi,0)
    circ1.rx(pi,1)
    circ1.barrier()
    
    circ2 = qk.QuantumCircuit(qk.QuantumRegister(2, 'regs_1'), name='circ2')
    circ2.barrier()

    
    # Append random measurements in two ways
    #   1 each circuit independently
    circ1m = append_random_unitaries(circ1, seed=10)
    circ2m = append_random_unitaries(circ2, seed=10)
    print(circ1m[0].name)
    print(circ1m[0].decompose())
    print(circ2m[0].name)
    print(circ2m[0].decompose())
    
    #   2 Or combine them into a single job (here same circuits)
    circ_all_m = append_random_unitaries([circ1, circ2], seed=100)
    
    # Also ensures right params at right points
    print(circ_all_m[0].name)
    print(circ_all_m[0].decompose())
    print(circ_all_m[5].name)
    print(circ_all_m[5].decompose())
    
    # note circuits now have different names:
    _print_names(circ1m);print()
    _print_names(circ2m);print()
    _print_names(circ_all_m);print()
    
    
    # run simulations
    try:
        res1
    except:
        res1 = instance.execute(circ1m)
        res2 = instance.execute(circ2m)
        res_all = instance.execute(circ_all_m)
        
    # If you used append_random_unitaries() to create the circuts, you don't need
    #   to worry about options (as long as same seed)
    
    # For different results objects
    F = cross_fidelity([res1, res1])
    print("CF from 2 qobj's (circ1 agains itsself): F = {}".format(F))
    
    # For the same res object, define circs be name
    F = cross_fidelity(res_all, prefix1='circ1', prefix2='circ2')
    print('CF from 1 qobj different namded circs in same job: F = {}'.format(F))
    
    # Or by place in list
    F = cross_fidelity(res_all, unitary_block1=range(5), unitary_block2=range(5,10))
    print("CF from 1 obj by list instead of names: F = {}".format(F))
    
    
    # Compare res object with saved result.to_dict()
    res1_dict = res1.to_dict()
    F = cross_fidelity([res1_dict, res2])
    print("CF between dict(1) and res oobject(2): F = {}".format(F))
    
    # Can also compare any and all parts:
    # compare res1 itsself from circs_all
    F = cross_fidelity([res1_dict, res_all], prefix2='circ1')
    print("CF between circ1 in res_all and circ1 in res1: F = {}".format(F))
