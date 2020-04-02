#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02/04/2020

@author: kiran

Functions work on a single anzatz circuit and append single qubit Haar
    random unitaries at the end. 
    + Can deal with two of the same size circuits running on singe device
      - Must be same size subcircuits
      - Jobs MUST be grouped at the QuantumRegister level
    + When one or more register grou"""



import qiskit as qk
import numpy as np
import scipy as sp
import copy
pi = np.pi

#%%
try:
    provider_free
except:
    provider_free = qk.IBMQ.load_account()

def _random_measurement_helper(circ, 
                               seeds, 
                               qregs_block=0):
    # get the registers to measure and add classical bits
    qregs = circ.qregs[qregs_block]
    nb_qubits = qregs.size
    cbits = qk.ClassicalRegister(nb_qubits, 'cl_'+qregs.name)
    circ.add_register(cbits)

    
    # generate Haar random (with known seeds)
   # if type(seeds) is not list:
   #     seeds = np.random.randint(0, 2**32, nb_qubits)
    u_random = [qk.quantum_info.random_unitary(2, seed=seeds[ii]) 
                for ii in range(nb_qubits)]
    
    # append unitaries and measure to right classical registers
    for ii in range(nb_qubits):
        circ.append(u_random[ii], [qregs[ii]])
    circ.measure(qregs, cbits)
    return circ
    
def random_measurement(circ, 
                       qregs_block_list=[0],
                       nb_random=5, 
                       seed=10):
    # input formatting 
    if type(qregs_block_list) == int:
        qregs_block_list = [qregs_block_list]
    np.random.seed(seed=seed)
    circ_list = []
    
    # Run over different numbers of circuits
    for ii in range(nb_random):
        # Fix circuit and set random seed (so each block has SAME unitaries)
        this_circ = copy.deepcopy(circ)
        nb_qubits = this_circ.qregs[qregs_block_list[0]].size
        seeds = np.random.randint(0, 2**32, nb_qubits)
        for qregs_block in qregs_block_list:
            this_circ = _random_measurement_helper(this_circ, 
                                                   qregs_block=qregs_block, 
                                                   seeds=seeds)
        circ_list.append(this_circ)
    return circ_list







def cross_fidelity_from_dual_res(results1, results2, d=2): # d=2 for qubit, 3 for qutrit
    """Impliments Eq. 2 from Cross-Platform Verification of Intermediate Scale 
        Quantum Devices
        Basic checks at the beginning to ensure inputs are valid (not exhaustive)
        
        NEEDED: Maybe update to work with results.results lists instead to 
                support multiple runs
        """
    # very basic checks on inputs
    nb_qubits1 = len(list(results1.get_counts(0).keys())[0])
    nb_qubits2 = len(list(results2.get_counts(0).keys())[0])
    nb_u1 = len(results1.results)
    nb_u2 = len(results2.results)
    if nb_qubits1 == nb_qubits2 and nb_u1 == nb_u2:
        nb_qubits = nb_qubits1
    else:
        return 'Error - results dont match dims for this method'
    
    # Generae all possible keys may, or may not have measurement realisation
    keys  = _gen_all_possilbe_keys(nb_qubits)
    nb_hilbert_dims = d**nb_qubits
    F = 0
    
    # Double sum in Eq. 2, coefficient times cross corrslation
    for k1 in keys:
        for k2 in keys: 
            hamming_distance = nb_qubits*sp.spatial.distance.hamming(list(k1), list(k2))
            hamming_distance = int(hamming_distance)
            coeff = (d)**nb_hilbert_dims * (-d)**(-hamming_distance)
            F += coeff * _correlation(results1, results2, k1, k2)
            
    return F


def _correlation(results1, results2, key1, key2):
    """Computs the corelation between at two measurement outcomes, across 
    all instances of the measurement results
    
    NEEDED: I should update this to handel results list e.g. (results.results) 
            so runs across mutiple jobs are suported """
    # get number of random unitaries 
    nb_u = len(results1.results)
    correlation = 0
    for ii in range(nb_u): # for each unitary
        counts1 = results1.get_counts(ii)
        counts2 = results2.get_counts(ii)
        
        keys1 = list(counts1.keys())
        keys2 = list(counts2.keys())
        
        norm = sum(counts1.values()) *  sum(counts2.values())
        
        # basic error handeling if not all measurements are realized
        if key1 in keys1 and key2 in keys2: 
            correlation += counts1[key1] * counts2[key2]
    correlation = correlation / nb_u / norm
    return correlation


def _gen_all_possilbe_keys(nb_qubits):
    vec = [bin(ii)[2:].zfill(nb_qubits) for ii in range(2**nb_qubits)]
    vec = [str(ii) for ii in vec]
    return vec


















if __name__ == '_main__':
    # quick checks to see if it wokrs
    backend = provider_free.get_backend('ibmq_qasm_simulator')
    instance = qk.aqua.QuantumInstance(backend, shots=2048, optimization_level=3)
    
    
    # Define simple circuit
    circ = qk.QuantumCircuit(qk.QuantumRegister(2, 'circ1'))
    circ.add_register(qk.QuantumRegister(2, 'circ2'))
    circ.ry(1,0)
    circ.ry(2,1)
    circ.barrier()
    circ.cx(0,1)
    circ.cx(2,3)
    circ.barrier()
    circ.draw()
    
    # Check measurements applied to different parts of the circuit, with SAME random U
    circs1 = random_measurement(circ, [0], nb_random=2)
    circs1[0].decompose().draw()
    
    circs2 = random_measurement(circ, [1], nb_random=2)
    circs2[0].decompose().draw()

    results1 = instance.execute(circs1)
    results2 = instance.execute(circs2)


    print(cross_fidelity_from_dual_res(results1, results2))

    keys =  _gen_all_possilbe_keys(2)
    check_norm = 0
    for k1 in keys:
        for k2 in keys:
            check_norm += _correlation(results1, results2, k1, k2)
    print(check_norm)







