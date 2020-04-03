#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02/04/2020

@author: kiran

Functions work on a single anzatz circuit and append single qubit Haar
    random unitaries at the end. 
    NEEDED: Have to add functionality to compare only subset of measurements!
    + Can deal with two of the same size circuits running on singe device
      - Must be same size subcircuits
      - Jobs MUST be grouped at the QuantumRegister level
    + When one or more register group
    """

import qiskit as qk
import numpy as np
import scipy as sp
import copy
import pdb
pi = np.pi

#%%
try:
    provider_free
except:
    provider_free = qk.IBMQ.load_account()

# Annoying to set this, but I don't want to hard code anything
NB_QUBITS_HAVE_DIM_2 = 2




#%% Appending random unitaries to circuit
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def append_random_unitaries(circ, 
                       qregs_block_list=[0],
                       nb_random=5, 
                       seed=10):
    """ Creates a list of n=nb_random circuits with Haar random unitaries. If
        mutiple qregs_block_list are passed, each block has the same unitaries
        - DOES NOT modify input circuit
        TO DO: Allow passing of qregs names in the list"""
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


def _random_measurement_helper(circ, 
                               seeds, 
                               qregs_block=0):
    """Appends (seeded) random unitaries for a circuit
        - Will modify input circ if the input circ is passes as object
        - Not made to be interacted with directly 
        - Adds measurement registers with same name as qregs_block
        - Assumes circuit construction is created with single /mutiple blocks"""
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



#%% Implimenting cross fidelity between measurement results
    # Handels results.to_dict() inputs OR direct result objctes
    # In Single obj is passed, user MUSH spesify which jobs to look at via the 
    #   unitary_blocks lists
    # Unfortinuately nb_qubits is nessassary (for now) to deal with the bullshit hex keys
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def cross_fidelity(results_in,
                   unitary_block1=None,
                   unitary_block2=None):
    """ Accepts as input results_in is a results object, a results.to_dict()
        or a list of EXACTLY two results objs in either format
        
        TO DO: add ability to pass miltiple blocks in 
        TO DO: add ability to pass list of measurement outcomes for spesific unitary"""
    results1, results2 = _gen_results_to_compare(results_in,
                                                 unitary_block1=unitary_block1,
                                                 unitary_block2=unitary_block2)       
    cross_overlap = density_matrix_overlap(results1, results2)
    first_overlap = density_matrix_overlap(results1, results1)
    secon_overlap = density_matrix_overlap(results2, results2)
    F = cross_overlap / max(first_overlap, secon_overlap)
    return F


def density_matrix_overlap(results1, results2): # d=2 for qubit, 3 for qutrit
    """Impliments Eq. 2 from Cross-Platform Verification of Intermediate Scale 
        Quantum Devices
        Basic checks at the beginning to ensure inputs are valid (not exhaustive)
        """
    # very basic checks on inputs
    nb_qubits1 = len(list(results1[0].keys())[0])
    nb_qubits2 = len(list(results2[0].keys())[0])
    nb_u1 = len(results1)
    nb_u2 = len(results2)
    if nb_qubits1 == nb_qubits2 and nb_u1 == nb_u2:
        nb_qubits = nb_qubits1
    else:
        assert False, 'Error - results dont match dims for this method'
    
    # Generae all possible keys may, or may not have measurement realisation
    keys1 = _gen_keys_from_lists(results1)
    keys2 = _gen_keys_from_lists(results2)
    keys = list(set.intersection(keys1, keys2))
    nb_hilbert_dims = NB_QUBITS_HAVE_DIM_2**nb_qubits
    
    # Double sum in Eq. 2, coefficient times cross corrslation paper
    Trace = 0
    for k1 in keys:
        for k2 in keys: 
            hamming_distance = int(nb_qubits*sp.spatial.distance.hamming(list(k1), list(k2)))
            coeff = nb_hilbert_dims * (-NB_QUBITS_HAVE_DIM_2)**(-hamming_distance)
            Trace += coeff * _correlation(results1, results2, k1, k2)
            
    return Trace


def _correlation(results1, results2, key1, key2):
    """Computs the corelation between at two measurement outcomes, assumes all inputs
        all contribute to correlation"""
    # get number of random U matricies
    if type(results1) != list: results1 = _gen_measurement_list(results1)
    if type(results2) != list: results2 = _gen_measurement_list(results2)

    nb_u = len(results1)
    correlation = 0
    for ii in range(nb_u): # for each unitary
        
        # Get counts and keys, and number of shots (robust to key errors)
        keys1 = list(results1[ii].keys())
        keys2 = list(results2[ii].keys())
        norm = sum(results1[ii].values()) *  sum(results2[ii].values())
        # basic error handeling if not all measurements are realized
        if key1 in keys1 and key2 in keys2: 
            correlation += results1[ii][key1] *results2[ii][key2]
    
    # Normalize ensemble mean to number of input shots, and number of unitariess
    correlation = correlation / nb_u / norm
    return correlation


# -----------------------------------------------------------------------------
# --------- B/S functions that deal with results.to_dict() keys----------------
# -----------------------------------------------------------------------------
    
def _gen_results_to_compare(results_in,
                            unitary_block1=None, 
                            unitary_block2=None):
    """ Helper function inputs results object OR dictionary, and returns a subset
        of measurement results corresponding to the unitary blocks
        - results_in can be two different objects or one single object
        - if results in is a single object, unitary blocks must be spesified
        """
    if type(results_in) == list:
        results1 = results_in[0]
        results2 = results_in[1]
        results1 = _gen_measurement_list(results1)
        results2 = _gen_measurement_list(results2)
        if unitary_block1 != None:
            results1 = [results1[ii] for ii in unitary_block1]
        if unitary_block2 != None:
            results2 = [results2[ii] for ii in unitary_block2]
    else:
        assert unitary_block1 != None, " Unitary blocks MUST be spesified for a single input"
        assert unitary_block2 != None, " Unitary blocks MUST be spesified for a single input"
        temp_res = _gen_measurement_list(results_in)
        results1 = [temp_res[ii] for ii in unitary_block1]
        results2 = [temp_res[ii] for ii in unitary_block2]
    return results1, results2
    


def _gen_measurement_list(results_in):
    if type(results_in) == dict:
        results_in = qk.result.result.Result.from_dict(results_in)
    ls = []
    for ii in range(len(results_in.results)):
        ls.append(results_in.get_counts(ii))
    return ls
    


def _gen_all_possilbe_keys(nb_qubits):
    """ Returns a list of all possible measurement outcomes for an nb_qubit 
        experiment
        - E.g for 2 qubits, it returns: ['00', '01', '10', '11']"""
    vec = [bin(ii)[2:].zfill(nb_qubits) for ii in range(2**nb_qubits)]
    vec = [str(ii) for ii in vec]
    return vec



def _gen_keys_from_lists(list_in):
    keys = set({})
    for rr in list_in:
        keys = keys.union(set(rr.keys()))
    return keys


#%% if main

if __name__ == '__main__':
    # quick checks to see if it all wokrs
    backend = provider_free.get_backend('ibmq_qasm_simulator')
    instance = qk.aqua.QuantumInstance(backend, shots=2048, optimization_level=3)
    
    
    # Define simple circuit
    circ = qk.QuantumCircuit(qk.QuantumRegister(2, 'circ1'))
    circ.add_register(qk.QuantumRegister(2, 'circ2'))
    circ.rx(pi,0)
    circ.rx(pi,1)
    circ.barrier()
    circ.draw()
    
    # Check measurements applied to different parts of the circuit, with SAME random U
    circs1 = append_random_unitaries(circ, qregs_block_list=[0], nb_random=10)
    print(circs1[0].decompose())
    
    circs2 = append_random_unitaries(circ, qregs_block_list=[1], nb_random=10)
    print(circs2[0].decompose())

    # Run simulations
    results1 = instance.execute(circs1)
    results2 = instance.execute(circs2)
    
    
    # Checking helper functions on single set of measurement results
    res_list_1, res_list_2 = _gen_results_to_compare([results1, results1], [1,2,3], [4,5,6])
    check_norm = 0
    for k1 in _gen_all_possilbe_keys(2):
        for k2 in _gen_all_possilbe_keys(2):
            check_norm += _correlation(res_list_1, res_list_2, k1, k2)
    
    # Is everything working 
    assert (check_norm) == 1, "Double sum over cross correlation should be 1"
    assert cross_fidelity([results1, results1]) == 1, "Cross fidelity of same object should be 1"
    assert round(5*cross_fidelity([results1, results2])) == 0, "Cross fidelity between orthogonal circuits should be < 5%"
    
    
    # Look at much larger circuit with redundent measurements
    qu = qk.QuantumRegister(5)
    cl = qk.ClassicalRegister(5)
    circ = qk.QuantumCircuit(qu, cl)
    print(circ)
    circ0 = append_random_unitaries(circ, nb_random=5)
    circ1 = append_random_unitaries(circ, nb_random=5)
    print(circ0[0])
    print(circ1[0])
    circ = circ0 + circ1
    rr = instance.execute(circ)
    
    print('CF from qobj: F = {}'.format(cross_fidelity(rr, unitary_block1=[0,1,2,3,4], unitary_block2=[5,6,7,8,9])))
    rr = rr.to_dict()
    print('CF from dict: F = {}'.format(cross_fidelity(rr, unitary_block1=[0,1,2,3,4], unitary_block2=[5,6,7,8,9])))
    
          
    
    
