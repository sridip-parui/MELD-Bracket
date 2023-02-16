#!/usr/bin/env python
# encoding: utf-8
# System
import glob
import argparse
import os
import itertools
import random
from pprint import pprint

# 3rd Party
import numpy as np
import mdtraj as md

# MELD
import simtk.openmm as mm
from meld.remd import ladder, adaptor, master_runner
from meld import comm, vault
from meld import system
from meld.system.restraints import TorsionRestraint
from meld.remd import reseed
from meld.system.restraints import LinearRamp,ConstantRamp

# Customize
import posescript
import ProteinAtomBank as atomBank




# Output the current working directory
P=os.getcwd()
# Number of replicas
N_REPLICAS = 18      

N_STEPS = 40000
BLOCK_SIZE = 100
N_LIGANDS = 6
def binding_restraint(s,R3=0.0,R4=0.0,scaler=None,Frac_Enforced=1.0,pairlist=None):
    pose_restraints = []
    dists = []
    active = int(np.floor(Frac_Enforced*len(pairlist)))
    k=0
    tmp_info = set()
    for ires,iname,jres,jname in pairlist:
        k+=1
        tmp_info.add((ires,jres))
        pose_restraints.append(
                s.restraints.create_restraint(
                        'distance',scaler,
                        ramp=LinearRamp(0,100,0,1),
                        r1= 0.0, r2=0.0,
                    r3=R3, r4=R4,
                        k = 250.,
                        atom_1_res_index=ires, atom_1_name=iname,
                        atom_2_res_index=jres, atom_2_name=jname,
                    use_pbc=True
                )
        )

    dists.append(s.restraints.create_restraint_group(pose_restraints, active))
    print(f"Restrained protein-ligand pair: {tmp_info}")

    return dists

def pulling_restraint(s,p_resid,p_atom_name,L_resid,L_atom_name,R2=0.0,R3=0.0,K=0.0,scaler=None):
    restraint = s.restraints.create_restraint(
                'distance', scaler,
                ramp=LinearRamp(0,100,0.5,1),
                r1=0.0, r2=R2 ,
            r3=R3, r4=R3+1.0,
            k=K,
                atom_1_res_index=L_resid , atom_1_name=L_atom_name,
                atom_2_res_index=p_resid , atom_2_name=p_atom_name,
            use_pbc=True
    )

    return restraint

def gen_state_templates(index,top,coord):
    print(index,top,coord)
    c = system.builder.load_amber_system('%s' %(top), '%s' %(coord))
    pos = c._coordinates
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    if index%(N_REPLICAS-1) == 0 and index != 0:
    	print("This is the box info:",c._box_vectors)
    return system.SystemState(pos, vel, alpha, energy, c._box_vectors)

def make_cartesian_collections(s, scaler, residues, ensemble=None,delta=None, restraint_atoms=None, k=None):
    cart = []
    assert delta is not None 
    assert k is not None 
    assert restraint_atoms is not None 

    if restraint_atoms is None or restraint_atoms not in ["backbone","all_heavy", "CA"]:
        raise Exception("Choose a proper cartesian restrain option: 'backbone' or 'all_heavy'.")

    if ensemble is None:
        raise Exception("You should provide structure for determing anchor position.")
    elif not isinstance(ensemble,list) or not isinstance(ensemble[0],system.System):
        raise Exception("You probably provide wrong ensemble.")
    else:
        if len(ensemble) == 1:
            print("Raise hand: Cartesian restraint position determined using single structure!!!")
        else: 
            print(f"Raise hand: Cartesian restraint position determined using ensemble of {len(ensemble)} structures!!!")

    print(f"Cartesian restraints applied on {restraint_atoms}")
    print(f"    We use {delta}A and {k}kJ/mol*nm-2")

    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)

    for i in residues:
        if restraint_atoms == "CA": # restrain only backbone atoms.
            rest_atoms = atomBank.AtomBankForProtein.CA()

        if restraint_atoms == "backbone": # restrain only backbone atoms.
            rest_atoms = atomBank.AtomBankForProtein.backbone()

        if restraint_atoms == "all_heavy": # restrain all heavy atoms of protein
            rest_atoms = atomBank.AtomBankForProtein.heavyAtoms(sequence[i])

            
        for b in rest_atoms:
            # For multiple structure cases, we assume the atom_index is the same.
            atom_index = s.index_of_atom(i,b) - 1
            # get (x,y,z) from all ensemble of structrues.
            pos = [ aSystem.coordinates[atom_index]/10. for aSystem in ensemble ]
            test_pos = random.sample(pos,1)[0]  # choose a random pos to check whether the average pos is a good choice. 
            pos = np.array(pos).mean(axis=0)   # average position 

            # check whether the average pos is a good choice.
            assert np.linalg.norm(test_pos-pos)<=delta and len(pos) == 3, "Something wrong with the input structure."
            rest = s.restraints.create_restraint('cartesian',scaler, res_index=i, atom_name=b,
                    x=pos[0], y=pos[1], z=pos[2],
                        delta=delta,force_const=k)
            cart.append(rest)
    print("Average backbone position is good.")
    assert len(cart) == len(residues) * len(rest_atoms), "Wrong number of cartesian restraints."
    return cart

def setup_system(prm,crd,lig_rest_files,lig_resid_idx):
    # load trajectory and get all the paramters
    trj = md.load('%s' %(crd[0]),top='%s' %(prm))
    prot_end= trj.topology.select("protein and name CA").shape[0] + 1 ## 1 residue more 

    # build simulation system
    s = system.builder.load_amber_system('%s' %(prm), '%s' %(crd[0]))

    # set temperature scaling
    s.temperature_scaler = system.ConstantTemperatureScaler(298.) ### Forget to change this one!!!!

    # Cartesian Restraint For Protein backbone
    const_scaler = s.restraints.create_scaler('constant')
    prot_delta = 0.5
    prot_residues = [135] ## Make sure this is the residue closest to the protein center of mass.

    # setup cartesian restraints
    # we are using average position from multiple structures
    ensembleOfSystem = [ 
            system.builder.load_amber_system(f"{prm}", f"{crdString}") for crdString in crd
            ]

    prot_Cartesian = make_cartesian_collections(
                s, const_scaler, prot_residues,
            ensemble=ensembleOfSystem,
            delta=prot_delta,
            restraint_atoms="backbone",
            k=250.    # kJ/mol*nm^(-2)
    )

    s.restraints.add_as_always_active_list(prot_Cartesian)

    

    # ligand unaggregate restraints
    # No need for this type of restraint for 2-competitor bracket
    R2 = 4.0  
    R3 = 20.0
    force_const = 250.
    lig_center_atoms = ["N1"]
    group_activate_num = 1
    lig_pairs = list(itertools.combinations(lig_resid_idx,2))
    lig_pair_push_rest = []

    for lig_i, lig_j in lig_pairs: # iterate every pair of ligands 
        for anchor in lig_center_atoms:
            lig_pair_push_rest.append(
                    pulling_restraint(s,
                        lig_i, anchor,
                        lig_j, anchor,
                        R2=R2, R3=R3,
                        K=force_const,
                        scaler=const_scaler
                    )
            )
        

    ###### create restraint group ##
    group_lig_pair_push_rest = s.restraints.create_restraint_group(
            lig_pair_push_rest,
            len(lig_pairs)*len(lig_center_atoms)
    )

    ###### create collection structure
    coll_lig_pair_push_rest = [group_lig_pair_push_rest]
    s.restraints.add_selectively_active_collection(coll_lig_pair_push_rest,group_activate_num)
    print(f"Number of unaggregating restraints: {len(lig_pair_push_rest)}")
    print(f"    We enforce: {len(lig_pairs)*len(lig_center_atoms)}")
    
    

    # ligand unbounded restraints
    R2 = 5.0     
    R3 = 20.0    # diagnol of the cubic box: 1.73*100/10=17.3nm
    force_const = 250.
    lig_push_atoms = ["C7","C8","N1"] ## Make sure those atoms appear in all ligands and contains in ligand's rigid scaffold. Pushing on atoms in floppy regions of ligand might introduce biases. 
    prot_push_atoms = (135,"CA") ## Make sure it is the closest residue to the protein center of mass. 
    group_activate_num = len(lig_resid_idx)-1 # N-1 unbounded ligand will be kept unbound during simulation.
    lig_resid_idx = lig_resid_idx


    lig_unbound_rests = []
    print(f"ligand residue: {lig_resid_idx}")
    for lig_resid in lig_resid_idx: ## N ligand residues
        one_lig_push_rest = []
        for anchor in lig_push_atoms:
            one_lig_push_rest.append(
                    pulling_restraint(s,
                        *prot_push_atoms,
                        lig_resid, anchor,
                        R2=R2, R3=R3,
                        K=force_const,
                        scaler=const_scaler
                    )
            )

        ## create restraint group ##
        group_lig_push_rest = s.restraints.create_restraint_group(one_lig_push_rest,len(lig_push_atoms))
        ## create collection structure
        lig_unbound_rests.append(group_lig_push_rest)

    assert len(lig_unbound_rests) == len(lig_resid_idx)
    print(f"Number of unbinding restraint groups: {len(lig_unbound_rests)}")
    print(f"    We enforce: {group_activate_num}")

    s.restraints.add_selectively_active_collection(lig_unbound_rests,group_activate_num)

    # setup ligand-protein distance restraints
    R3 = 0.7
    R4 = 0.9    
    enforce = 0.7 # original testing uses 0.9. To prevent underestimating entropic contribution, reduce it to 0.7!!!
    group_activate_num = 1 # only 1 ligand will be kept bound through out the simulation.
    binding_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.3, alpha_max=0.8, factor=4.0)

    lig_prot_rests = []
    ligand_rest_pairs = posescript.get_constraints(lig_rest_files)
    for ligand in ligand_rest_pairs:
          lig_prot_rests.extend(
                  binding_restraint(s,
                      R3=R3, R4=R4, 
                      Frac_Enforced=enforce,
                      scaler=binding_scaler,
                      pairlist=ligand
                  )
          )
    print(f"Number of binding restraint groups: {len(lig_prot_rests)}")
    print(f"    We enforce: {group_activate_num}")
    s.restraints.add_selectively_active_collection(lig_prot_rests, group_activate_num)

    # rest2_region 
    binding_site_residues_tmp = itertools.chain.from_iterable([lig_resid_idx]) # 1-based. ONLY two binders here.
    binding_site_residues = set([resid-1 for resid in binding_site_residues_tmp ]) # 0-based.
    if len(binding_site_residues) != N_LIGANDS:
        raise ValueError("The rest2 region is not set properly.")
    print(f"binding_site_residues (0-based): {list(binding_site_residues)}")

    # Simulation Options
    options                                  =        system.RunOptions(solvation='explicit')
    options.use_bigger_timestep              =        True                      # 4.5fs
    options.cutoff                           =        0.8						
    options.use_amap                         =        False						# ff14SB+TIP3P
    options.timesteps                        =        11112						# every 50ps ouput traj
    options.use_rest2                        =        True
    options.rest2_scaler                     =        system.REST2Scaler(298.,system.GeometricTemperatureScaler(0., 1, 298., 600.)) 
    options.rest2_region                     =        binding_site_residues
    store                                    =        vault.DataStore(s.n_atoms, N_REPLICAS, None, block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)
    l                                        =        ladder.NearestNeighborLadder(n_trials=48*48)
    policy                                   =        adaptor.AdaptationPolicy(2.0, 50, 50)
    a                                        =        adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy)
    remd_runner                              =        master_runner.MasterReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)
    # create and store the communicator
    c                                        =         comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)
    states                                   =         [gen_state_templates(i,prm,crd[i%len(crd)]) for i in range(N_REPLICAS)]
    store.save_states(states, 0)
    # save data_store
    store.save_data_store()

    return s.n_atoms

# MAIN CALL
parser = argparse.ArgumentParser()

parser.add_argument("--LigandRestraints",help="list of files containing ligand restraints",type=str,nargs='+',default = [])
parser.add_argument("--LigandResidueIdx",help="list of files containing ligand restraints",type=int,nargs='+',default = [])
parser.add_argument("--Top",help="Topology file",type=str,nargs='+',default = [])
parser.add_argument("--Coord",help="Well Equilibrated System's Coordinates",type=str,nargs='+',default = [])

args = parser.parse_args()
prm = args.Top[0]
crd = args.Coord
lig_rest_files = args.LigandRestraints
lig_resid_idx = args.LigandResidueIdx

print(crd)
setup_system(prm,crd,lig_rest_files,lig_resid_idx)

