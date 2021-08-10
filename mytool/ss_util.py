#!/usr/bin/env python
##### CODE BY YOUNG JAE CHOI #####

from ase import Atoms, Atom
import random
from ase.build import make_supercell
import numpy as np

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

def slice2str(s):
    string = ''
    if s.start is not None:
        string += str(s.start)
    string += ':'
    if s.stop is not None:
        string += str(s.stop)
    string += ':'
    if s.step is not None:
        string += str(s.step)
    return string

def parse_slice(s):
    # Slice format
    if ':' in s:
        a = [int(e) if e.strip() else None for e in s.split(":")]
    # Int format
    else:
        if int(s) == -1:
            a = [-1, None]
        else:
            a = [int(s), int(s)+1]
    return slice(*a)

def get_number_of_lines(f_obj):
    f_obj.seek(0)
    for i, l in enumerate(f_obj):
        pass
    f_obj.seek(0)
    return i + 1

def str_slice_to_list(
    str_slice,
    obj_length=None,
    ):
    """
    str_slice (str)  = String type slice. e.g. 3:-1:2
    obj_length (int) = Length of object that will use this slice object. Set to None if unknown.
    """

    #
    str_slice = str(str_slice)
    if obj_length:
        obj_length = int(obj_length)

    ## Main
    if ':' in str_slice:
        slice_list = [int(e) if e.strip() else None for e in str_slice.split(":")]
        if len(slice_list) == 2:
            slice_list.append(None)
        elif len(slice_list) != 3:
            raise ValueError('String slice option is unreadable. --> {}'.format(str_slice))
    elif int(str_slice) == -1:
        slice_list = [-1, None, None]
    else:
        slice_list = [int(str_slice), int(str_slice)+1, None]

    ## Post-process : To achieve unity.
    if obj_length:
        slice_list = slice(*slice_list).indices(obj_length)
    else:
        # 0
        if slice_list[0] == None:
            slice_list[0] = 0
        # 2
        if slice_list[2] == None:
            slice_list[2] = 1
    return slice_list

def rectify_curve(
    curve,
    rectify_cut,
    ):
    iter = True
    while True:
        test = curve[1:] - curve[:-1]
        peak_bool = np.array(list(test[:,1] > (-1 * rectify_cut)) + [True], dtype=np.bool)
        if False not in peak_bool:
            break
        curve = curve[peak_bool]
    return curve

def get_chem_ind_arr(
    chemical_symbols,
    ):
    chem_arr = np.array(chemical_symbols)
    unique_chem = np.unique(chem_arr)
    ind = np.arange(len(chem_arr))
    chem_ind_arr = []
    for chem in unique_chem:
        chem_ind_arr.append(ind[chem_arr==chem])
    return unique_chem, np.array(chem_ind_arr)

def read_txt_matrix(
    file_name,
    start_line=0,
    end_line=-1,
    ):
    mat = []
    with open(file_name, 'r') as f:
        for i in range(start_line):
            f.readline()
        if end_line != -1:
            for i in range(end_line - start_line + 1):
                line = f.readline()
                seg = line.split()
                mat.append(seg)
        else:
            while True:
                line = f.readline()
                if not line: break
                seg = line.split()
                mat.append(seg)
    return mat

def nospace(string): return string.replace(' ', '_')
def E_fromAlist(alist):
    """ Return potential energy array """
    energies = []
    for atoms in alist:
        energies.append(atoms.get_potential_energy())
    return np.array(energies)
def Eperatom_fromAlist(alist):
    """ Return potential energy per atom array """
    energies = []
    for atoms in alist:
        energies.append(atoms.get_potential_energy()/len(atoms))
    return np.array(energies)
def F_fromAlist(alist):
    """ Return force array array """
    forces = []
    for atoms in alist:
        forces.append(atoms.get_forces())
    return np.array(forces)
def v_fromAlist(alist):
    """ Return velocity array array """
    velocities = []
    for atoms in alist:
        velocities.append(atoms.get_velocities())
    return np.array(velocities)
def P_fromAlist(alist):
    """ Return momentum array array """
    momenta = []
    for atoms in alist:
        momenta.append(atoms.get_momenta())
    return np.array(momenta)
def S_fromAlist(alist):
    """ Return stress array (xx, yy, zz, yz, zx, xy) array """
    stresses = []
    for atoms in alist:
        stresses.append(atoms.get_stress())
    return np.array(stresses)
def AtomicK_fromAlist(alist):
    """ Return atomic kinetic energy array array """
    kinetic_e=[]
    for atoms in alist:
        kinetic_e.append(0.5 * atoms.get_masses() * np.square(np.linalg.norm(atoms.get_velocities(), axis=-1)))
    return np.array(kinetic_e)

def column2np(
    txtfile,
    column,
    start_line = 1,
    end_line   = -1,
    interval   = 1,
    ):
    txt = open(txtfile, "r")
    lines = txt.readlines()
    if end_line == -1:
        end_line = len(lines)
    nums = []
    for l in range(start_line, end_line, interval):
        nums.append(float(lines[l].split()[column-1]))
    array = np.array(nums)
    return array

def RanPoAtoms(cut_off_radius,
               symbols=None,
               positions=None, numbers=None,
               tags=None, momenta=None, masses=None,
               magmoms=None, charges=None,
               scaled_positions=None,
               cell=None, pbc=None, celldisp=None,
               constraint=None,
               calculator=None,
               info=None):
    if positions is not None:
        print("\npositions must not be given\n")
        exit(1)
    if scaled_positions is not None:
        print("\nscaled_positions must not be given\n")
        exit(1)
    else:
        atoms = Atoms(symbols=symbols,
                      positions=positions, numbers=numbers,
                      tags=tags, momenta=momenta, masses=masses,
                      magmoms=magmoms, charges=charges,
                      scaled_positions=None,
                      cell=cell, pbc=pbc, celldisp=celldisp,
                      constraint=constraint,
                      calculator=calculator,
                      info=info)    
        l = 0
        while True:
            l+=1
            print("trying step :: "+str(l))
            scaled_posis = []
            for i in range(len(atoms)):
                scaled_posi = []
                for j in range(3):
                    scaled_posi.append(random.random())
                scaled_posis.append(scaled_posi)
            atoms.set_scaled_positions(scaled_posis)
            supercell = make_supercell(atoms,[[2,0,0],[0,2,0],[0,0,2]])
            dist = supercell.get_all_distances()
            coll = []
            for i in range(len(supercell)):
                for j in range(len(supercell)):
                    if i is not j:
                        coll.append(dist[i][j])
            if min(coll) > cut_off_radius:
                break
        return atoms

def RanPoAtoms_2(cut_off_radius,
                random_degree,
                symbols=None,
                positions=None, numbers=None,
                tags=None, momenta=None, masses=None,
                magmoms=None, charges=None,
                scaled_positions=None,
                cell=None, pbc=None, celldisp=None,
                constraint=None,
                calculator=None,
                info=None):
    if positions is None:
        if scaled_positions is None:
            print("\nNo backbone structure is given.\n")
            exit(1)
    else:
        atoms = Atoms(symbols=symbols,
                      positions=positions, numbers=numbers,
                      tags=tags, momenta=momenta, masses=masses,
                      magmoms=magmoms, charges=charges,
                      scaled_positions=scaled_positions,
                      cell=cell, pbc=pbc, celldisp=celldisp,
                      constraint=constraint,
                      calculator=calculator,
                      info=info)    

        ############### shuffle positions ################
        array_scaled_positions = atoms.get_scaled_positions()
        shuffled_scaled_posis = array_scaled_positions.tolist()
        random.shuffle(shuffled_scaled_posis)
        atoms.set_scaled_positions(shuffled_scaled_posis)
        
        ############### get local random distribution radius ################
        supercell = make_supercell(atoms,[[2,0,0],[0,2,0],[0,0,2]])
        dist = supercell.get_all_distances()
        coll = []
        for i in range(len(supercell)):
            for j in range(len(supercell)):
                if i is not j:
                    coll.append(dist[i][j])
        ran_radi = min(coll)

        ############### shuffled position list ################
        array_shuffled_posis = atoms.get_positions()

        l = 0
        while True:
            l+=1
            print("trying step :: "+str(l))
            shuffled_posis = array_shuffled_posis.tolist()
            for i in range(len(atoms)):
                for j in range(3):
                    shuffled_posis[i][j]+=((random.random()-0.5)*2*ran_radi*random_degree)
            tmp = atoms
            tmp.set_positions(shuffled_posis)
            supercell = make_supercell(tmp,[[2,0,0],[0,2,0],[0,0,2]])
            dist = supercell.get_all_distances()
            coll = []
            for i in range(len(supercell)):
                for j in range(len(supercell)):
                    if i is not j:
                        coll.append(dist[i][j])
            if min(coll) > cut_off_radius:
                break
        atoms = tmp
        return atoms

def count2list(dict_in):
    list_out = []
    for i in range(len(dict_in.keys())):
        key = list(dict_in.keys())[i]
        num = 0
        itera = dict_in[key]
        for j in range(itera):
            list_out.append(key)
            num += 1
    return list_out

def list2count(list_inp):
    """
    input: list.
    output: dict of counts of elements.
    """
    keys = list(set(list_inp))
    dict_out = dict()
    for i in keys:
        dict_out[i] = 0
    for i in list_inp:
        dict_out[i] += 1
    return dict_out

def list2numlist(list_inp):
    """
    input: list of any element types.
    output: list of integers. Numbers are assigned sequently to elements of first appearing.
    ex)
    inp: ['b', 'a', 'a', 'd', 'a', 'c']
    out: [ 0 ,  1 ,  1 ,  2 ,  1 ,  3 ]
    """
    numlist = []
    keys=dict()
    i=0
    for ele in list_inp:
        if ele in keys:
            numlist.append(keys[ele])
        else:
            keys[ele]=i
            numlist.append(keys[ele])
            i+=1
    return numlist

def ind_dict(list_inp):
    items = np.unique(list_inp)
    ind_dict = {}
    for item in items:
        ind_dict[item] = np.where(np.array(list_inp) == item)[0]
    return ind_dict

def ordered_unique(in_list):
    tmp = set()
    return [x for x in in_list if not (x in tmp or tmp.add(x))]

def covalent_expect(input):
    """ Returns covalent bond expectation value of system 
    input : dict or list
    e.g. {'Si':1, 'H':2} or ['Si', 'H', 'H']

    """

    from ase.atoms import symbols2numbers as s2n
    from ase.data import covalent_radii
    
    if isinstance(input, list):
        num_spec_dict = list2count(input)
    elif isinstance(input, dict):
        num_spec_dict = input
    else:
        raise TypeError("input is not list nor dict")
    
    tot_num = np.sum(list(num_spec_dict.values()))
    r_sum = 0
    for key, value in num_spec_dict.items():
        r_sum += covalent_radii[s2n([key])[0]] * value
    expect_value = r_sum / tot_num
    return expect_value

def random_atoms_gen(
    backbone,
    num_spec_dict = None,
    fix_ind_dict  = None,
    pin_the_fixed = False,
    cutoff_radi   = None,
    cutoff_frac   = None,
    random_radi   = None,
    random_frac   = None,
    strain        = None,
    strain_ratio  = [1.,1.,1.],
    vacuum        = None,
    vacuum_ratio  = None,
    max_trial_sec = 5,
    log           = True,
    ):
    """ 
    RAG generates randomly positioned image based on backbone structure.

    backbone : An ASE atoms object
        All atomic positions of the backbone structure will be the lattice sites.

    num_spec_dict : dict or None
        Number of atoms for each species. If None, set to be identical to the backbone structure.
        "V" correspond to vacancy (Not placing any atom at the lattice site).
        Note that the following condition must be satisfied.
            --> np.sum(num_spec_dict.values()) == len(backbone)
        E.g.) {'Ge': 3, 'Te': 4, 'V': 1}
           --> Condition: the backbone structure must have 8 atoms, totally.

    fix_ind_dict : dict or list or None
        Dict of atomic indices in the backbone for each species whose sites will not be shuffled.
        It can be a list. The site will not be shuffled for the atoms with indices included in the list.
        Note that every fixed atom will also have positional deviation, unless 'pin_the_fixed' set to be 'True'.
        Set to "None" if all atoms' positions must be shuffled.
        E.g.) {'Te': [0, 2, 4, 6]}
            * I.e. The lattice sites at the positions, backbone.positions[[0,2,4,6]], will be filled by Te atoms, and others will be shuffled.

    pin_the_fixed : Boolean
        If true, atoms included in the fix_ind_dict will not be deviated from their lattice sites (even if the cutoff_radi and cutoff_frac are not zero).

    Note)
        If both cutoff_radi and cutoff_frac are provided simultaneously, occurs error.
        If both cutoff_radi and cutoff_frac are not provided, cutoff_radi is set to be zero by default. i.e. no constraint

        cutoff_radi : Float
            Cutoff radius for the minimum distance between two atoms.
            If there is a pair of atoms closer than the cutoff_radi, the position of the atom will be randomized again.

        cutoff_frac : Float
            Set cutoff_radi as {cutoff_frac *2 *(expected value of covalent radius of the material)}.
            i.e. When cutoff_frac == 1, all atoms are spaced with at least two times of expected covalent radius.

    Note)
        If both random_radi and random_frac are provided simultaneously, occurs error.
        If both random_radi and random_frac are not provided, random_radi is set to be zero by default. i.e. no deviation

        random_radi : float
            Maximum magnitude of the random deviation vector.
            Each atom will deviated from its lattice point by a random vector shorter than random_radi.

        random_frac : Float
            Set random_radi as {random_frac *0.5 *(position of RDF first peak)}
            RDF = radial distribution function or pair correlation function
            !!Caution!!
                When this option is switched on, RDF calculation will be carried out, which is a time consuming process for big systems.

    strain : List of three floats e.g. [0.,0.,5.]
        Values specify how much you magnify the provided backbone cell.
        Cell gets longer along lattice vectors.
        positions will stay scaled positions of backbone atoms.

    strain_ratio : List of three floats e.g. [1.,1.1,1.]
        Values specify how much you magnify the provided backbone cell.
        Cell gets longer along lattice vectors.
        positions will stay scaled positions of backbone atoms.

    vacuum : List of three floats e.g. [0.,0.,5.]
        Values specify how much you magnify the provided backbone cell with vacuum.
        Cell gets longer along lattice vectors.
        positions will stay absolute positions of backbone atoms.
        insert vacuum after strain (if provided)

    vacuum_ratio : List of three floats e.g. [1,1.1,1]
        Values specify how much you magnify the provided backbone cell with vacuum.
        Cell gets longer along lattice vectors.
        positions will stay absolute positions of backbone atoms.
        insert vacuum after strain (if provided)

    max_trial_sec : float
        Maximum number of trials to get a new atomic position.
        If fails for more than "max_trial_sec" seconds, start to find new atoms from scratch.

    """

    ##
    from ss_util import count2list, list2count
    backbone = backbone.copy()
    # Get spec_list and num_spec_dict
    if num_spec_dict is None:
        spec_list = backbone.get_chemical_symbols()
        num_spec_dict = list2count(spec_list)
        num_vacancy = 0

    # Get spec_list.
    else:
        if 'V' in num_spec_dict.keys():
            num_vacancy = num_spec_dict['V']
            del(num_spec_dict['V'])
        else:
            num_vacancy = 0

        spec_list = count2list(num_spec_dict)

    # Get num_fix_dict and num_shffl_spec_dict
    num_fix_dict = {}
    from copy import deepcopy
    num_shffl_spec_dict = deepcopy(num_spec_dict)
    from numpy import ndarray
    if isinstance(fix_ind_dict, list) or isinstance(fix_ind_dict, np.ndarray):
        fix_ind_arr = np.array(deepcopy(fix_ind_dict))
        fixed_atoms = backbone.copy()[fix_ind_arr]
        fixed_spec = np.array(fixed_atoms.get_chemical_symbols())
        fix_ind_dict = {}
        for spec in np.unique(fixed_spec):
            fix_ind_dict[spec] = list(fix_ind_arr[fixed_spec==spec])
    if isinstance(fix_ind_dict, dict): 
        fix_ind_dict = deepcopy(fix_ind_dict)
        for key, value in fix_ind_dict.items():
            num_fix_dict[key] = len(value)
            num_shffl_spec_dict[key] -= len(value)
            fix_ind_dict[key] = np.array(value, dtype=np.int32).tolist()
            if key == 'V' and num_vacancy == 0:
                raise ValueError('The fix_ind_dict can not have "V", if num_spec_dict do not have "V"')
        if 'V' not in fix_ind_dict.keys():
            fix_ind_dict['V'] = []
            num_fix_dict['V'] = 0
            num_shffl_spec_dict['V'] = 0
    elif fix_ind_dict is None:
        fix_ind_dict = {}
        for key in num_spec_dict.keys():
            fix_ind_dict[key] = []
            num_fix_dict[key] = 0
        fix_ind_dict['V'] = []
        num_fix_dict['V'] = 0
    else:
        raise ValueError('Unknown type of fix_ind_dict.')

    # Get shffl_spec_list.
    shffl_spec_list = count2list(num_shffl_spec_dict)

    # Covalent bond length expectation value
    from ss_util import covalent_expect
    coval_expect = covalent_expect(spec_list)
                
    ## Cell strain adjustment.
    if strain_ratio is not None and strain is not None:
        raise ValueError("strain_ratio & strain parameters provided simultaneously. \
            Just provide one.")
    if strain is not None:
        strain = np.array(strain)
        if strain.shape != (3,):
            raise ValueError("Somethings wrong with strain parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        strain_ratio = strain / norm + 1
    if strain_ratio is not None:
        strain_ratio = np.array(strain_ratio)
        if strain_ratio.shape != (3,):
            raise ValueError("Somethings wrong with strain_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms = True,
            )
    if strain_ratio is None and strain is None:
        strain_ratio = [1.,1.,1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms = True,
            )

    ## Vacuum layer adjustment.
    if vacuum_ratio is not None and vacuum is not None:
        raise ValueError("vacuum_ratio & vacuum parameters provided simultaneously. \
            Just provide one.")
    if vacuum is not None:
        vacuum = np.array(vacuum)
        if vacuum.shape != (3,):
            raise ValueError("Somethings wrong with vacuum parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        vacuum_ratio = vacuum / norm + 1
    if vacuum_ratio is not None:
        vacuum_ratio = np.array(vacuum_ratio)
        if vacuum_ratio.shape != (3,):
            raise ValueError("Somethings wrong with vacuum_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms = False,
            )
    if vacuum_ratio is None and vacuum is None:
        vacuum_ratio = [1.,1.,1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms = True,
            )

    ## Determine cutoff radius.
    if cutoff_radi is not None and cutoff_frac is not None:
        raise ValueError("cutoff_radi & cutoff_frac parameters provided simultaneously. \
            Just provide one.")
    if cutoff_radi is not None:
        cutoff_r = cutoff_radi
    elif cutoff_frac is not None:
        cutoff_r = coval_expect * 2 * cutoff_frac
    else:
        cutoff_r = 0.

    ## Get random adjust radius
    if random_frac is not None and random_radi is None:
        supercell = make_supercell(backbone,[[2,0,0],[0,2,0],[0,0,2]])
        from ase.optimize.precon.neighbors import estimate_nearest_neighbour_distance as rNN
        rdf_1st_peak = rNN(supercell)
        ran_radi = rdf_1st_peak / 2 * random_frac
        if log:
            print("")
            print("********* Please check carefully !!!! ***********".center(120))
            print("RDF 1st peak / 2 == {:.2f}".format(rdf_1st_peak/2).center(120))
            print("Positional deviation degree == {:.2f}".format(random_frac).center(120))
            print("==> Random deviation radius == {:.2f}".format(ran_radi).center(120))
            print("it is {:.2f} % of covalent-bond-length expectation value.".format(ran_radi / coval_expect * 100).center(120))
            print("")
            print("C.f. ) covalent-bond-length expectation value == {:.2f}".format(coval_expect).center(120))
            print("C.f. ) cutoff radius == {:.2f}".format(cutoff_r).center(120))
            print("C.f. ) cutoff radius / covalent bond expectation *2 == {:.2f} %".format(cutoff_r / coval_expect / 2 * 100).center(120))
            print("")
    elif random_radi is not None and random_frac is None:
        ran_radi = float(random_radi)
    else:
        raise ValueError('Check random_radi or random_frac parameters.')

    ## Main
    if num_vacancy != 0:
        # Choose vacancy indices.
        vacancy_ind = np.random.permutation(
            np.setdiff1d(
                range(len(backbone)),
                np.concatenate(list(fix_ind_dict.values())),
                True,
                ),
            )[:num_vacancy - num_fix_dict['V']]
        # Add fixed-vacancy indicies to the array.
        vacancy_ind = np.concatenate([vacancy_ind, fix_ind_dict['V']]).astype(int)

        # Remove vacancies from the backbone.
        vacancy_bool = np.array([True] *len(backbone))
        vacancy_bool[vacancy_ind] = False
        backbone = backbone[vacancy_bool]

        # Update fix_ind_dict.
        del(fix_ind_dict['V'])
        for key, value in fix_ind_dict.items():
            for i in range(len(value)):
                value[i] -= np.sum(vacancy_ind < value[i])
            fix_ind_dict[key] = value
    fix_ind_arr = np.concatenate(list(fix_ind_dict.values())).astype(np.int32)

    from time import time
    len_atoms = len(backbone)
    old_posi  = backbone.get_positions()
    cell      = backbone.get_cell()
    cell_inv  = np.linalg.inv(cell)

    new_posi = []
    while len(new_posi) < len_atoms:
        # Start time of this loop.
        time_i = time()
        while True:
            # Attach one more atom.
            new_posi.append(old_posi[len(new_posi)].copy())
            # Give positional deviation to the lastest atom.
            if not pin_the_fixed or len(new_posi)-1 not in fix_ind_arr:
                direc_vec = np.random.rand(3)-0.5
                direc_vec /= np.linalg.norm(direc_vec)
                new_posi[-1] += direc_vec * ran_radi
            # Get minimum distance from latest atom.
            rel_new_posi = np.array(new_posi) @ cell_inv
            if len(new_posi) != 1:
                min_dist = np.min(np.linalg.norm(((rel_new_posi[:-1] - rel_new_posi[-1] + np.array([0.5]*3)) % 1.0 - np.array([0.5]*3)) @ cell, axis=1))
            else:
                min_dist = cutoff_r + 1.
            # Get elapsed time of this loop.
            time_f = time()
            time_d = time_f - time_i
            # If the latest atom is properly positioned, break this loop.
            if min_dist > cutoff_r:
                if log:
                    if len(new_posi) % 100 == 0:
                        print("( {} th / {} ) new atom position found".format(len(new_posi), len_atoms))
                break
            # If the latest atom is too close to another atom, remove the latest one.
            elif time_d < max_trial_sec:
                new_posi.pop()
            # If failed for more than "max_trial_sec" seconds, restart from scratch.
            else:
                new_posi = []
                break
    new_atoms = backbone.copy()
    new_atoms.set_positions(new_posi, apply_constraint = False)
    # Shuffle positions
    shuffle_ind = np.setdiff1d(range(len(new_atoms)), np.concatenate(list(fix_ind_dict.values())), True)
    new_positions = new_atoms.get_positions().copy()
    new_positions[shuffle_ind] = np.random.permutation(new_positions[shuffle_ind])
    new_atoms.set_positions(new_positions, apply_constraint = False)
    
    # Correct chemical symbols
    new_species = np.array(['XX']*len(new_atoms))
    new_species[shuffle_ind] = shffl_spec_list
    if len(fix_ind_arr):
        new_species[fix_ind_arr] = count2list(num_fix_dict)

    ## Set the chemical symbols
    new_atoms.set_chemical_symbols(new_species)
    # Sort by chemical numbers
    new_atoms = new_atoms[np.argsort(new_atoms.get_atomic_numbers())]

    return new_atoms

class Logger(object):

    """
    **** this class is totally copied from AMP's ****
    **** ref) https://amp.readthedocs.io/en/latest/ ****
    Logger that can also deliver timing information.

    Parameters
    ----------
    file : str
        File object or path to the file to write to.  Or set to None for
        a logger that does nothing.
    """
    def __init__(self, file):
        if file is None:
            self.file = None
            return
        if isinstance(file, str):
            self.filename = file
            file = open(file, 'a')
        self.file = file
        self.tics = {}

    def tic(self, label=None):
        """Start a timer.

        Parameters
        ----------
        label : str
            Label for managing multiple timers.
        """
        import time
        if self.file is None:
            return
        if label:
            self.tics[label] = time.time()
        else:
            self._tic = time.time()

    def __call__(self, message, toc=None, tic=False, no_space=None):
        """Writes message to the log file.

        Parameters
        ---------
        message : str
            Message to be written.
        toc : bool or str
            If toc=True or toc=label, it will append timing information in
            minutes to the timer.
        tic : bool or str
            If tic=True or tic=label, will start the generic timer or a timer
            associated with label. Equivalent to self.tic(label).
        """
        import time
        if self.file is None:
            return
        dt = ''
        if toc:
            if toc is True:
                tic = self._tic
            else:
                tic = self.tics[toc]
            dt = (time.time() - tic) # ssrokyz start
            dt = ' {:.1f} sec.'.format(dt) # ssrokyz end
        if self.file.closed:
            self.file = open(self.filename, 'a')
        if no_space or (no_space == None and toc):
            self.file.write(nospace(message + dt) + '\n')
        else:
            self.file.write(message + dt + '\n')
        self.file.flush()
        if tic:
            if tic is True:
                self.tic()
            else:
                self.tic(label=tic)
