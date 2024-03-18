import pandas as panda
import numpy as np
import os
from ase.db import connect
from ase.io import write, read
import shutil
from os import path
import re
import math
#Utility functions

#MP Potentials
#Ta and Os are VASP recomended (More may be missing)
MP_pot=[("Ac","Ac"),("Ag","Ag"),("Al","Al"),("Ar","Ar"),("As","As"),("Au","Au"),("B","B"),("Ba","Ba_sv"),("Be","Be_sv"),("Bi","Bi"),("Br","Br"),("C","C"),("Ca","Ca_sv"),("Cd","Cd"),("Ce","Ce"),("Cl","Cl"),("Co","Co"),("Cr","Cr_pv"),("Cs","Cs_sv"),("Cu","Cu_pv"),("Dy","Dy_3"),("Er","Er_3"),("Eu","Eu"),("F","F"),("Fe","Fe_pv"),("Ga","Ga_d"),("Gd","Gd"),("Ge","Ge_d"),("H","H"),("He","He"),("Hf","Hf_pv"),("Hg","Hg"),("Ho","Ho_3"),("I","I"),("In","In_d"),("Ir","Ir"),("K","K_sv"),("Kr","Kr"),("La","La"),("Li","Li_sv"),("Lu","Lu_3"),("Mg","Mg_pv"),("Mn","Mn_pv"),("Mo","Mo_pv"),("N","N"),("Na","Na_pv"),("Nb","Nb_pv"),("Nd","Nd_3"),("Ne","Ne"),("Ni","Ni_pv"),("Np","Np"),("O","O"),("Os","Os_pv"),("P","P"),("Pa","Pa"),("Pb","Pb_d"),("Pd","Pd"),("Pm","Pm_3"),("Pr","Pr_3"),("Pt","Pt"),("Pu","Pu"),("Rb","Rb_sv"),("Re","Re_pv"),("Rh","Rh_pv"),("Ru","Ru_pv"),("S","S"),("Sb","Sb"),("Sc","Sc_sv"),("Se","Se"),("Si","Si"),("Sm","Sm_3"),("Sn","Sn_d"),("Sr","Sr_sv"),("Ta","Ta_pv"),("Tb","Tb_3"),("Tc","Tc_pv"),("Te","Te"),("Th","Th"),("Ti","Ti_pv"),("Tl","Tl_d"),("Tm","Tm_3"),("U","U"),("V","V_pv"),("W","W_pv"),("Xe","Xe"),("Y","Y_sv"),("Yb","Yb_2"),("Zn","Zn"),("Zr","Zr_sv")]


def identify_MP_pot(ele_string=None, index=0): 
    elements, potentials = zip(*MP_pot)
    for idx, row in enumerate(elements): 
        if row == ele_string: 
            index = idx
        else: 
            continue
    return (elements[index], potentials[index])


def SG_number_checker(number): 
    #return sg 
    sg = ["P1","P-1","P2","P21","C2","Pm","Pc","Cm","Cc","P2/m","P21/m","C2/m","P2/c","P21/c","C2/c","P222","P2221","P21212","P212121","C2221","C222","F222","I222","I212121","Pmm2","Pmc21","Pcc2","Pma2","Pca21","Pnc2","Pmn21","Pba2","Pna21","Pnn2","Cmm2","Cmc21","Ccc2","Amm2","Aem2","Ama2","Aea2","Fmm2","Fdd2","Imm2","Iba2","Ima2","Pmmm","Pnnn","Pccm","Pban","Pmma","Pnna","Pmna","Pcca","Pbam","Pccn","Pbcm","Pnnm","Pmmn","Pbcn","Pbca","Pnma","Cmcm","Cmce","Cmmm","Cccm","Cmme","Ccce","Fmmm","Fddd","Immm","Ibam","Ibca","Imma","P4","P41","P42","P43","I4","I41","P-4","I-4","P4/m","P42/m","P4/n","P42/n","I4/m","I41/a","P422","P4212","P4122","P41212","P4222","P42212","P4322","P43212","I422","I4122","P4mm","P4bm","P42cm","P42nm","P4cc","P4nc","P42mc","P42bc","I4mm","I4cm","I41md","I41cd","P-42m","P-42c","P-421m","P-421c","P-4m2","P-4c2","P-4b2","P-4n2","I-4m2","I-4c2","I-42m","I-42d","P4/mmm","P4/mcc","P4/nbm","P4/nnc","P4/mbm","P4/mnc","P4/nmm","P4/ncc","P42/mmc","P42/mcm","P42/nbc","P42/nnm","P42/mbc","P42/mnm","P42/nmc","P42/ncm","I4/mmm","I4/mcm","I41/amd","I41/acd","P3","P31","P32","R3","P-3","R-3","P312","P321","P3112","P3121","P3212","P3221","R32","P3m1","P31m","P3c1","P31c","R3m","R3c","P-31m","P-31c","P-3m1","P-3c1","R-3m","R-3c","P6","P61","P65","P62","P64","P63","P-6","P6/m","P63/m","P622","P6122","P6522","P6222","P6422","P6322","P6mm","P6cc","P63cm","P63mc","P-6m2","P-6c2","P-62m","P-62c","P6/mmm","P6/mcc","P63/mcm","P63/mmc","P23","F23","I23","P213","I213","Pm-3","Pn-3","Fm-3","Fd-3","Im-3","Pa-3","Ia-3","P432","P4232","F432","F4132","I432","P4332","P4132","I4132","P-43m","F-43m","I-43m","P-43n","F-43c","I-43d","Pm-3m","Pn-3n","Pm-3n","Pn-3m","Fm-3m","Fm-3c","Fd-3m","Fd-3c","Im-3m","Ia-3d"]
    return sg[int(number)-1]

def return_period_number(element): 
    period_PT = {
     '1': ["H", "He"],
     '2': ["Li", "Be", "B", "C", "N", "O", "F", "Ne"],
     '3': ["Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"],
     '4': ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
           "Ge", "As", "Se", "Br", "Kr"],
     '5': ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
           "Sn", "Sb", "Te", "I", "Xe"],
     '6': ["Cs", "Ba", "La", "Cr", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
           "Er", "Tm", "Yb", "Lb", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
           "Tl", "Pb", "Bi", "Po", "At", "Rn"],
     '7': ["Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
           "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
            "Nh", "Fl", "Mc", "Lv", "Ts", "Os"],
                } 
    for key in period_PT: 
        if element in period_PT[key]: 
            return key
# Needs to verify all the electronegativity values
# Based on Pauling : The chemical bond, 1967
def return_electronegativity(element_symbol):
    electronegativities = {
        "H": 2.20,  "He": None, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98, "Ne": None,
        "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": None, "K": 0.82, "Ca": 1.00,
        "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
        "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96, "Kr": None, "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33,
        "Nb": 1.60, "Mo": 2.16, "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69, "In": 1.78, "Sn": 1.96,
        "Sb": 2.05, "Te": 2.10, "I": 2.66, "Xe": None, "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14,
        "Pm": 1.13, "Sm": 1.17, "Eu": 1.20, "Gd": 1.20, "Tb": 1.20, "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25, "Yb": 1.10,
        "Lu": 1.27, "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.28, "Pt": 2.20, "Au": 2.54, "Hg": 2.00,
        "Tl": 1.62, "Pb": 2.33, "Bi": 2.02, "Po": 2.00, "At": 2.20, "Rn": None, "Fr": 0.70, "Ra": 0.90, "Ac": None, "Th": 1.30,
        "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.13, "Cm": 1.28, "Bk": 1.30, "Cf": 1.30, "Es": 1.30, "Fm": 1.30,
        "Md": 1.30, "No": 1.30, "Lr": None, "Rf": None, "Db": None, "Sg": None, "Bh": None, "Hs": None, "Mt": None, "Ds": None,
        "Rg": None, "Cn": None, "Nh": None, "Fl": None, "Mc": None, "Lv": None, "Ts": None, "Og": None,
    }

    element_symbol = element_symbol.capitalize()

    if element_symbol in electronegativities:
        electronegativity = electronegativities[element_symbol]
        if electronegativity is not None:
            return electronegativity
        else:
            return "Electronegativity not available for this element."
    else:
        return "Element symbol not found."

def return_atomic_number(element): 
    atomic_numbers = {
        "H": 1,  "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
        "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
        "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
        "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
        "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
        "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
        "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
        "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
        "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
        "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
        "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
    }

    element_symbol = element.capitalize()
    if element_symbol in atomic_numbers:
        return atomic_numbers[element_symbol]
    else:
        return None

def return_group_number(element): 
    group_PT = {
     '1': ["H", "Li", "Na", "K", "Rb", "Cs", "Fr"],
     '2': ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"],
     '3': ["Sc", "Y", "La", "Cr", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
           "Er", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
           "Cf", "Es", "Fm", "Md", "No", "Lr"],
     '4': ["Ti", "Zr", "Hf", "Rf"],
     '5': ["V", "Nb", "Ta", "Db"],
     '6': ["Cr", "Mo", "W", "Sg"],
     '7': ["Mn", "Tc", "Re", "Bh"],
     '8': ["Fe", "Ru", "Os", "Hs"],
     '9': ["Co", "Rh", "Ir", "Mt"],
     '10': ["Ni", "Pd", "Pt", "Ds"],
     '11': ["Cu", "Ag", "Au", "Rg"],
     '12': ["Zn", "Cd", "Hg", "Cn"],
     '13': ["B", "Al", "Ga", "In", "Tl", "Nh"],
     '14': ["C", "Si", "Ge", "Sn", "Pb", "Fl"],
     '15': ["N", "P", "As", "Sb", "Bi", "Mc"],
     '16': ["O", "S", "Se", "Te", "Po", "Lv"],
     '17': ["F", "Cl", "Br", "I", "At", "Ts"],
     '18': ["He", "Ne", "Ar", "Kr", "Xe", "Rn", "Os"],
            } 
    for key in group_PT: 
        if element in group_PT[key]: 
            return key

def return_valence(element): 
   
    valence_electrons = {
    '1': ["H", "Li", "Na", "K", "Rb", "Cs", "Fr"],
    '2': ["Be", "Mg", "Ca", "Sr", "Ba", "Ra", "Es", "Fm", "Md", "No"],
    '3': ["Sc", "Y", "La", "Ac", "B", "Al", "Ga", "In", "Tl", "Lw"],
    '4': ["Ti", "Zr", "Hf", "C", "Si", "Ge", "Sn", "Pb", "Ce", "Th"],
    '5': ["V", "Nb", "Ta", "N", "P", "As", "Sb", "Bi", "Pr", "Pa"],
    '6': ["Cr", "Mo", "W", "O", "S", "Se", "Te", "Po", "Nd", "U"],
    '7': ["Mn", "Tc", "Re", "F", "Cl", "Br", "I", "At", "Pm", "Np"],
    '8': ["Fe", "Ru", "Os", "Ne", "Ar", "Kr", "Xe", "Rn", "Og", "Sm", "Pu"],
    '9':  ["Co", "Rh", "Ir", "Eu", "Am"],
    '10': ["Ni", "Pd", "Pt", "Gd", "Cm"],
    '11': ["Cu", "Ag", "Au", "Tb", "Bk"],
    '12': ["Zn", "Cd", "Hg", "Dy", "Cf"],
    '13': ["Ho"],
    '14': ["Er"],
    '15': ["Tm"],
    '16': ["Yb", "Lu"],
                        }

    for key in valence_electrons: 
        if element in valence_electrons[key]: 
            return key


def entropy_summation(mixing=None):
    B_const = 0.00008617
    n_tot = sum(mixing)
    #print(n_tot)
    #print(mixing[0])
    #print(mixing[1])
    #print(mixing[0]/n_tot*np.log(mixing[0]/n_tot))
    #print(mixing[1]/n_tot*np.log(mixing[1]/n_tot))
        
    return B_const*(mixing[0]/n_tot*np.log(mixing[0]/n_tot) + mixing[1]/n_tot*np.log(mixing[1]/n_tot)) 
 

def combinations(system=None):
#extract all the combinations
#from the material system
    if len(system)==0:
        return [[]]
    combos = []
    for row in combinations(system[1:]):
        combos +=  [row, row+[system[0]]]
    return combos

def connect_higher_order(lst):
#fits the components to the db notation
#sorts and adds '-' for higher order systems
    tmp = []
    for idx, row, in enumerate(lst):
        if len(row) == 1:
            tmp.append(row)
        elif len(row) == 2:
            row.sort()
            tmp.append([row[0]+'-'+row[1]])
        elif len(row) == 3:
            row.sort()
            tmp.append([row[0]+'-'+row[1]+'-'+row[2]])
        elif len(row) == 4:
            row.sort()
            tmp.append([row[0]+'-'+row[1]+'-'+row[2]+'-'+row[3]])
    return tmp

def check_gcd(lst1, lst2): #Check the global common denominator
    gcd1 = np.gcd.reduce(list(lst1.values()))
    gcd2 = np.gcd.reduce(list(lst2.values()))
    lst1 = {key: value / gcd1 for key, value in lst1.items()}
    lst2 = {key: value / gcd2 for key, value in lst2.items()}
    if lst1 == lst2:
        return True
    else:
        return False

def check_criterias(row, composition, exclude): 
    rowcomp = row.count_atoms()
    sg = row.get('SG')
    if 'temperature' in exclude: 
        T = row.get('temperature')
        if T == None:
            pass
        else: 
            if T == exclude['temperature']: 
                pass
            else:
            #    print('TEMP: ', sg, row.get('elements'), row.get('formula'), row.get('DFT'), T, exclude['temperature'], composition, exclude)
                return True
            
#    print('passed temperature')
    if exclude['composition'] == True and 'sg' in exclude: 
        if len(rowcomp) == len(composition): 
            if check_gcd(rowcomp, composition)==True and sg == exclude['sg']:
#                print('SG+COMP: ', sg, row.get('elements'), row.get('formula'), row.get('DFT'), composition, exclude)
                return True

    elif exclude['composition'] == True and 'sg' not in exclude: 
        if len(rowcomp) == len(composition): 
            if check_gcd(rowcomp, composition): 
                return True

    if 'sg' in exclude and 'composition' not in exclude: 
        if sg == exclude['sg']: 
            return True

    if 'distance_from_hull' in exclude: 
        variable = row.get(str(list(exclude['distance_from_hull'].keys())[0]))
        if variable >  list(exclude['distance_from_hull'].values())[0]:
            return True 

def write_to_csv(filename, dct): 
    df = panda.DataFrame.from_dict(dct)
    df.to_csv(filename, index=False, header=True)


def make_check_dir(name): #make a new directory
    if not os.path.exists(name): 
        os.mkdir(name) 

def replace_string_in_file(filename, target_string, new_string):
    with open(filename) as f: 
        s = f.read()
        if target_string not in s: 
            print('"{target_string}" not found in {filename}.'.format(**locals()))
            return 
    with open(filename, 'w') as f:
        s = s.replace(target_string, new_string)
        f.write(s)

def create_poscar_from_id(mpid=None, name=None):
    db = connect('/proj/materialsdesigndivision/MP_database/MP.db')
    for row in db.select('MPid='+str(mpid)):
        write(format='vasp', filename='POSCAR_'+mpid, images=row.toatoms())
        if name:
            shutil.move(str('POSCAR_'+mpid), 'POSCARS/POSCAR_'+str(name))
        else: 
            shutil.move(str('POSCAR_'+mpid), 'POSCAR')


def read_and_write_potcar(pp_source=None, pp=None, path=None, append=None):
    with open(pp_source+'/'+pp+'/POTCAR', 'r+') as read_pp:
        lines_potcar = read_pp.readlines()
    read_pp.close()
    if append > 0:
        with open('POTCAR', 'a') as write_pp:
            for line in lines_potcar:
                write_pp.write(line)
    else:
        with open('POTCAR', 'w') as write_pp:
            for line in lines_potcar:
                write_pp.write(line)
    write_pp.close()

def check_path_and_energy(source):
    if not path.exists(source):
        return False, 'no energy'

    with open(source, 'r') as output:
        line = output.readline()
        line = line.split()

    if len(line) <=3:
        return False, 'check energy'

    return True, 'ok', line

def grep_delta_mag(fname, last_steps=None):
    E0_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if 'POSCAR found : ' in line:
                n = int(line.split()[6])
            if 'E0' in line:
                E0_list.append(line)
    if last_steps:
        E0_list = E0_list[-2:]
    e1 = float(E0_list[0].split()[9])
    e2 = float(E0_list[-1].split()[9])
    dE0=(e2-e1)/n*1000
    return dE0

def grep_delta_E0(fname, last_steps=None, mag=False):
    E0_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if 'POSCAR found : ' in line:
                n = int(line.split()[6])
            if 'E0' in line:
                E0_list.append(line)
    if last_steps:
        E0_list = E0_list[-2:]
    e1 = float(E0_list[0].split()[4])
    e2 = float(E0_list[-1].split()[4])
    dE0=(e2-e1)/n*1000
    return dE0

def grep_e0(fname, software='VASP'):
    E0 = ""
    if software=='VASP':
        tag = 'free  energy   TOTEN' 
    if software=='QE':
        tag = '!    total energy' 
    with open(fname, 'r') as fin:
        for line in fin:
            if tag in line:
                E0 = line
    return float(E0.split()[4])



def export_to_db(db_name='tmp.db', path_to_db=None, struct_name=None, outcar_name=None, symprec=1e-1, sg=None, sqs=None, focus_elements=None, count_atoms=False, count_elements=None):
    from ase.spacegroup import get_spacegroup
    db = connect(path_to_db+'/'+db_name)

    if not struct_name or not outcar_name:
        for fname in os.listdir('.'):
            if 'CONTCAR' in fname:
                struct_name = fname
            if 'OUTCAR' in fname:
                outcar_name = fname

    atoms = read(struct_name)

    if not sg:
        sg = get_spacegroup(atoms, symprec=symprec)

    e0 = grep_e0(outcar_name)
    n = atoms.get_global_number_of_atoms()
    elements = list(set(atoms.get_chemical_symbols()))
    elements.sort()
    kvp={
            'SG': sg,
            'elements': "-".join(elements),
            'DFT': e0,
            'DFT_per_atom': e0/n,
            'nsites': n,
#            'm1': focus_elements[0],
#            'm2': focus_elements[1],
            'dS': 0,
            }
    if focus_elements: 
        for idx, row in enumerate(focus_elements): 
            kvp['m'+str(idx+1)]=row

    if sqs:
        dct = {}
        for symbol in atoms.get_chemical_symbols(): 
            dct[symbol] = dct.get(symbol, 0) + 1 

        mixing_data = [], [] 
        for ele in dct: 
            if ele not in focus_elements: 
                continue 
            mixing_data[0].append(str(dct[ele]))
            mixing_data[1].append(ele)
        mixing = mixing_data[1][0]+':'+mixing_data[1][1]+' '+mixing_data[0][0]+':'+mixing_data[0][1]
        kvp['SG'] =  sg #P1
        kvp['SQS'] =  True
        kvp['mixing'] =  mixing
#        kvp['orgi_SG'] = sg 
        
        ds = sum([int(i) for i in mixing_data[0]])/sum(dct.values())*entropy_summation(mixing = [int(i) for i in mixing_data[0]])       
#        print(ds) 
        kvp['dS'] = ds
    if count_atoms == True:
   
        kvp['nM1'] = count_elements[0]
        kvp['nM2'] = count_elements[1]
        kvp['nBoron'] = count_elements[2]

    db.write(atoms=atoms, key_value_pairs=kvp) 

def combine_dbs(db1, db2, new_db=None):
    if new_db: 
        shutil.copy(db1, new_db)
        final_db = connect(new_db)
    else: 
        final_db = connect(db1)
    db2 = connect(db2)

    for row in db2.select(): 
        a = row.toatoms()
        kvp = row.get('key_value_pairs')
        final_db.write(atoms=a, key_value_pairs=kvp)
   


def create_latin(cellvalues, matrix, positions, elements, mixing, conc): 
    cellvalues = [ str(int(i)) for i in cellvalues ]
    cellvalues = " ".join(cellvalues)

    f = open('lat.in', 'a')
    
    f.write(cellvalues+' \n')
    for row in matrix: 
        f.write(row+' \n')
    
    for idx, row in enumerate(elements): 
        coords = [str(i) for i in positions[idx]]
        coords = " ".join(coords)
      
        if row in mixing: 
            f.write(coords+' '+mixing[0]+'='+str(conc)+','+mixing[1]+"="+str(1-conc)+' \n')
        else: 
            f.write(coords+' '+row+' \n')
    f.close()
def create_sqs_run(time='01:00:00', name='tmp', n=10): 
    f = open('RUN.sh', 'a')
    f.write("#!/bin/bash \n")
    f.write("#SBATCH -N 1 \n")
    f.write("#SBATCH -t "+time+'\n')
    f.write("#SBATCH --ntasks-per-node=1 \n")
    f.write("#SBATCH -J "+name+' \n')
    f.write("export PATH=/software/sse/manual/ATAT/3.34/nsc1-gcc-2018a-eb/bin:$PATH \n")
    f.write("corrdump -l=lat.in -ro -noe -nop -clus -2=10 -3=3.4 \n")
    f.write("mcsqs -l=lat.in -n="+n+" -2=10 -3=3.4 ; getclus \n")
    f.close()


def sort_elements(count, elements): 
    count, elements = zip(*sorted(zip(count, elements)))  
    return count, elements 

def count_atoms_obj(atoms):
        dct = {}
        for symbol in atoms.get_chemical_symbols(): 
            dct[symbol] = dct.get(symbol, 0) + 1
        return dct 

def append_or_equal(dct, kvp): 
    denoter = kvp['m1']+'-'+kvp['m2']
    if denoter not in dct: 
        dct[denoter] = kvp
    else: 
        if dct[denoter]['DFT_per_atom'] > kvp['DFT_per_atom']: 
            dct[denoter] = kvp
         
def extract_MP_pot(element): 
    MP_list = list(zip(*MP_pot))
    ele_index = MP_list[0].index(element)
    return MP_list[1][ele_index]

def generate_lobsterin(path=None, startenergy=None, endenergy=None, basisset=None, basisfunctions=None, bondlengths=None): 
    if path:  
        f = open(str(path)+'lobsterin', 'a')
    else: 
        f = open('lobsterin', 'a')
    f.write("COHPstartEnergy "+str(startenergy)+"\n")
    f.write("COHPendEnergy "+str(endenergy)+"\n")
    f.write("basisSet "+str(basisset)+"\n")
    for key in basisfunctions: 
        f.write("basisfunctions "+str(key)+" "+str(basisfunctions[key])+"\n")
    f.write("cohpGenerator from "+str(bondlengths[0])+" to "+str(bondlengths[1])+"\n")
    f.close()

def split_string_and_int(string): 
    match = re.match(r"([a-z]+)([0-9]+)", string, re.I)
    if match: 
        items = match.groups()
    return items


def split_and_extract_lines(lines, split_string, number_of_lines_after):
#    lines = lines.split('split_string')
    for idx, line in enumerate(lines): 
        if split_string in line: 
            split_index = idx
    lines = lines[split_index: split_index + number_of_lines_after]
    return lines

def check_magnetism(phase):
    magnetic = ['Mn', 'Cr', 'Fe', 'Co']
    for ele in magnetic:
        if ele in phase:
            return True
    return False

def number_of_cores_tetralith(n): 
    if n <= 32: 
        return 1
    elif n > 32 and n <= 64: 
        return 2
    elif n > 64 and n <= 96: 
        return 3
    elif n > 96 and n <= 128: 
        return 4
    elif n > 128 and n <= 160: 
        return 5
    elif n > 160 and n <= 192: 
        return 6
    elif n > 192 and n <= 224: 
        return 7
    elif n > 224 and n <= 256: 
        return 8
    elif n > 256 and n <= 288: 
        return 9
    else: 
        return 10

def bash_grep(pattern, file): 
    with open(file, 'r') as f: 
        lines = f.readlines()
    
    matched_lines = [line for line in lines if re.search(pattern, line)]
    return matched_lines
