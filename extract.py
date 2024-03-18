import sys
sys.dont_write_bytecode = True
from ase.phasediagram import PhaseDiagram
from utilities import combinations as comb
from utilities import connect_higher_order as cho
from utilities import make_check_dir as mcd
from utilities import check_criterias
from utilities import entropy_summation
from ase.db import connect
from ase.io import write
from ase.formula import Formula
from pathlib import Path

import pandas as panda
import shutil
import csv 

class data_handler:
    def __init__(self, db_name=None, material_system=None, composition=None, exclude=None, additional_db=None, temperature=0):
         #db_name path to ase.db. if notthing els use MP.db
         #material_syste string denoting which system ('Mo-Al-B') 
         #composition dictionary of the material composition {'Mo': 3, 'Al': 1, 'B': 4}
         #exclude dictionary with exclude criterias, such as
           #composition, bool, if True exclude all of same composition
           #sg, string, if added exclude all of the following sg
           #if both active: exclude all compositions with the given sg
           #distance_from_hull: limits the phases with distance from hull 
         #additional_db used wehn attaching external db

         self.composition = composition
         self.components = comb(material_system.split('-'))[1:]
         self.components = cho(self.components)
         self.temperature = temperature
         self.phases = {}
         if not db_name: 
             db = connect("/proj/materialsdesigndivision/MP_database/MP.db")
         else:
             db = connect(db_name)
         for component in self.components:
             self.phases[component[0]] = []
             for row in db.select('elements='+str(component[0])):
##
 
                 if row.formula == 'FeW2B2':
                     print(row)
                     print(row.SG, row.get('SG'))
#
                 if exclude: 
                     if check_criterias(row, composition, exclude) == True: 
                        continue
                 self.phases[component[0]].append(row)
            
             if additional_db:
                 adb = connect(additional_db)
                 for adb_row in adb.select('elements='+str(component[0])):
                     if exclude: 
                         if check_criterias(adb_row, composition, exclude) == True: 
                            continue
                     self.phases[component[0]].append(adb_row)
    
    def return_all_phases_dict(self):
        return self.phases

    def create_csv(self, data=None, name=None):
        if not name:
            name='phases.csv'
        with open(name, 'w', encoding='UTF8') as f:
            writer = csv.writer(f)
            writer.writerow(['phases', 'energy'])
            for key in self.phases:
                for row in self.phases[key]:
                    if row.get('SQS') == True:
                        text = [row.get('formula'), row.get('DFT')+row.get('natoms')*self.temperature*row.get('dS')]
                    else: 
                        text = [row.get('formula'), row.get('DFT')]
                    writer.writerow(text)
        f.close()

    def calculate_ecp(self, competing_phases=None, temperature=None, mixing=None, coeffs_check=True):
        nan=False
        for row in self.components: 
            if "N-Na" in row: 
                nan=True
        df = panda.read_csv(competing_phases)
        if nan==True:
            inds = list(df.loc[panda.isna(df["phases"]), :].index)
            #inds = list(df.loc[panda.isnull(df["phases"]), :].index)
            df.drop(index=inds, axis=0, inplace=True)
        tuples = [tuple(x) for x in df.values]
        if nan==True: 
            tuples.append(('NaN', -7.100157))
        pd = PhaseDiagram(tuples)
        phase=''
        for key in self.composition: 
            phase = phase+key+str(self.composition[key])
        
        resulted_phases = "" 
        energy, indices, coefs = pd.decompose(formula=phase)
        coefs_tmp = []

        for idx, row in enumerate(indices):
            if coeffs_check == True:
                if round(coefs[idx], 3) <= 0:
                    continue
                coefs_tmp.append(coefs[idx])
                phases = str(Formula(str(df.iloc[row]['phases'])).reduce()[0])
                resulted_phases += phases + ' + '
                    
            else:  
                phases = df.iloc[row]['phases']
                resulted_phases += phases + ' + '
                coefs_tmp = coefs
        if resulted_phases:
            phases = resulted_phases[:-3]

 
        if temperature and mixing: 
           mix_elements = mixing.split('-')[0].split(':')
           mix_ratios = [ int(i) for i  in mixing.split('-')[1].split(':')]
            
           tot_mixing = sum(mix_ratios)/sum(self.composition.values())
           energy = energy + tot_mixing*entropy_summation(temperature=temperature, mixing=mix_ratios)
           

        return energy, indices, coefs_tmp, phases


    def create_competing_poscars(self, limit=None): 
        lst = list(self.phases.values())
        for row in lst:
            if not row: 
                continue 
            else: 
                for atomsrow in row:
                    idn = atomsrow.get('MPid')
                    write(format='vasp', filename='POSCAR_'+idn, images=atomsrow.toatoms(), direct=True)      
        mcd('POSCARS')
        for POSCAR_file in Path('.').glob('POSCAR_*'):
            shutil.move(str(POSCAR_file), 'POSCARS/')


    def create_poscar_from_id(self, mpid=None):
        mcd('POSCARS')
        for row in db.select('MPid='+str(mpid)):
            write(format='vasp', filename='POSCAR_'+mpid, images=atomsrow.toatoms())      
            shutil.move(str('POSCAR_'+mpid), 'POSCARS/')
