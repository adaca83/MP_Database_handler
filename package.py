from utilities import split_string_and_int as ssai
import pandas as pd
import os
#USPEX handling

class read_USPEX_data():
    def __init__(self, files=[], file_path=None, limit=None, elements=None):
        self.data = {}
        for row in files:
            if not row in self.data:
                self.data[row] = {'energy': [], 'composition': [], 'gen': [], 'id': [], 'sg': [], 'phase': [], 'natoms': []}
            with open(file_path+'/'+row, 'r') as f:
                lines = f.read().splitlines()
                del lines[0] #tar bort första raden
                del lines[0] #tar bort andra raden
                for line in lines:
                    ' '.join(line.split()) #Tar bort överflöd av alla " "
                    line=line.split()
                    n=int(line[4])+int(line[5])+int(line[6])
                    phase = elements[0]+line[4]+elements[1]+line[5]+elements[2]+line[6]
                    self.data[row]['energy'].append(float(line[8])/n)
                    self.data[row]['gen'].append(line[0])
                    self.data[row]['id'].append(line[1])
                    self.data[row]['sg'].append(line[16])
                    self.data[row]['phase'].append(phase)
                    self.data[row]['natoms'].append(n)
                    self.data[row]['composition'].append([int(line[4]), int(line[5]), int(line[6])])

    def return_USPEX_data(self):
        return self.data


class COHP_analyser(): 
    def __init__(self, path=None, POSCAR_name=None):
        self.phase = {}
        if path: 
            self.path = path
        else: 
          self.path = '.'
        files = os.listdir(self.path) 
        if 'lobsterout' not in files: 
            print('missing lobsterout')

        if not POSCAR_name: 
            POSCAR_name = 'POSCAR'
        with open(str(self.path)+'/'+str(POSCAR_name), 'r') as POSCAR_file:
     
            lines = POSCAR_file.readlines()
            elements = lines[5].split()
            print(elements)
            n_atoms = lines[6].split()
            for idx, n in enumerate(elements):
                if n not in self.phase:  
                    self.phase[n] = int(n_atoms[idx])
                else: 
                    self.phase[n] += int(n_atoms[idx])
                #self.dos_dct['interactions'][n] = []
        print(self.phase)



        
    def check_spilling(self):
        self.spillings = []
        with open(str(self.path)+'/lobsterout', 'r') as lobsterout:
            for line in lobsterout: 
                if 'spilling' in line:
                    self.spillings.append(line)
        return self.spillings           

    def ICOHP(self, limits=None, elements=None):
        self.ICOHP = {'cohp_nr': [], 'atomMU': [], 'atomNU': [], 'distance': [], 'translation': [], 'ICOHP_at_ef': []}
        with open(str(self.path)+'/ICOHPLIST.lobster', 'r') as ICOHP_file: 
            for idx, line in enumerate(ICOHP_file):
                if idx == 0: continue
                line_str=line.split() 
                self.ICOHP['cohp_nr'].append(int(line_str[0]))
                self.ICOHP['atomMU'].append(line_str[1])
                self.ICOHP['atomNU'].append(line_str[2])
                self.ICOHP['distance'].append(float(line_str[3]))
                self.ICOHP['translation'].append([int(line_str[4]), int(line_str[5]), int(line_str[6])])
                self.ICOHP['ICOHP_at_ef'].append(float(line_str[7]))
       
    def return_ICOHP(self):
        return self.ICOHP

    def pCOHP(self, magnetic=False, prec=3, limits={'energy': [-30, 20]} ):
        print(self.phase)
        self.cohp_dct = {
                         'energy': [],
                         'avg_COHP': [],
                         'total_COHP': [],
                         'interactions': {},
                        }
        with open(str(self.path)+'/COHPCAR.lobster', 'r') as COHP_file: 
            lines = COHP_file.readlines()
            self.n_int = int(lines[1].split()[0])-1

        for idx, idx_int in enumerate(lines[3:4+self.n_int-1]): 
            
            row = idx_int.split(':')[1].split('(')
            distance = row[1]
            elements = row[0].split('->')
            elements = [ssai(elements[0]), ssai(elements[1])]
            int_num = idx + 1

            interaction = str(elements[0][0])+'-'+str(elements[1][0])

            if interaction not in self.cohp_dct['interactions']: 
                self.cohp_dct['interactions'][interaction] = {
                                         'interaction_numbers' : [], 
                                         'pCOHP' : [], 
                                         }
            self.cohp_dct['interactions'][interaction]['interaction_numbers'].append(int_num)
        
        for idx, idx_int in enumerate(lines[int_num+4:]): 
            row = idx_int.split()

            energy = round(float(row[0]), prec)
            cohp_total = round(float(row[1]), prec)
            if 'energy' in limits: 
                if energy > limits['energy'][0] and energy < limits['energy'][1]: 
                    self.cohp_dct['energy'].append(energy)
                    self.cohp_dct['avg_COHP'].append(cohp_total)
               
            
                    pcohp_interactions = row[3::2]
                    total_cohp = 0
                    for interaction in self.cohp_dct['interactions']:
                         pcohp = sum([ float(pcohp_interactions[int(i)-1]) for i in self.cohp_dct['interactions'][interaction]['interaction_numbers']])
                         n = len(self.cohp_dct['interactions'][interaction]['interaction_numbers'])
                         self.cohp_dct['interactions'][interaction]['pCOHP'].append(pcohp/n)
                         total_cohp += pcohp/n
                    self.cohp_dct['total_COHP'].append(total_cohp) 
                        


    def COHP_avg(self, magnetic=False, prec=3, limits={'energy': [-30, 20]}):
        self.cohp_dct = {
                         'energy': [],
                         'avg_COHP': [],
                        }
        with open(str(self.path)+'/COHPCAR.lobster', 'r') as COHP_file: 
            lines = COHP_file.readlines()
            self.n_int = int(lines[1].split()[0])-1

        for idx, idx_int in enumerate(lines[3:4+self.n_int-1]): 
            
            row = idx_int.split(':')[1].split('(')
            distance = row[1]
            elements = row[0].split('->')
            elements = [ssai(elements[0]), ssai(elements[1])]
            int_num = idx + 1
                
        for idx, idx_int in enumerate(lines[int_num+3:]): 
            row = idx_int.split()

            energy = round(float(row[0]), prec)
            cohp = round(float(row[1]), prec)
            if 'energy' in limits: 
                if energy > limits['energy'][0] and energy < limits['energy'][1]: 
                    self.cohp_dct['energy'].append(energy)
                    self.cohp_dct['avg_COHP'].append(cohp)
  
    def return_COHP_avg(self): 
        return self.cohp_dct


    def DOS_avg(self, magnetic=False, prec=3, limits={'energy': [-30, 20]}):
        self.dos_dct = {
                         'energy': [],
                         'avg_DOS': [],
                        }

        with open(str(self.path)+'/DOSCAR.lobster', 'r') as DOS_file: 
            lines = DOS_file.readlines()
            self.dos_rows = int(lines[5].split()[2])
        
        for idx, idx_int in enumerate(lines[6:self.dos_rows+6]):
  
            row = idx_int.split()
            self.dos_dct['energy'].append(round(float(row[0]), prec))
            self.dos_dct['avg_DOS'].append(round(float(row[1]), prec))
 

#    def return_DOS_avg(self): 
#        return self.dos_dct

    def pDOS(self, magnetic=False, prec=3, limits={'energy': [-30, 20]}):
        self.dos_dct = {
                         'energy': [],
                         'avg_DOS': [],
                         'total_DOS': [],
                         'interactions': {},
                        }
        print(self.phase)
        for key in self.phase:
            self.dos_dct['interactions'][key] = []             


#        with open(str(self.path)+'/'+str(POSCAR_name), 'r') as POSCAR_file: 
#            lines = POSCAR_file.readlines()
#            elements = lines[5].split()
#            n_atoms = lines[6].split()
#            for idx, n in enumerate(elements):
#                if n not in self.phase:  
#                    self.phase[n] = int(n_atoms[idx])
#                else: 
#                    self.phase[n] += int(n_atoms[idx])
#                self.dos_dct['interactions'][n] = []
#        print(self.phase) 


        with open(str(self.path)+'/DOSCAR.lobster', 'r') as DOS_file: 
#            sep_lines = DOS_file.readlines()
#            self.dos_rows = int(sep_lines[5].split()[2])

#            all_lines = DOS_file.read()
#            all_lines = all_lines.split('; Z =')
            
            sep_lines = DOS_file.readlines()
            self.dos_rows = int(sep_lines[5].split()[2])
            print(self.dos_rows)
#        print(self.phase)       
        for idx, idx_int in enumerate(sep_lines[6:self.dos_rows+6]):
            row = idx_int.split()
            self.dos_dct['energy'].append(round(float(row[0]), prec))
            self.dos_dct['avg_DOS'].append(round(float(row[1]), prec))

        counter = 1
        for item in self.phase.items(): 
            n = item[1]
            element = item[0]
            for idx, line in enumerate(sep_lines[self.dos_rows+7: 2*self.dos_rows]):
               pdos = 0

               line = line.split()
               energy = float(line[0])
               if 'energy' in limits: 
                   if energy < limits['energy'][0] or energy > limits['energy'][1]:
                   
                       continue
#               if '-2.00200' in line: 
#                   print(line, counter, n+1)
               for i in range(counter, n+counter):
#                   if i == 1: 
#                   else: 
              #     print(sep_lines[idx + i*(self.dos_rows) + i-1].split(), counter, n+1, i)

                   if '; ' in sep_lines[idx + i*self.dos_rows] or "; " in sep_lines[idx + i*self.dos_rows + i-1]:
                       continue
                   else:
                       #if '-2.00200' in line:
                       if i == 1:
                           pdos += sum( [float(j) for j in sep_lines[idx + 7 + i*self.dos_rows].split()][1:])
#                           print(sep_lines[idx + 7 + i*self.dos_rows].split(), pdos) 
                       else:
#                           print(sep_lines[idx +7 + i*self.dos_rows + i-1 ].split()) 
                           pdos += sum( [float(j) for j in sep_lines[idx+7 + i*self.dos_rows + i-1 ].split()[1:]])
#                           print(sep_lines[idx +7 + i*self.dos_rows + i-1 ].split(), pdos) 
#               if '-2.0' in line: 
#                   print(pdos/n) 
#          rit(pdos/n)   
               self.dos_dct['interactions'][element].append(pdos)
            counter += n
                       
                   
#            n = self.phase[ele]
#            for 

#        print(all_lines)
#        print(len(all_lines))

        

#        total_value = 7
#        for row in enumerate(self.phase):
#            pdos = 0
#            for n in range(total_value, total_value+self.phase[row[1]]):
#                if n == 0: 
#                    index_list = [1009]
#                else: 
#                    index_list.append(n*1000 + 1010 + n*2-1)
#   
#            for idx, i in enumerate(sep_lines[n]): 
#                    pdos += sum( [ float(j) for j in i[1:]] ) 
#            if self.dos_dct['interactions'][row[1]] == None:
#                self.dos_dct['interactions'][row[1]].append(pdos)
#            else: 
#                self.dos_dct['interactions'][row[1]][idx] += pdos
#
#            total_value  += self.phase[row[1]]
#               

    def return_DOS_avg(self): 
        return self.dos_dct

 
