# MP_database_handler Manual

This is a manual for the package written by [Adam Carlsson](https://liu.se/medarbetare/adaca83) which goal is to ease the pathway of thermodynamical phase stability predictions by utilizing the [Materials Project database](https://next-gen.materialsproject.org/) and the package [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/). A brief introduction to the concepts of thermodynamical phase stability predictions is given below

## Thermodynamical phase stability prediction

The thermodynamical phase stability of an arbitrary material may be evaluated by correlating the energy with a set of competing phases. Different thermodynamical phase stability descriptors are often destinguished by the included materials in the set of competing phases. For example, the **formation energy** descriptor compares the energy of an arbitrary materials with a linear combination of the elements, in bulk form, which constitues the material in focus, the **iso-structural formation enthalpy** compares the energy of an arbitrary material with a set of competing phases all with the same strucrue, and the **formation enthalpy** which compares the energy of an arbitrary material with a set of most competing phases of the material system. The set of most competing phases is basically the linear combination of compounds within the studied materials system that produces the lowest energy and forms the same composition as the material in focus. 

Below follows brief definitions of the different descriptors for an arbitrary material system A<sub>2</sub>BC<sub>2</sub> with elements A, B, and C

**Formation energy**<br>
compares the energy with a linear combination of the signle elements in bulk form which constitues the material in focus.<br><br>
<div align="center">E<sub>f</sub> = E(A<sub>2</sub>BC<sub>2</sub>) - (2A + B + 2C)<div align="left"><br>

**Iso-structural formation enthalpy**<br>
compares the energy with a set of competing phases whith the same structures. Commonly used to describe mixing tendencies in material systems such as A<sub>1-x</sub>B<sub>x</sub>C<sub>2</sub>. The iso-structural formation enthalpy in this case compares the energy of A<sub>1-x</sub>B<sub>x</sub>C<sub>2</sub> with its endstructures (x = 0, 1).<br><br>
<div align="center">H<sub>iso</sub> = E(A<sub>1-x</sub>B<sub>x</sub>C<sub>2</sub>) - (xE(AC<sub>2</sub>) + (1-x)E(BC<sub>2</sub>))<div align="left"><br>

**Formation enthalpy**<br>
compares the energy op an arbitrary material, A<sub>2</sub>BC<sub>2</sub>, with the linear combination of competing phases within the A-B-C material system which yields the lowest energy and forms the same composition as the material in focus. <br><br>
<div align="center">H_{cp} = E(A<sub>2</sub>BC<sub>2</sub>) - minE<sub>cp</sub>(2<sup>A</sup>, 1<sup>B</sup></sup>, 2<sup>C</sup>)<div align="left"><br>

# extract.py

```data_handler``` is a python class for evaluating the **formation enthalpy** using the Materials Project database, in addition to alternative external ASE databases, as set of competing phases. It utilizes ```ase.phasediagram``` to derive the decomposition energy and phases. In addition to correlating the energy to the materials stored in the Materials Project database, the script alowes users to also correlate them to external ASE databases including both chemically ordered and disordered materials with an entropic contribution which may be evaluated at elevated temperatures. <br><br>

The ```data_handler``` class may take the following inputs <br>
```python
data_handler( db_name = None,          # The name of the database
              material_system = None,  # The specific material system in focus, i.e., A₂BC₂ would have A-B-C
              composition = None,      # The composition in focus in a dictionary, i.e., {A: 2, B: 1, C: 2}
              exclude = None,          # Criterias what and what not to include in the set of competing phases
              additional_db = None,    # If an external database should be included in addition to Materials Project
              temperature = 0,         # If the entropy contribution should be evaluated at elevated temperatures
            )
```
<br><br>
**Example** <br>
Lets evaluate the phase stability for Mo<sub>2</sub>AlB<sub>2</sub> with only Materials Project as the set of competing phases. For this, we need to initiate the ```data_handler``` class and use ```data_handler.create_csv``` followed by the ```data_handler.calculate_ecp``` function. 
```python
from extract import data_handler     # first we import the class data_handler from the package estract.py

system = "Mo-Al-B"                    # define the material system
phase = {"Mo": 2, "Al": 1, "B": 2}    # define the material as a dictionary

dbh = data_handler( material_system = system                            # initiate the data_handler class
                    composition = phase                                 # by only giving it the specific system
                  )                                                     # and composition

dbh.create_csv() # create a csv file with the competing phases

ecp_energy, indices, coefs, ecp_phases = dbh.calculate_ecp(competing_phases="phases.csv")     # evaluate the energy of the set of most competing phases

```
```data_handler.create_csv``` in this case create a csv file with the set of competing phases titles "phases.csv". The content of the file looks something accordingly: 
```
phases,energy
Mo,-10.41927976
Mo,-9.97480741
Mo4,-41.61408967
Al4,-13.77135701
Al4,-14.93651943
Al,-3.74557583
Al50Mo10,-310.18030452
Al2Mo6,-74.49671068
Al12Mo4,-93.46687211
Al10Mo2,-62.00682951
Al16Mo24,-317.62690291
B28,-186.25465764
B50,-322.78944777
B,-5.9718089
Mo7B24,-240.57687211
Mo3B12,-112.70558374
Mo4B2,-58.64037758
Mo2B6,-64.27747392
MoB2,-25.03127307
Al23B50,-423.83275243
AlB2,-17.2333181
Al2Mo2B2,-45.24105473
.
.
.
```

This file is essentially the input requiered when evaluating the phase stability as it contains the complete set of competing phases. The more competing phases included means the better representation of the complete phase space. The linear combination of these phases which produces the lowest energy is made through the 

```python
ecp_energy, indices, coefs, ecp_phases = dbh.calculate_ecp(competing_phases="phases.csv")
```

Here, ```ecp_energy``` is the energy of the linear combination of set of most competing phases§ for the given formula unit (in this case 5 atoms since we are studying Mo<sub>2</sub>AlB<sub>2</sub>). ```indices``` helps users keep track of which competing phase in phases.csv is included in the set of most competing phases. ```coefs``` is the coefficients, or the fraction, of how much each competing phase contributes. ```ecp_phases``` is a string of the set of most competing phases. 












