# MP_database_handler Manual

This is a manual for the package written by [Adam Carlsson](https://liu.se/medarbetare/adaca83) which goal is to ease the pathway of thermodynamical phase stability predictions by utilizing the [Materials Project database](https://next-gen.materialsproject.org/) and the package [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/). A brief introduction to the concepts of thermodynamical phase stability predictions is given below

## Thermodynamical phase stability prediction

The thermodynamical phase stability of an arbitrary material may be evaluated by correlating the energy with a set of competing phases. Different thermodynamical phase stability descriptors are often destinguished by the included materials in the set of competing phases. For example, the **formation energy** descriptor compares the energy of an arbitrary materials with a linear combination of the elements, in bulk form, which constitues the material in focus, the **iso-structural formation enthalpy** compares the energy of an arbitrary material with a set of competing phases all with the same strucrue, and the **formation enthalpy** which compares the energy of an arbitrary material with a set of most competing phases of the material system. The set of most competing phases is basically the linear combination of compounds within the studied materials system that produces the lowest energy and forms the same composition as the material in focus. 

Below follows brief definitions of the different descriptors for an arbitrary material system A$_{2}$BC$_{2}$ with elements A, B, and C

**Formation energy**
compares the energy with a linear combination of the signle elements in bulk form which constitues the material in focus
$E_{f}$ $=$ $E(A_{2}BC_{2})$ - ($2A$ + $B$ + $2C$)

**Iso-structural formation enthalpy**
compares the energy with a set of competing phases whith the same structures. Commonly used to describe mixing tendencies in material systems such as A$_{1-x}$B$_{x}$C$_{2}$. The iso-structural formation enthalpy in this case compares the energy of A$_{1-x}$B$_{x}$C$_{2}$ with the endstructures ($x$ = 0, 1)
$H_{iso}$ $=$ $E(A_{1-x}B_{x}C_{2})$ - ($x$E(AC$_{2}$) + (1-x)E(BC$_{2}$))

**Formation enthalpy**
compares the energy op an arbitrary material, A$_{2}$BC$_{2}$, with the linear combination of competing phases within the A-B-C material system which yields the lowest energy and forms the same composition as the material in focus. 
$H_{cp}$ $=$ $E(A_{2}B_§{x}C_{2})$ - minE$_{cp}$(2$^{A}$, 1$^{B}$, 2$^{C}$) 
