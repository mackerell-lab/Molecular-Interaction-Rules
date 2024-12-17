<h1 align="center">Molecular Interaction Rules</h1>

<p align="center">
<img width="784" alt="Screenshot 2024-05-29 at 10 14 28 PM" src="https://github.com/mackerell-lab/Non-Covalent-Molecular-Interaction-Rules/assets/11812946/880e237a-f9a3-43d5-bb75-c7aeb756f28a">
</p>

Welcome to the Non-Covalent Molecular Interaction Rule Database. Molecules are recorded in their internal coordinate system and quantum mechanically optimized with `mp2/aug-cc-pvdz` geometry with manual community edits as needed on visual inspection. Monomers and Dimers for NCIs can be formed readily to the user.  

![Downloads](https://pepy.tech/badge/molecular-interaction-rules)


<h2 align="center">Quickstart</h2>

#### Install

<h3 align="center">Pip Installation</h3>


```bash

pip install molecular-interaction-rules 

```

<h3 align="center">Local Installation</h3>

Clone the repository and run

```
python -m pip install -e .
```


#### Get Atom Names

```python

from molecular_interaction_rules import MoleculerDatabase

molecules = MoleculerDatabase()
benzene_atom_names = molecules.get_atom_names('benzene')

print(benzene_atom_names)

```

Output:

```
['RC1', 'H1']
```

#### Get Monomer Coordinates

```python

from molecular_interaction_rules import MoleculerDatabase

molecules = MoleculerDatabase()
benzene_monomer = molecules.get_monomer_coordinates('benzene', 'RC1')

print (benzene_monomer)

```

Output:

```
X11
C11  X11  1.3940
C12  C11  1.3774 X11   60.0000
C13  C12  1.3774 C11  120.0000 X11    0.0000
C14  C13  1.3774 C12  120.0000 C11    0.0000
C15  C14  1.3774 C13  120.0000 C12    0.0000
C16  C15  1.3774 C14  120.0000 C13    0.0000
H11  C11  1.0756 C12  120.0000 C13  180.0000
H12  C12  1.0756 C11  120.0000 C13  180.0000
H13  C13  1.0756 C12  120.0000 C11  180.0000
H14  C14  1.0756 C13  120.0000 C12  180.0000
H15  C15  1.0756 C14  120.0000 C13  180.0000
H16  C16  1.0756 C15  120.0000 C11  180.0000
0 1

```

#### Get Dimer Coordinates - Benzene T-Shape

```python

from molecular_interaction_rules import MoleculerDatabase

molecules = MoleculerDatabase()
benzene_dimer = molecules.form_dimer_coordinates('benzene', 'RC1', 'benzene', 'RC1')

print (benzene_dimer)

```

Output
```
X11
C11  X11  1.3940
C12  C11  1.3774 X11   60.0000
C13  C12  1.3774 C11  120.0000 X11    0.0000
C14  C13  1.3774 C12  120.0000 C11    0.0000
C15  C14  1.3774 C13  120.0000 C12    0.0000
C16  C15  1.3774 C14  120.0000 C13    0.0000
H11  C11  1.0756 C12  120.0000 C13  180.0000
H12  C12  1.0756 C11  120.0000 C13  180.0000
H13  C13  1.0756 C12  120.0000 C11  180.0000
H14  C14  1.0756 C13  120.0000 C12  180.0000
H15  C15  1.0756 C14  120.0000 C13  180.0000
H16  C16  1.0756 C15  120.0000 C11  180.0000
0 1
--
X21   X11  DISTANCE  C11   ANGLE  C12   DIHEDRAL
C21  X21  1.3940  X11   90.0000  C11  180.0000
C22  C21  1.3774 X21   60.0000  X11   90.0000
C23  C22  1.3774 C21  120.0000 X21    0.0000
C24  C23  1.3774 C22  120.0000 C21    0.0000
C25  C24  1.3774 C23  120.0000 C22    0.0000
C26  C25  1.3774 C24  120.0000 C23    0.0000
H21  C21  1.0756 C22  120.0000 C23  180.0000
H22  C22  1.0756 C21  120.0000 C23  180.0000
H23  C23  1.0756 C22  120.0000 C21  180.0000
H24  C24  1.0756 C23  120.0000 C22  180.0000
H25  C25  1.0756 C24  120.0000 C23  180.0000
H26  C26  1.0756 C25  120.0000 C21  180.0000
0 1
```

#### Get Dimer Coordinates - Methylacetamide - Acetone

```python

from molecular_interaction_rules import MoleculerDatabase

molecules = MoleculerDatabase()
methylacetamide_acetone_dimer = molecules.form_dimer_coordinates('methylacetamide', 'H1', 'acetone', 'O1')

print (benzene_dimer)

```

Output
```
```


<h2 align="center">Moleculer Database</h2>


| Functional Group Class | Molecules  |
|-|-|
| Aromatic      | Azulene, Benzene, Bipyrrole, Bromobenzene, Chlorobenzene, Cytosine, Fluorobenzene, Four Pyridinone, Furan, Imidazole, Imidazolium, Indole, Indolizine, Iodobenzene, Isoxazole, Methylene Oxide, Nitrobenzene, 1 Phenyl-4-Pyridinone, Phenol, Phenoxazine, Pyridine, Pyridinium, Pyrimidine, Pyrrolidine, Thiophene, 3-Aminopyridine, 2-H-Pyran, Uracil |    | Alcohols      | Methanol |  
| Alkanes       | Cyclobutane, Cyclohexane, Cyclopropane, Neopentane, Propane |  
| Alkenes       | Cyclohexene, Cyclopentene, Methoxyethene, 1,3-dibutene, Propene, 2-Pyrroline |  
| Alkynes       | Propyne |  
| Amides        | Acetamide, Amidinium, Azetidinone, DimethylFormamide, Methylacetamide, Prolineamide, 2-pyrrolidinone |  
| Amines        | Ammonia, Dimethylamine, Ethyl Ammonium, Hydrazine, Methylamine, Piperidine, (Z)-N-methylethanimine, Tetramethylammonium, Trimethylamine, Triethylammonium |  
| Carbonyls     | Acetaldehyde, Acetate, Acetic Acid, Acetone, Carbon Dioxide, Formaldehyde, Methylacetate, Urea |  
| Ethers        | Dimethyl ether, Epoxide, Oxetane, Tetrahydrofuran, Tetrahydropyran |  
| Imines        | Ethenamine |  
| Halogens      | Bromoethane, Chloroethane, Dibromoethane, Dichloroethane, Fluoroethane, Difluoroethane, Tribromoethane, Trichloroethane, Trifluoroethane |  
| Nitriles      | Acetonitrile |  
| Organophosphorus      | Methyl Phosphate, Dimethyl Phosphate |  
| Organosulfur      | Dimethyl sulfone, Dimethyl Sulfoxide, Dimethyl trithiocarbonate, Dimethyl Disulfide, Ethylsulfanyl Phosphonic Acid, Methanethiol, Methylthiolate |  

<h2 align="center">Contact</h2>

Lead Developer: Suliman Sharif
Co-Authors: Anmol Kumar, Alexander D. MacKerell Jr.

© Copyright 2024 – University of Maryland School of Pharmacy, Computer-Aided Drug Design Center All Rights Reserved

