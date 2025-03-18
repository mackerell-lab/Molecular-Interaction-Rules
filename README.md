<h1 align="center">Molecular Interaction Rules</h1>

<p align="center">
<img width="784" alt="Screenshot 2024-05-29 at 10 14 28 PM" src="https://github.com/mackerell-lab/Non-Covalent-Molecular-Interaction-Rules/assets/11812946/880e237a-f9a3-43d5-bb75-c7aeb756f28a">
</p>

Welcome to the Non-Covalent Molecular Interaction Rule Database. Molecules are recorded in their internal coordinate system and quantum mechanically optimized with `mp2/aug-cc-pvdz` geometry with manual community edits as needed on visual inspection. Monomers and Dimers for NCIs can be formed readily to the user.  

![Downloads](https://pepy.tech/badge/molecular-interaction-rules)

#### Pip Installation


```bash

pip install molecular-interaction-rules 

```

#### Local Installation

Clone the repository and run

```
python -m pip install -e .
```

#### Import MoleculerDatabase

```python
from molecular_interaction_rules import MoleculerDatabase

molecules = MoleculerDatabase()
```

#### Get Full List of Molecules in Current Database

```python
# Get Full List of Molecules in Current Database
all_molecules = molecules.get_molecule_list()
print (all_molecules)
```

#### Get Functional Group Family of Desired Molecule

```python
# Get Functional Group Family of Desired Molecule
mol1= "benzene"
mol1_fg_family = molecules.get_molecule_fg_family(mol1)
print(mol1_fg_family)
```

#### Get Site Names in a Molecule for Interaction with Another Molecule

```python
# Get Site Names in a Molecule for Interaction with Another Molecule
mol1_atom_names = molecules.get_atom_names(mol1)
print(mol1_atom_names)
```

Output:

```
['RC1', 'H1']
```

#### Get Monomer Z-Matrix for Site of Interest

```python
# Get Monomer Z-Matrix for Site RC1
mol1_site_of_interest = "RC1"
mol1_zmat_for_site_of_interest = molecules.get_monomer_coordinates(mol1, mol1_site_of_interest)
print (mol1_zmat_for_site_of_interest)
```

#### Get Dimer Z-Matrix for Desired Interaction Type 

##### Interaction Type 1 (pi-stacking)

```python
# Get Benzene-Dimer Z-Matrix for pi- stacking
mol2 = 'benzene'
mol2_site_of_interest = 'RC1'
benzene_dimer_pistack = molecules.form_dimer_coordinates(mol1, mol1_site_of_interest, mol2, mol2_site_of_interest)
print (benzene_dimer_pistack)
```

##### Interaction Type 2 (C-H--pi interaction; T-shaped)

```python
# Get Benzene-Dimer Z-Matrix for C-H--pi interaction mode
mol2 = 'benzene'
mol2_site_of_interest = 'H1'
benzene_dimer_chpi = molecules.form_dimer_coordinates(mol1, mol1_site_of_interest, mol2, mol2_site_of_interest)
print (benzene_dimer_chpi)
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
