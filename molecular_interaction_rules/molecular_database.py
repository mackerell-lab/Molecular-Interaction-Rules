#!/usr/bin/env python3
#
# Molecular Interation Rules: MoleculerDatabase
# ---------------------------------------------

# Imports
# -------

import os, sys

# Reconfigurations
# ----------------

if os.name == 'nt':
   sys.stdin.reconfigure(encoding='utf-8')
   sys.stdout.reconfigure(encoding='utf-8')

# Molecules

## Alcohol

from molecular_interaction_rules.molecules.non_aromatic.alcohols.methanol import Methanol

## Alkanes

from molecular_interaction_rules.molecules.non_aromatic.alkanes.cyclobutane import Cyclobutane
from molecular_interaction_rules.molecules.non_aromatic.alkanes.cyclohexane import Cyclohexane
from molecular_interaction_rules.molecules.non_aromatic.alkanes.cyclopropane import Cyclopropane
from molecular_interaction_rules.molecules.non_aromatic.alkanes.neopentane import Neopentane
from molecular_interaction_rules.molecules.non_aromatic.alkanes.propane import Propane

## Alkenes

from molecular_interaction_rules.molecules.non_aromatic.alkenes.cyclohexene import Cyclohexene
from molecular_interaction_rules.molecules.non_aromatic.alkenes.cyclopentene import Cyclopentene
from molecular_interaction_rules.molecules.non_aromatic.alkenes.methoxyethene import Methoxyethene
from molecular_interaction_rules.molecules.non_aromatic.alkenes.one_three_dibutene import OneThreeDibutene
from molecular_interaction_rules.molecules.non_aromatic.alkenes.propene import Propene
from molecular_interaction_rules.molecules.non_aromatic.alkenes.two_pyrroline import TwoPyrroline

## Alkynes

from molecular_interaction_rules.molecules.non_aromatic.alkynes.propyne import Propyne

## Amides

from molecular_interaction_rules.molecules.non_aromatic.amides.acetamide import Acetamide
from molecular_interaction_rules.molecules.non_aromatic.amides.amidinium import Amidinium
from molecular_interaction_rules.molecules.non_aromatic.amides.azetidinone import Azetidinone
from molecular_interaction_rules.molecules.non_aromatic.amides.dimethylformamide import Dimethylformamide
from molecular_interaction_rules.molecules.non_aromatic.amides.methylacetamide import MethylAcetamide
from molecular_interaction_rules.molecules.non_aromatic.amides.acetamide import Acetamide
from molecular_interaction_rules.molecules.non_aromatic.amides.methylacetamide import MethylAcetamide
from molecular_interaction_rules.molecules.non_aromatic.amides.prolineamide import Prolineamide
from molecular_interaction_rules.molecules.non_aromatic.amides.prolineamide_charged import ProlineamideCharged
from molecular_interaction_rules.molecules.non_aromatic.amides.two_pyrrolidinone import TwoPyrrolidinone

## Amines

from molecular_interaction_rules.molecules.non_aromatic.amines.ammonia import Ammonia
from molecular_interaction_rules.molecules.non_aromatic.amines.dimethylamine import Dimethylamine
from molecular_interaction_rules.molecules.non_aromatic.amines.ethoxy_guanidine import EthoxyGuanidine
from molecular_interaction_rules.molecules.non_aromatic.amines.ethyl_ammonium import EthylAmmonium
from molecular_interaction_rules.molecules.non_aromatic.amines.hydrazine import Hydrazine
from molecular_interaction_rules.molecules.non_aromatic.amines.methylamine import Methylamine
from molecular_interaction_rules.molecules.non_aromatic.amines.piperidine import Piperidine
from molecular_interaction_rules.molecules.non_aromatic.amines.tetramethylammonium import Tetramethylammonium
from molecular_interaction_rules.molecules.non_aromatic.amines.trimethylamine import Trimethylamine
from molecular_interaction_rules.molecules.non_aromatic.amines.trimethylammonium import Trimethylammonium

## Carbonyls

from molecular_interaction_rules.molecules.non_aromatic.carbonyls.acetaldehyde import Acetaldehyde
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.acetate import Acetate
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.acetic_acid import AceticAcid
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.acetone import Acetone
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.carbon_dioxide import CarbonDioxide
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.formaldehyde import Formaldehyde
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.methylacetate import MethylAcetate
from molecular_interaction_rules.molecules.non_aromatic.carbonyls.urea import Urea

## Ethers

from molecular_interaction_rules.molecules.non_aromatic.ethers.dimethylether import Dimethylether
from molecular_interaction_rules.molecules.non_aromatic.ethers.epoxide import Epoxide
from molecular_interaction_rules.molecules.non_aromatic.ethers.oxetane import Oxetane
from molecular_interaction_rules.molecules.non_aromatic.ethers.tetrahydrofuran import Tetrahydrofuran
from molecular_interaction_rules.molecules.non_aromatic.ethers.tetrahydropyran import Tetrahydropyran

## Halogens

from molecular_interaction_rules.molecules.non_aromatic.halogens.bromoethane import Bromoethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.chloroethane import Chloroethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.dibromoethane import Dibromoethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.difluoroethane import Difluoroethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.fluoroethane import Fluoroethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.tribromoethane import Tribromoethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.trichloroethane import Trichloroethane
from molecular_interaction_rules.molecules.non_aromatic.halogens.trifluoroethane import Trifluoroethane

## Imines

from molecular_interaction_rules.molecules.non_aromatic.imines.ethenamine import Ethenamine

## Nitriles

from molecular_interaction_rules.molecules.non_aromatic.nitriles.acetonitrile import Acetonitrile

## Organophosphorous

from molecular_interaction_rules.molecules.non_aromatic.organophosphorous.dimethyl_phosphate import DimethylPhosphate
from molecular_interaction_rules.molecules.non_aromatic.organophosphorous.methyl_phosphate import MethylPhosphate

## Organosulfur

from molecular_interaction_rules.molecules.non_aromatic.organosulfur.dimethyl_sulfone import DimethylSulfone
from molecular_interaction_rules.molecules.non_aromatic.organosulfur.dimethyl_sulfoxide import DimethylSulfoxide
from molecular_interaction_rules.molecules.non_aromatic.organosulfur.dimethyl_trithiocarbonate import DimethylTrithiocarbonate
from molecular_interaction_rules.molecules.non_aromatic.organosulfur.dimethyldisulfide import DimethylDisulfide
from molecular_interaction_rules.molecules.non_aromatic.organosulfur.ethylsulfanyl_phosphonic_acid import EthylSulfanylPhosphonicAcid
from molecular_interaction_rules.molecules.non_aromatic.organosulfur.methanethiol import MethaneThiol
from molecular_interaction_rules.molecules.non_aromatic.organosulfur.methylthiolate import MethylThiolate

## Aromatics

from molecular_interaction_rules.molecules.aromatic.azulene import Azulene
from molecular_interaction_rules.molecules.aromatic.benzene import Benzene
from molecular_interaction_rules.molecules.aromatic.bipyrrole import Bipyrrole
from molecular_interaction_rules.molecules.aromatic.bromobenzene import BromoBenzene
from molecular_interaction_rules.molecules.aromatic.chlorobenzene import ChloroBenzene
from molecular_interaction_rules.molecules.aromatic.fluorobenzene import FluoroBenzene
from molecular_interaction_rules.molecules.aromatic.four_pyridinone import FourPyridinone
from molecular_interaction_rules.molecules.aromatic.furan import Furan
from molecular_interaction_rules.molecules.aromatic.imidazole import Imidazole
from molecular_interaction_rules.molecules.aromatic.imidazolium import Imidazolium
from molecular_interaction_rules.molecules.aromatic.indole import Indole
from molecular_interaction_rules.molecules.aromatic.indolizine import Indolizine
from molecular_interaction_rules.molecules.aromatic.iodobenzene import IodoBenzene
from molecular_interaction_rules.molecules.aromatic.isoxazole import Isoxazole
from molecular_interaction_rules.molecules.aromatic.methyleneoxindole import Methyleneoxindole
from molecular_interaction_rules.molecules.aromatic.nitrobenzene import Nitrobenzene
from molecular_interaction_rules.molecules.aromatic.one_phenyl_four_pyridinone import OnePhenylFourPyridinone
from molecular_interaction_rules.molecules.aromatic.phenol import Phenol
from molecular_interaction_rules.molecules.aromatic.phenoxazine import Phenoxazine
from molecular_interaction_rules.molecules.aromatic.pyridine import Pyridine
from molecular_interaction_rules.molecules.aromatic.pyridinium import Pyridinium
from molecular_interaction_rules.molecules.aromatic.pyrimidine import Pyrimidine
from molecular_interaction_rules.molecules.aromatic.pyrrolidine import Pyrrolidine
from molecular_interaction_rules.molecules.aromatic.thiophene import Thiophene
from molecular_interaction_rules.molecules.aromatic.three_amino_pyridine import ThreeAminoPyridine
from molecular_interaction_rules.molecules.aromatic.two_h_pyran import TwoHPyran
from molecular_interaction_rules.molecules.aromatic.uracil import Uracil

class MoleculerDatabase(object):

    def __init__(self):

        self.master_list = self.build_molecular_database()

    def build_molecular_database(self):

        molecules = [
          
            # Alcohols

            Methanol,
          
            # Alkanes

            Cyclobutane,
            Cyclohexane,
            Cyclopropane,
            Neopentane,
            Propane,
          
            # Alkenes

            Cyclohexene,
            Cyclopentene,
            Methoxyethene,
            OneThreeDibutene,
            Propene,
            TwoPyrroline,
          
            # Alkynes

            Propyne,
            
            # Amides

            Acetamide,
            Amidinium,
            Azetidinone,
            Dimethylformamide,
            MethylAcetamide,
            Acetamide,
            MethylAcetamide,
            Prolineamide,
            ProlineamideCharged,
            TwoPyrrolidinone,
          
            # Amines

            Ammonia,
            Dimethylamine,
            EthoxyGuanidine,
            EthylAmmonium,
            Hydrazine,
            Methylamine,
            Piperidine,
            Tetramethylammonium,
            Trimethylamine,
            Trimethylammonium,
          
            # Carbonyls

            Acetaldehyde,
            Acetate,
            AceticAcid,
            Acetone,
            CarbonDioxide,
            Formaldehyde,
            MethylAcetate,
            Urea,
            
            # Ethers

            Dimethylether,
            Epoxide,
            Oxetane,
            Tetrahydrofuran,
            Tetrahydropyran,
            
            # Halogens

            Bromoethane,
            Chloroethane,
            Dibromoethane,
            Difluoroethane,
            Fluoroethane,
            Tribromoethane,
            Trichloroethane,
            Trifluoroethane,
            
            # Imines
  
            Ethenamine,
            
            # Nitriles
  
            Acetonitrile,
            
            # Organophosphorous
  
            DimethylPhosphate,
            MethylPhosphate,
            
            # Organosulfur
  
            DimethylSulfone,
            DimethylSulfoxide,
            DimethylTrithiocarbonate,
            DimethylDisulfide,
            EthylSulfanylPhosphonicAcid,
            MethaneThiol,
            MethylThiolate,
          
            # Aromatics
          
            Azulene,
            Benzene,
            Bipyrrole,
            BromoBenzene,
            ChloroBenzene,
            FluoroBenzene,
            FourPyridinone,
            Furan,
            Imidazole,
            Imidazolium,
            Indole,
            Indolizine,
            IodoBenzene,
            Isoxazole,
            Methyleneoxindole,
            Nitrobenzene,
            OnePhenylFourPyridinone,
            Phenol,
            Phenoxazine,
            Pyridine,
            Pyridinium,
            Pyrimidine,
            Pyrrolidine,
            Thiophene,
            ThreeAminoPyridine,
            TwoHPyran,
            Uracil

        ]

        return molecules

    def get_atom_names(self, monomer):

        '''

        Retrieve Atom Names of Molecule

        '''

        atom_names = []

        for molecule in self.master_list:

            name = str(molecule.__name__).lower()

            if name == monomer:
                molecule = molecule()
                atom_names = list(molecule.get_monomer_a_species().keys())

        return atom_names

    def get_monomer_coordinates(
      self,
      monomer,
      atom_name,
      monomer_b=False,
      functional_group_family=False
    ):


      '''

      Get the Coordinates of a Monomer

      '''

      coordinates = ''

      for molecule in self.master_list:

        name = str(molecule.__name__).lower()

        if name == monomer:

            molecule = molecule()
            fg_family = molecule.__module__.split('.')[-2]

            if functional_group_family:
                return fg_family
            try:
              if monomer_b:
                coordinates = molecule.get_monomer_b_species()[atom_name][0]
                rank_order = molecule.get_monomer_b_species()[atom_name][1]
              else:
                coordinates = molecule.get_monomer_a_species()[atom_name][0]
                rank_order = molecule.get_monomer_a_species()[atom_name][1]

            except KeyError as e:
                print ('Missing Entry: %s | %s' % (name, atom_name))
                return coordinates

            return coordinates

    def form_dimer_coordinates(
      self,
      molecule_name_1,
      atom_name_1,
      molecule_name_2,
      atom_name_2,
    ):

      '''

      Get the Dimer Coordinates

      '''

      monomer_a = self.get_monomer_coordinates(
        molecule_name_1,
        atom_name_1,
      )

      monomer_b = self.get_monomer_coordinates(
        molecule_name_2,
        atom_name_2,
        monomer_b=True,
      )

      # Fetch First 3 Atoms from Monomer A

      monomer_atom_one = monomer_a.split('\n')[0]
      monomer_atom_two = monomer_a.split('\n')[1].split()[0]
      monomer_atom_three = monomer_a.split('\n')[2].split()[0]

      # Replace Monomer B Connections with Atoms from Monomer A

      monomer_b = monomer_b.replace(':1', monomer_atom_one)
      monomer_b = monomer_b.replace(':2', monomer_atom_two)
      monomer_b = monomer_b.replace(':3', monomer_atom_three)

      dimer = monomer_a + '--\n' + monomer_b

      return dimer

