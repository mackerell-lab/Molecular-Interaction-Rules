import psi4
import numpy as np

# configurations
# --------------

psi4.core.set_num_threads(8)
psi4.set_memory('8000mb')
psi4.set_options({
  'scf_type': 'df',
  'g_convergence': 'gau_tight',
  'reference': 'rhf',
  'freeze_core': 'true',
})

monomer_a = '''\
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
'''

monomer_b = '''\
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
'''

dimer = '''\
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
X21   X11  DISTANCE  C11   180.0000  C12   90.0000
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
'''

if __name__ == '__main__':

  distances = np.arange(1.0, 8.1, 0.1)
  
  # Load Monomers

  monomer_a_molecule = psi4.geometry(monomer_a)
  monomer_b_molecule = psi4.geometry(monomer_b)
  
  # Calculate Monomer Energies

  monomer_energy_a = psi4.energy(
                      theory='mp2/aug-cc-pvdz',
                      molecule=monomer_a_molecule,
                      bsse='cp'
                    )
  
  monomer_energy_b = psi4.energy(
                      theory='mp2/aug-cc-pvdz',
                      molecule=monomer_b_molecule,
                      bsse='cp'
                    )
  
  monomer_energy_a = monomer_energy_a * psi4.constants.hartree2kcalmol
  monomer_energy_b = monomer_energy_b * psi4.constants.hartree2kcalmol

  interaction_energies = []    

  # Run Potential Energy Scan
  
  for distance in distances:

      dimer_zmatrix = dimer.replace('DISTANCE', str(distance))
      dimer_molecule = psi4.geometry(monomer_b)
      dimer_molecule.update_geometry()
      dimer_molecule.print_out()

      dimer_energy = psi4.energy(
                theory='mp2/aug-cc-pvdz',
                molecule=dimer_molecule,
                bsse='cp'
      )

      dimer_energy = dimer_energy * psi4.constants.hartree2kcalmol

      interaction_energy = dimer_energy - monomer_energy_a - monomer_energy_b
      interaction_energies.append(interaction_energy)

  print ('Interaction Energies: %s' % interaction_energies) 
