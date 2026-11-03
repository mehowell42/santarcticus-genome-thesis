#!/usr/bin/env python3
"""
run_synteny_part1.py

No-Jupyter version of Synteny_main.ipynb, part 1 of 3.
Runs setup through write_proteomes(), then stops -- deliberately, so you can
check everything looks right before kicking off the slow OrthoFinder step
(part 2), which you'll likely want to submit as a proper batch job rather
than run inline.

Edit the SETTINGS block below, then:
    python run_synteny_part1.py
"""

import pickle
from Synteny_functions import Synteny

##########
## Settings -- EDIT THESE
##########

root_directory = './my_synteny_run'   # <- whatever you passed as --root earlier
species_codes = ['Bfl', 'Pma', 'Aam', 'Aja', 'Abi', 'Afi', 'Spu', 'San']
run_name = 'my_run'                   # <- choose a name

##########
## Start
##########

synteny = Synteny(root_directory, run_name, species_codes, proteome_ext='.fasta')

##########
## Karyotype files
##########

chromosomes = {
    'Bfl': 19, 'Pma': 19, 'Aam': 22, 'Aja': 23, 'Abi': 11, 'Afi': 20, 'Spu': 21, 'San': 20,
}
synteny.build_karyotype(chromosomes)

columns = ['Chromosome', '1', 'Length', 'SPECIES', '12', '25252']
labels = ['Chr', 'Start', 'End', 'species', 'size', 'color']
synteny.clean_karyotype(columns, labels)
synteny.write_karyotype()

##########
## Proteomes
##########

synteny.read_proteomes()
synteny.proteome_add_species()
synteny.write_proteomes()

##########
## Save progress so part 2 can pick up from exactly here
##########

with open('synteny_checkpoint_1.pkl', 'wb') as f:
    pickle.dump(synteny, f)

print("\nPart 1 complete. Check the output above for errors, then run:")
print("    python run_synteny_part2.py")
