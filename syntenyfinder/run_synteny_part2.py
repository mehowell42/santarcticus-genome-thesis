#!/usr/bin/env python3
"""
run_synteny_part2.py

No-Jupyter version of Synteny_main.ipynb, part 2 of 3: runs OrthoFinder.
This is the slow step -- submit this as a proper batch job (sbatch/qsub,
whatever your HPC uses) rather than running it on a login node, so it
doesn't die if your connection drops.

Edit orthofinder_path and threads below, then submit/run:
    python run_synteny_part2.py
"""

import pickle

with open('synteny_checkpoint_1.pkl', 'rb') as f:
    synteny = pickle.load(f)

##########
## Settings -- EDIT THESE
##########

orthofinder_path = '/fast/jd139945/miniconda3/envs/genespace_env/bin/orthofinder'   # <- full path if it's not on your $PATH
threads = 24                        # <- match your job's allocation

##########
## Run OrthoFinder
##########

synteny.run_orthofinder(orthofinder_path, threads)
synteny.incorporate_orthofinder()
synteny.single_copy_orthologues()

with open('synteny_checkpoint_2.pkl', 'wb') as f:
    pickle.dump(synteny, f)

print("\nPart 2 complete. Next: python run_synteny_part3.py")
