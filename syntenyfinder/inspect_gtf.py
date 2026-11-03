#!/usr/bin/env python3
"""
inspect_gtf.py

For each species, gtf_to_dataframe() needs to know which attribute holds
the gene ID (annotation_type, e.g. 'gene_id' or 'transcript_id') and which
feature rows to keep (feature, e.g. 'CDS'). These differ between RefSeq,
Ensembl, and custom annotations, so check each species before committing to
settings in run_synteny_part3.py.

Usage (one species at a time):
    python inspect_gtf.py Bfl
    python inspect_gtf.py Aam
    ...
"""

import pickle
import sys

if len(sys.argv) != 2:
    sys.exit("Usage: python inspect_gtf.py <species_code>")

species = sys.argv[1]

with open('synteny_checkpoint_2.pkl', 'rb') as f:
    synteny = pickle.load(f)

df = synteny.gtf_to_dataframe(species)  # no extra args -- raw attribute column
print(f"\n--- {species}: feature types present ---")
print(df['feature'].value_counts())
print(f"\n--- {species}: first 3 raw attribute strings ---")
for a in df['attribute'].head(3):
    print(a)
print(f"\nLook for the ID you want (gene_id=\"...\" or ID=... etc) in the "
      f"attribute strings above, and whether '=' or ' ' separates key/value. "
      f"Use that to fill in this species' entry in run_synteny_part3.py.")
