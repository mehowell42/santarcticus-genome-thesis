#!/usr/bin/env python3
"""
run_synteny_part3_bfl_pma.py

Part 3 for a run containing BOTH Bfl and Pma as ALG-label sources.

  - Bfl labels come from the Simakov table (bfl_alg_lookup.tsv), broad
    coverage but NO R genes.
  - Pma labels come from the bryozoan paper's Supplemental Data S1
    (pma_alg_lookup.tsv), 995 genes incl. 11 R genes.

Both anchor proteomes must have been rewritten with headers
"index|ALG|accession" (build_bfl_alg_proteome.py works for either lookup)
BEFORE running setup/part1/part2, exactly as was done for Bfl.

For each orthogroup, the final label is:
   Bfl label if present, else Pma label if present, else NA.
Where BOTH have a label and they disagree, the orthogroup keeps the Bfl
label and is counted as a conflict (printed at the end) -- worth looking at.

    python run_synteny_part3_bfl_pma.py
"""

import pickle
from collections import Counter
from Bio import SeqIO

with open('synteny_checkpoint_2.pkl', 'rb') as f:
    synteny = pickle.load(f)

##########
## Settings
##########

ALG_ANCHORS = ['Bfl', 'Pma']   # species whose proteome headers carry ALG labels

GTF_SETTINGS = {
    'Bfl': dict(annotation_type='protein_id', feature='CDS', equivalence='='),
    'Pma': dict(annotation_type='protein_id', feature='CDS', equivalence='='),
    'Aam': dict(annotation_type='ID', feature='mRNA', equivalence='='),
    'Aja': dict(annotation_type='ID', feature='mRNA', equivalence='='),
'Abi': dict(annotation_type='ID', feature='mRNA', equivalence='='),
    'Afi': dict(annotation_type='ID', feature='mRNA', equivalence='='),
    'Spu': dict(annotation_type='ID', feature='mRNA', equivalence='='),
    'San': dict(annotation_type='ID', feature='mRNA', equivalence='='),
}

##########
## Parse gene coordinates for every species
##########

for species, kwargs in GTF_SETTINGS.items():
    print(f"Parsing {species} with {kwargs} ...")
    df = synteny.gtf_to_dataframe(species, **kwargs)
    print(f"  -> {len(df)} rows with a matched {kwargs['annotation_type']}")

##########
## ANCHOR FIX -- for each anchor, remap the extracted protein_id back to the
## exact "index|ALG|accession" header OrthoFinder saw, BEFORE
## gtf_add_species() prefixes everything with "<species>|".
##########

for anchor in ALG_ANCHORS:
    accession_to_header = {}
    with open(synteny.dirs[f'{anchor}_proteomes']) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            parts = record.id.split('|')
            if len(parts) == 3:
                accession_to_header[parts[2]] = record.id

    rows = synteny.species_data[anchor]['gene_rows']
    n_hit = rows['annotation'].isin(accession_to_header).sum()
    print(f"\n{anchor} anchor fix: {n_hit} of {len(rows)} extracted accessions "
          f"matched a proteome header")
    if n_hit < 0.5 * len(rows):
        print(f"  WARNING: under half matched for {anchor}. Check that its GFF "
              f"and proteome came from the same annotation.")
    rows['annotation'] = rows['annotation'].map(accession_to_header).fillna(rows['annotation'])

synteny.gtf_add_species()
synteny.merge_gtf()

##########
## Build ALG labels from each anchor, then combine
##########

def labels_from(anchor):
    """Species-prefixed annotation is 'Sp|index|ALG|accession' -> ALG is field 3."""
    out = []
    for v in synteny.species_data[anchor]['orthogene_coords']['annotation']:
        parts = str(v).split('|')
        out.append(parts[2] if len(parts) == 4 else 'NA')
    return out

per_anchor = {a: labels_from(a) for a in ALG_ANCHORS}
for a, labs in per_anchor.items():
    n = sum(1 for x in labs if x != 'NA')
    print(f"\n{a}: {n} of {len(labs)} orthogroups carry a real ALG label")

n_orthogroups = len(per_anchor[ALG_ANCHORS[0]])
combined = []
stats = Counter()
for i in range(n_orthogroups):
    labs = [per_anchor[a][i] for a in ALG_ANCHORS]
    real = [x for x in labs if x != 'NA']
    if not real:
        combined.append('NA')
        stats['unlabelled'] += 1
    elif len(set(real)) == 1:
        combined.append(real[0])
        stats['labelled (anchors agree or only one labelled)'] += 1
    else:
        combined.append(real[0])          # first anchor listed wins (Bfl)
        stats['CONFLICT (anchors disagree)'] += 1

synteny.trace = combined

print("\nCombined ALG labels:")
for k, v in stats.items():
    print(f"  {k}: {v}")
print(f"  R-labelled orthogroups: {combined.count('R')}")

##########
## Write final output
##########

columns = ['Complete', 'sequence', 'start', 'end', 'TRACE']
synteny.clean_coords(columns)
synteny.write_coords()

print("\nDone. Check the output/ folder for <Code>_coordinates.tsv files.")
