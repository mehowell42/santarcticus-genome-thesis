#!/usr/bin/env python3
"""
Synteny_functions.py

Plain-Python conversion of dependencies/Synteny_functions.ipynb (identical
class, no notebook magic). Import from this instead of using Jupyter's
%run to load the same Synteny class on the command line.
"""

import subprocess
import os
import re
import pickle
import pandas as pd
from Bio import SeqIO


class Synteny:
    def __init__(self, root_directory, run_name, species_codes, proteome_ext = '.fa'):
        self.run_name = run_name
        self.dirs = dict()
        self.dirs['root'] = root_directory
        self.species_codes = species_codes
        self.species_data = {species:dict() for species in species_codes}
        self.check(proteome_ext)
        self.build()

    ###########
    ### Initialising functions
    ###########
    # Check to see if all the necessary files are correctly named and can be found.
    # Returns a dictionary dirs containing the appropriate file paths containing the input data
    def check(self, proteome_ext):
        self.add_dir('root', 'input_data')
        for type, extension in [('proteomes', proteome_ext),
                                ('genomes', '.fna'),
                                ('gene_rows', '.gtf')]:
            self.add_dir('input_data', type)

            # print(f"Checking files in {self.dirs[type]}...")
            contents = os.listdir(self.dirs[type])
            contents = [filename for filename in contents if filename.endswith(extension)]
            for species in self.species_codes:
                files = [filename for filename in contents if filename.startswith(species)]
                # print(files[0])
                if len(files) != 1:
                    print(f"Incorrect number of {type} files found for {species} \n{files}")
                else:
                    identifier = f"{species}_{type}"
                    self.dirs[identifier] = os.path.join(self.dirs[type], files[0])
        print('File check complete')
        return

    # Builds neccesary directories to run
    def build(self):
        try:
            self.make_dir('root', self.run_name)
            self.make_dir(self.run_name, 'run_files')
            self.make_dir(self.run_name, 'output')
            self.make_dir('run_files', 'run_proteomes')
            print('Directories built for run...')
        except:
            print('Adding existing directories...')
            self.add_dir('root', self.run_name)
            self.add_dir(self.run_name, 'run_files')
            self.add_dir(self.run_name, 'output')
            self.add_dir('run_files', 'run_proteomes')
        try:
            self.incorporate_orthofinder()
            print('OrthoFinder results found at ' + self.dirs['orthofinder_results'])
        except:
            print('No orthofinder results found .')
        self.dirs['save_synteny'] = os.path.join(self.dirs['run_files'], f"{self.run_name}.pkl")
        return

    ###########
    ### Preparing proteomes
    ###########
    # Read proteome files and save the SeqRecords into a list.
    # Return a dictionary sending each species code to its respective SeqRecords list.
    def read_proteomes(self):
        print('Reading proteomes...')
        for species in self.species_codes:
            filepath = self.dirs[f'{species}_proteomes']
            records = []
            with open(filepath, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    records.append(record)
            self.species_data[species]['proteomes'] = records
        return self.species_data[species]['proteomes'][:10]

    # For algs provided through adapted gene names in the proteome, parse the
    # algs and simplify the gene names prior to running orthofinder.
    def proteome_id_trim(self, species, position, delimiter = '_'):
        index = position - 1

        species_proteome = self.species_data[species]['proteomes']
        gene_ids = {record.id : record.id.split(delimiter)[index] for record in species_proteome}

        for gene in self.species_data[species]['proteomes']:
            gene.id = gene_ids[gene.id]

        return self.species_data[species]['proteomes'][:10]


    # Adds the three letter species code before each gene to make orthofinder output
    # easer to understand.
    # If no list of species codes is given, all proteomes are modified.
    def proteome_add_species(self, species_list = None):
        if species_list is None:
            species_list = self.species_codes
        try:
            for species in species_list:
                for gene in self.species_data[species]['proteomes']:
                    gene.id = species + '|' + gene.id
                    gene.description = ''
                print('Proteome gene names modified for ' + species)
            return self.species_data[species]['proteomes']
        except:
            print('Error occured in adding species codes to proteome.')
            print('Check that read_proteomes() has been run on all species')
            return

    # Writes proteomes selected for use to run_proteomes folder in run_files
    def write_proteomes(self):
        try:
            for species in self.species_codes:
                new_proteome_file = os.path.join(self.dirs['run_proteomes'], f"{species}.fa")
                self.dirs[f"{species}_run_proteomes"] = new_proteome_file
                with open(new_proteome_file, 'w') as output_handle:
                    SeqIO.write(self.species_data[species]['proteomes'], output_handle, 'fasta')
                print(f'{species} proteomes successfully written to {new_proteome_file}.')
            return self.species_data[species]['proteomes'][:10]
        except:
            print('Error occured in writing proteomes. Check that read_proteomes() has been run.')
            return

    # Function to primary_transcript.py as specified in the OrthoFinder tutorial?

    ###########
    ### Orthofinder
    ###########
    # Running orthofinder using nohup
    def run_orthofinder(self, orthofinder_path, threads):
        self.dirs['orthofinder_executable'] = orthofinder_path
        orthofinder_output = os.path.join(self.dirs['run_files'], 'orthofinder_output')
        command = ['nohup',
                   orthofinder_path,
                   '-f', self.dirs['run_proteomes'],
                   '-o', orthofinder_output,
                   '-t', str(threads)]
        print(command)

        result = subprocess.Popen(command, stdout=subprocess.PIPE)
        for line in result.stdout:
            print(line.decode(), end = '')
        result.wait()
        return

    # Incorporating Orthofinder results
    def incorporate_orthofinder(self):
        self.add_dir('run_files', 'orthofinder_output')
        results_dir = os.listdir(self.dirs['orthofinder_output'])[-1]
        self.dirs['orthofinder_results'] = os.path.join(self.dirs['orthofinder_output'], results_dir)
        self.add_dir('orthofinder_results', 'Orthogroups')
        self.dirs['orthogroups'] = os.path.join(self.dirs['Orthogroups'], 'Orthogroups.tsv')
        self.orthogroups = pd.read_csv(self.dirs['orthogroups'], sep = '\t')
        self.dirs['sco'] = os.path.join(self.dirs['Orthogroups'], 'Orthogroups_SingleCopyOrthologues.txt')
        return

    # Return a table of single copy orthologues for a subset of species
    # While it is possible to use the list of single copy orthologues generated by orthofinder,
    # this method allows for more flexibility and will result in a larger number of shared genes
    # if run with only a subset of the species used for orthofinder.
    def single_copy_orthologues(self, species_list = None):
        if species_list is None:
            species_list = self.species_codes
        orthogroup_subset = (self.orthogroups.copy()
                                   [['Orthogroup'] + species_list].dropna())
        multiple_copies = orthogroup_subset.map(lambda x: ',' in x)
        self.single_copy_orthogroups  = orthogroup_subset[~multiple_copies.any(axis = 1)]
        return self.single_copy_orthogroups

    ###########
    ### Building the karyotype file from genome
    ###########
    def build_karyotype(self, chromosomes_per_species):
        print('Reading karyotype files: this may take a moment...')
        for species in self.species_codes:
            species_karyotype = chromosomes_per_species[species]

            # Using SeqIO to read the genome file, get the scaffold id and length for each entry.
            genome_path = self.dirs[f'{species}_genomes']
            genome_records = list(SeqIO.parse(genome_path, 'fasta'))
            scaffolds = [{'Scaffold': record.id, 'Length': len(record.seq)} for record in genome_records]

            # Sort the scaffolds from longest to shortest and get the top scaffolds
            scaffolds = (pd.DataFrame(scaffolds).sort_values('Length', ascending = False)
                                                .reset_index(drop = True))
            chromosomes = scaffolds.head(species_karyotype)
            chromosomes.columns = ['Chromosome', 'Length']
            self.species_data[species]['karyotype'] = chromosomes
        return chromosomes

    def clean_karyotype(self, columns, labels = None):
        for species in self.species_codes:
            karyotype = self.species_data[species]['karyotype']
            cleaned_karyotype = pd.DataFrame(index = range(karyotype.shape[0]))
            for header in columns:
                if header in karyotype.columns:
                    cleaned_karyotype[header] = karyotype[header]
                elif header == 'SPECIES':
                    cleaned_karyotype[header] = species
                else:
                    cleaned_karyotype[header] = header
            if labels is not None:
                cleaned_karyotype.columns = labels
            self.species_data[species]['cleaned_karyotype'] = cleaned_karyotype
        return cleaned_karyotype

    def write_karyotype(self, index = False):
        for species in self.species_codes:
            file_path = os.path.join(self.dirs['output'], f"{species}_karyotype.txt")
            self.species_data[species]['cleaned_karyotype'].to_csv(file_path,
                                                                   sep = '\t',
                                                                   index = index)

    ###########
    ### Incorporating gene positions from GTF
    ###########
    # Query to keep only rows with certain key words
    def gtf_to_dataframe(self, species, annotation_type = None, feature = None, equivalence = ' '):
        gtf_filepath = self.dirs[f'{species}_gene_rows']
        columns = ["sequence", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        gtf_dataframe = pd.read_csv(gtf_filepath, sep = '\t',
                                    comment = '#',
                                    header = None,
                                    names = columns,
                                    dtype={'start': int, 'end': int})

        # Keep only rows of a given feature
        if feature is not None:
            gtf_dataframe = gtf_dataframe[gtf_dataframe.feature == feature]

        # Extract the specified annotation type from the 'attribute' column
        if annotation_type is None:
            gtf_dataframe['annotation'] = gtf_dataframe['attribute']
        else:
            gtf_dataframe['annotation'] = (gtf_dataframe['attribute'].str
                                          .extract(f'{annotation_type}{equivalence}"?([^";]+)"?', expand=False))
            gtf_dataframe = gtf_dataframe.dropna()
        gtf_dataframe = gtf_dataframe.drop_duplicates(subset = 'annotation', keep = 'first')
        self.species_data[species]['gene_rows'] = gtf_dataframe
        return gtf_dataframe

    # Add a species code with a pipe to the beginning of the chosen annotation names to match proteome gene name format
    def gtf_add_species(self):
        for species in self.species_codes:
            all_files_read = True
            if 'gene_rows' not in self.species_data[species]:
                print(f'gtf file for {species} not yet read.')
                all_files_read = False
        if all_files_read == True:
            for species in self.species_codes:
                new_annotations = species + '|' + self.species_data[species]['gene_rows']['annotation']
                self.species_data[species]['gene_rows']['annotation'] = new_annotations
        return new_annotations

    # Merges sco genes with the associated GTF information
    def merge_gtf(self):
        for species in self.species_codes:
            gene_rows = self.species_data[species]['gene_rows'].copy()
            orthogenes = self.single_copy_orthogroups[[species]].copy()
            orthogene_coords = pd.merge(orthogenes, gene_rows,
                               left_on = species, right_on = 'annotation',
                               how = 'left')
            self.species_data[species]['orthogene_coords'] = orthogene_coords
        return orthogene_coords


    ###########
    ### Incorporating lineage groups
    ###########
    # After merging the dataframes and creating the chromosome files, trace the genes based on
    # their chromosome in a given species
    def trace_chromosomes(self, species):
        groups = self.species_data[species]['cleaned_karyotype']['Chr'].tolist()
        print(groups)
        group_ids = [chr(ord('A') + i) for i in range(len(groups))]
        genes_to_groups = dict(zip(groups, group_ids))
        genes = self.species_data[species]['orthogene_coords']['sequence']
        gene_groups = [genes_to_groups[gene] if gene in genes_to_groups.keys() else gene for gene in genes ]
        self.trace = gene_groups
        return

    # For algs provided through adapted gene names in the proteome, parse the
    # algs and simplify the gene names prior to running orthofinder.
    def parse_algs(self, species, delimiter = '_', columns = ['index', 'alg', 'gene_id'],
                                                               mapping_add_species = True):
        try:
            alg_index = columns.index('alg')
        except:
            raise Exception('alg column not specified')
        try:
            gene_id_index = columns.index('gene_id')
        except:
            raise Exception('gene_id column not specified')

        filepath = self.dirs[f'{species}_proteomes']
        proteome_headers = []
        with open(filepath, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                proteome_headers.append(record.id)
        gene_headers = pd.DataFrame([record.split(delimiter) for record in proteome_headers])
        gene_headers.columns = columns
        if mapping_add_species:
            gene_headers['gene_id'] = species + '|' + gene_headers['gene_id']
        self.alg_mapping = dict(zip(gene_headers['gene_id'], gene_headers['alg']))
        return gene_headers

    def trace_algs(self, species):
        if not hasattr(self, 'alg_mapping'):
            raise Exception('No alg_mapping found. Please ensure parse_alg() has been run.')
        genes = self.species_data[species]['orthogene_coords']['annotation']
        algs = [self.alg_mapping[gene] if gene in self.alg_mapping.keys() else gene for gene in genes ]
        self.trace = algs
        return

    ###########
    ### Format and write coordinate files
    ###########
    # Reformat the merged coordinates dataframe.
    def clean_coords(self, columns, labels = None):
        for species in self.species_codes:
            coords = self.species_data[species]['orthogene_coords']
            cleaned_coords = pd.DataFrame(index = range(coords.shape[0]))
            for header in columns:
                if header == "TRACE":
                    cleaned_coords[header] = self.trace
                elif header in coords.columns:
                    cleaned_coords[header] = coords[header]
                else:
                    cleaned_coords[header] = header
            self.species_data[species]['cleaned_coords'] = cleaned_coords
            if labels is not None:
                cleaned_coords.columns = labels
        return cleaned_coords

    # Write the cleaned coordinates to a csv file
    def write_coords(self):
        for species in self.species_codes:
            file_path = os.path.join(self.dirs['output'], f"{species}_coordinates.tsv")
            self.species_data[species]['cleaned_coords'].to_csv(file_path,
                                                                sep = '\t',
                                                                header = False)

    ###########
    ### Compatibility helpers
    ###########
    # Given a species code and a substring which indicates a suffix, truncate the gene names in
    # the single_copy_orthogroups dataframe to omit the suffix.
    def truncate_sco(self, species, suffix):
        self.single_copy_orthogroups = (self.truncate(
                                                 self.single_copy_orthogroups,
                                                      species, suffix))
        return self.single_copy_orthogroups

    # Given a species code and a substring which indicates a suffix, truncate the annotations in
    # the gene_rows dataframe to omit the suffix.
    def truncate_gene_rows(self, species, suffix):
        self.species_data[species]['gene_rows'] = (self.truncate(
                                                                 self.species_data[species]['gene_rows'],
                                                                     'annotation', suffix))
        return self.species_data[species]['gene_rows']


    # Save the object to a file so that the analysis can be continuted from a later point
    def save_synteny(self):
        with open(self.dirs['save_synteny'], 'wb') as file:
            pickle.dump(self, file)

    ###########
    ### Helper functions
    ###########
    # Add a directory to dirs, assuming that the parent is already in dirs
    def add_dir(self, parent, add):
        new_path = os.path.join(self.dirs[parent], add)
        if not os.path.exists(new_path):
            raise Exception(f"Expected directory at {new_path} was not found.")

        else:
            self.dirs[add] = new_path
        return

    # Create a new directory and add it to dirs
    def make_dir(self, parent, make):
        new_dir = os.path.join(self.dirs[parent], make)
        os.makedirs(new_dir)
        self.dirs[make] = new_dir
        return

    # Given a dataframe, a column, and a suffix, trim the suffix from the column entries and
    # return the altered dataframe.
    def truncate(self, dataframe, column, suffix):
        dataframe = dataframe.copy()
        entries_list = dataframe[column].tolist()
        truncated = [self.truncate_gene(entry, suffix) for entry in entries_list]
        dataframe[column] = truncated
        return dataframe

    def truncate_gene(self, gene, suffix):
        try:
            index = gene.find(suffix)
            return(gene[:index] if index != -1 else gene)
        except:
            print('Issue detected in trucating gene.')
            return(gene)

def load_synteny(root_directory, run_name):
    saved_file = os.path.join(root_directory, run_name, 'run_files', f'{run_name}.pkl')
    with open(saved_file, 'rb') as file:
        synteny = pickle.load(file)
    return(synteny)
