# Institut Pasteur Internship - IMOD Team

## Description
This repository contains all the scripts developed during my internship 
at **Institut Pasteur** within the **IMOD Team**.

The goal of this project is to evaluate the DrugCLIP program on _Escherichia coli data_.

---

## Author
- **Name :** Sixte Guillet de Chatellus
- **Internship period :** March - July 2026
- **Supervisor :** Dr. Christophe Zimmer
- **Team :** IMOD Team, Institut Pasteur, Paris

---

## 📁 Repository Structure

## Bash

### Preprocessing 

DrugCLIP env :

- **program_generate_files.sh :** generates folders of each PDB complex given a CSV of PDB complexes
- **program_conversion_PDB_smiles.sh :** extracts SMILES sequences from PDB
- **clean_smiles_files.sh :** cleans the SMILES sequences extracted from PDB

DeepCoy env:

- **program_prepare_data_DeepCoy.sh :** preprocessing for DeepCoy
- **program_generate_decoys.sh :** generates decoys given actives

DrugCLIP env:

- **program_conversion_actives_mol.sh :** converts actives files into mol.2 files
- **program_conversion_decoys_mol.sh :** converts decoys files into mol.2 files
- **script_clean_decoys.sh :** clean decoys files
- **program_conversion_ligand_crystalmol2.sh :** converts ligands files into crystalmol.2 files
- **script_remove_water.sh :** remove water molecules from crystal_ligand files
- **script_numerote_atoms.sh :** renumbers atoms of crystal_ligands files after removing water
- **program_write_dude_multi.sh :** preprocessing for DrugCLIP

### Run DrugCLIP

- **test_drugclip.sh :** Test DrugCLIP, specify the weights and the base folder for data (For DUD-E for example)
- **test_drugclip_fold.sh :** Test DrugCLIP on your own data, runs a 6 fold experiment
- **test_embedding.sh :** Test DrugCLIP on your own data, extracts the embeddings and creates a UMAP.

### Other bash

- **clean_actives_names:** clean the names of actives
- **drugclip.sh :** program of DrugCLIP
- **encode_mols.sh :** program used by DrugCLIP
- **encode_pocket.sh :** program used by DrugCLIP
- **lit_pcba.sh :** program used by DrugCLIP
- **lit_pcba_pockets.sh :** program used by DrugCLIP
- **pocket_size.sh :** outputs the average and maximum size of pockets
- **program_extract_file_multiple_ligands.sh :** Gathers the ligands from a same protein
- **run_compute_zscore_rank.sh :** study of the rank of each active
- **run_create_mols_lmdb_lit_pcba.sh :** preprocessing for LIT-PCBA
- **run_umap.sh :** create a UMAP given embeddings
- **script_count_actives.sh :** count the number of actives per protein
- **script_create_mols_alphafold.sh :** create mols.lmdb for the alphafold dataset
- **script_read_lmdb.sh :** script to analyse lmdb files

## Scripts
This file contains the scripts used by the bash files. Gather them in the same file for easier use.

---

## Githubs used in this work

- **DrugCLIP :** https://github.com/THU-ATOM/Drug-The-Whole-Genome/
- **DeepCoy :** https://github.com/fimrie/DeepCoy/
Both githubs required the installation of a dedicated environment in order to run the programs. Most programs can be run with the DrugCLIP environment. However the decoy generation requires the DeepCoy environment
