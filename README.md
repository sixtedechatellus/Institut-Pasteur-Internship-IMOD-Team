# Institut Pasteur Internship - IMOD Team

## Description
This repository contains all the scripts developed during my internship 
at **Institut Pasteur** within the **IMOD Team**.

The goal of this project is to evaluate the DrugCLIP program on _Escherichia coli data_ 

---

## Author
- **Name :** Sixte Guillet de Chatellus
- **Internship period :** March - July 2026
- **Supervisor :** Dr. Christophe Zimmer
- **Team :** IMOD Team, Institut Pasteur, Paris

---

## 📁 Repository Structure

# Bash
- **clean_actives_names.sh : clean the names of actives
- **clean_smiles_files.sh : cleans the SMILES sequences extracted from PDB
- **drugclip.sh : program of DrugCLIP
- **encode_mols.sh : program used by DrugCLIP
- **encode_pocket.sh : program used by DrugCLIP
- **lit_pcba.sh : program used by DrugCLIP
- **lit_pcba_pockets.sh : program used by DrugCLIP
- **pocket_size.sh : outputs the average and maximum size of pockets
- **program_conversion_PDB_smiles.sh : extracts SMILES sequences from PDB
- **program_conversion_actives_mol.sh : converts actives files into mol.2 files
- **program_conversion_decoys_mol.sh : converts decoys files into mol.2 files
- **program_conversion_ligand_crystalmol2.sh : converts ligands files into crystalmol.2 files
- **program_extract_file_multiple_ligands.sh : extracts ligands PDB from the complex PDB
- **program_generate_decoys.sh : generates decoys given actives
- **program_generate_files.sh : generates folders of each PDB complex given a CSV of PDB complexes
- **program_prepare_data_DeepCoy.sh : preprocessing for DeepCoy
- **program_write_dude_multi.sh : preprocessing for DrugCLIP
- **run_compute_zscore_rank.sh : study of the rank of each active
- **run_create_mols_lmdb_lit_pcba.sh : preprocessing for LIT-PCBA
- **run_umap.sh : create a UMAP given embeddings
- **script_clean_decoys.sh : clean decoys files
- **script_count_actives.sh : count the number of actives per protein
- **script_create_mols_alphafold.sh : 
- **script_numerote_atoms.sh
- **script_read_lmdb.sh
- **script_remove_water.sh
- **test_drugclip.sh
- **test_drugclip_fold.sh
- **test_drugclip_fold_dedicated.sh
- **test_embedding.sh
