#!/bin/bash
python3 manage.py populate_db --delete_model PPI,PDI,Drug,Protein,Tissue
python3 manage.py populate_db -p .
python3 manage.py populate_db --data_dir . --exp_file gene_tissue_expression.gct
python3 manage.py populate_db --data_dir . --drug_file drug-file.txt
python3 manage.py populate_db -pp .
python3 manage.py populate_db -pd .