 #!/bin/bash
python3 manage.py populate_db --delete_model PPI,PDI,Drug,Protein,Tissue,Disorder,PDiAssociations

python3 manage.py populate_db --data_dir . -p
python3 manage.py populate_db --data_dir . -exp gene_tissue_expression.gct

python3 manage.py populate_db --data_dir . -dr
#python3 manage.py populate_db --data_dir . -pdr drug-protein-interaction.txt
python3 manage.py populate_db -di
#python3 manage.py populate_db --data_dir . -pdi "" -ddi ""
#python3 manage.py populate_db -pp protein_protein_interaction_file.txt