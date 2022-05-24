#The new data folder is data_drugstone All the necessary files are there.

```bash
python3 manage.py migrate
python3 manage.py createfixtures

python3 manage.py populate_db --delete_model PPI,PDI,Drug,Protein,Tissue,Disorder,PDiAssociations

python3 manage.py populate_db --data_dir . -p protein-file.txt

python3 manage.py populate_db --data_dir . -exp gene_tissue_expression.gct

python3 manage.py populate_db --data_dir . -dr drug-file.txt -pdr drug-protein-interaction.txt

python3 manage.py populate_db --data_dir . -di "" -pdi "" -ddi ""

python3 manage.py populate_db --data_dir . -pp protein_protein_interaction_file.txt

python3 manage.py make_graphs

```

### Docker PROD environment (building is optional)
``docker-compose up --build``


### Docker DEV environment (building is optional)
``docker-compose -f docker-compose.yml up -d --build``
