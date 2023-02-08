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

### Data folder
Static datasets are mounted from a directory now, instead of fusing them into the image. Download them from the following link and put them into the data folder that is mounted by the docker-compose.yml:
https://wolken.zbh.uni-hamburg.de/index.php/s/gywnL3HP26CWrgA
