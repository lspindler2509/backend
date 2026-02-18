import os
import requests
from drugstone.util.query_db import query_proteins_by_identifier
from urllib.parse import urlparse, parse_qs
from drugstone.settings import DEBUG

def download_geneset(url, save_path):
    os.makedirs(save_path, exist_ok=True)
    # libraryName aus URL extrahieren
    parsed_url = urlparse(url)
    query_params = parse_qs(parsed_url.query)
    library_name = query_params.get("libraryName", ["geneset"])[0]
    filename = f"{library_name}.txt"

    file_path = os.path.join(save_path, filename)
    response = requests.get(url)
    response.raise_for_status()
    with open(file_path, 'w') as f:
        f.write(response.text)
    print(f"Downloaded {filename} to {file_path}")
    return file_path, library_name

def create_file(filename, data, path_genesets):
    os.makedirs(path_genesets, exist_ok=True)
    file_path = os.path.join(path_genesets, filename)
    with open(file_path, 'a') as f:
        print("Create file: ", file_path)
        for pathway, genes in data.items():
            f.write("{}\t{}\n".format(pathway, '\t'.join(genes)))

# the files to be parsed have to be in data/gene_sets
# result files will be in data/gene_sets/new_genesets to not overwrite the original files
# source of the files: https://maayanlab.cloud/Enrichr/#libraries
def parse_genesets(kegg_url, reactome_url, wiki_url, reviewed):
    root_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    path_genesets = os.path.join(root_path, 'data/gene_sets')

    kegg_file, _ = download_geneset(kegg_url, path_genesets)
    reactome_file, _ = download_geneset(reactome_url, path_genesets)
    wiki_file, _ = download_geneset(wiki_url, path_genesets)

    gene_sets = []

    for file_path in [kegg_file, reactome_file, wiki_file]:
        pathway_dict = {}
        with open(file_path, "r") as f:
            for line in f:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = [gene for gene in parts[1:] if gene]
                pathway_dict[pathway] = genes
        gene_sets.append(pathway_dict)

    genesets_new = []
    for geneset in gene_sets:
        if DEBUG:
            print("Query proteins for set: ", len(geneset))
        entrez, symbol, uniprot, ensembl = {}, {}, {}, {}
        for i, pathway in enumerate(geneset.keys(), 1):
            if DEBUG:
                print("Query proteins for pathway: ", len(pathway), " pathway: ", i)
            nodes_mapped, _ = query_proteins_by_identifier(geneset[pathway], "symbol", reviewed)
            entrez[pathway] = set()
            symbol[pathway] = set()
            uniprot[pathway] = set()
            ensembl[pathway] = set()
            for node in nodes_mapped:
                entrez[pathway].update(node["entrez"])
                symbol[pathway].update(node["symbol"])
                uniprot[pathway].update(node["uniprot"])
                if "ensg" in node:
                    ensembl[pathway].update(node["ensg"])
        genesets_new.append({"symbol": symbol, "entrez": entrez, "uniprot": uniprot, "ensembl": ensembl})

    print("Create pathway files")
    for i, d in enumerate(genesets_new, 1):  # 1 → kegg, 2 → reactome, 3 → wiki
        for key, value in d.items():
            for pathway, genes in value.items():
                if i == 1:
                    geneset = "kegg"
                elif i == 2:
                    geneset = "reactome"
                elif i == 3:
                    geneset = "wiki"
                if reviewed:
                    filename = f"{geneset}_{key}_reviewed.txt"
                else:
                    filename = f"{geneset}_{key}.txt"
                create_file(filename, {pathway: genes}, path_genesets)
    
    
    