
import os
from drugstone.util.query_db import (
    query_proteins_by_identifier,
)

def create_file(filename, data, path_genesets):
    path_new_genesets = os.path.join(path_genesets, "new_genesets")
    if not os.path.exists(path_new_genesets):
        os.makedirs(path_new_genesets)

    file_path = os.path.join(path_new_genesets, filename)
    with open(file_path, 'a') as f:
        print("Create file: ", file_path)
        for pathway, genes in data.items():
            f.write("{}\t{}\n".format(pathway, '\t'.join(genes)))

# the files to be parsed have to be in data/gene_sets
# result files will be in data/gene_sets/new_genesets to not overwrite the original files
def parse_genesets(kegg_filename, reactome_filename, wiki_filename, reviewed):
    root_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    path_genesets = os.path.join(root_path, 'data/gene_sets')
    gene_sets = []
    
    # parse KEGG
    pathway_kegg = {}
    with open(os.path.join(path_genesets, kegg_filename), "r") as file:
        for line in file:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = parts[1:]
            genes = [gene for gene in genes if gene]
            pathway_kegg[pathway] = genes
    gene_sets.append(pathway_kegg)
    
    # parse reactome
    pathway_reactome = {}
    with open(os.path.join(path_genesets, reactome_filename), "r") as file:
        for line in file:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = parts[1:]
            genes = [gene for gene in genes if gene]
            pathway_reactome[pathway] = genes
    gene_sets.append(pathway_reactome)
    
     # parse wiki
    pathway_wiki = {}
    with open(os.path.join(path_genesets, wiki_filename), "r") as file:
        for line in file:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = parts[1:]
            genes = [gene for gene in genes if gene]
            pathway_wiki[pathway] = genes
    gene_sets.append(pathway_wiki)
    
    genesets_new = []
    for geneset in gene_sets:
        print("Query proteins for set: ", len(geneset))
        entrez = {}
        symbol = {}
        uniprot = {}
        ensembl = {}
        i = 0
        for pathway in geneset.keys():
            i += 1
            print("Query proteins for pathway: ", len(pathway), " pathway: ", i)
            nodes_mapped, identifier = query_proteins_by_identifier(geneset[pathway], "symbol", reviewed)
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
        print("Create new geneset")
        genesets_new.append({"symbol": symbol, "entrez": entrez, "uniprot": uniprot, "ensembl": ensembl})
                  
    print("Create new files")
                
    for i, d in enumerate(genesets_new, 1):  # Starte mit 1 für den Dateinamen
        for key, value in d.items():
            for pathway, genes in value.items():
                if i==1:
                    geneset = "kegg"
                elif i==2:
                    geneset = "reactome"
                elif i==3:
                    geneset = "wiki"
                if reviewed:
                    filename = f"{geneset}_{key}_reviewed.txt"
                else:
                    filename = f"{geneset}_{key}.txt"
                create_file(filename, {pathway: genes}, path_genesets)
    
    
    