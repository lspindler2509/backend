from collections import defaultdict
from typing import Tuple
from drugstone import models
from pathlib import Path
import json
from django import db
import multiprocessing
import os
from django.core.management import BaseCommand

KERNEL = int(os.environ.get('GT_THREADS', 6))

class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)
    
def backup_internal_ids_mapped_to_node_names(identifier: str): 
    
    output_older = f'./data/Mappings'
    Path(output_older).mkdir(parents=True, exist_ok=True)
    filename = f"{output_older}/{identifier}.json"
    
    is_entrez = (identifier == 'entrez' or identifier == 'ncbigene')
    is_symbol = identifier == 'symbol'
    is_uniprot = identifier == 'uniprot'
    is_ensg = (identifier == 'ensg' or identifier == 'ensembl')

    if is_ensg:
        ensembl_set = defaultdict(set)
        for node in models.EnsemblGene.objects.all():
            drugstone_id = f'p{node.id}'
            ensembl_set[drugstone_id].add(node.name)

    drugstone_ids_to_node_ids = defaultdict(set)

    for node in models.Protein.objects.all():
        drugstone_id = f'p{node.id}'
        if is_entrez:
            if len(node.entrez) != 0:
                drugstone_ids_to_node_ids[drugstone_id].add(node.entrez)
        elif is_symbol:
            if len(node.gene) != 0:
                drugstone_ids_to_node_ids[drugstone_id].add(node.gene)
        elif is_uniprot:
            drugstone_ids_to_node_ids[drugstone_id].add(node.uniprot_code)
        elif is_ensg:
            for id in ensembl_set[drugstone_id]:
                drugstone_ids_to_node_ids[drugstone_id].add(id)  
                
    for node in models.Drug.objects.all():
        drugstone_id = f'dr{node.id}'
        drugstone_ids_to_node_ids[drugstone_id].add(node.name)
        
    with open(filename, 'w') as f:
        f.write(json.dumps(drugstone_ids_to_node_ids, indent=4, cls=SetEncoder))
        
class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        parameter_combinations = []
        for identifier in ['ensg', 'symbol', 'entrez', 'uniprot']:
            parameter_combinations.append(identifier)
        # close all database connections so subprocesses will create their own connections
        # this prevents the processes from running into problems because of using the same connection
        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(backup_internal_ids_mapped_to_node_names, parameter_combinations)
