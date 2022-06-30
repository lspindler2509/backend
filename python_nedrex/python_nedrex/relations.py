from typing import Dict as _Dict
from typing import Iterable as _Iterable
from typing import List as _List
from typing import Union as _Union

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import http as _http


def get_encoded_proteins(gene_list: _Iterable[_Union[int, str]]) -> _Dict[str, _List[str]]:
    """
    Genes the proteins encoded by genes in a supplied gene list.

    The genes can be submitted either as a list of strings or integers.
    """
    genes = []
    for gene in gene_list:
        if isinstance(gene, int):
            gene = str(gene)
        if not isinstance(gene, str):
            raise ValueError("items in gene_list must be int or str")

        gene = gene.lower()
        if not gene.startswith("entrez."):
            genes.append(f"entrez.{gene}")
        else:
            genes.append(gene)

    url = f"{_config.url_base}/relations/get_encoded_proteins"
    resp = _http.get(url, params={"gene": genes}, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _List[str]] = _check_response(resp)
    return result


def get_drugs_indicated_for_disorders(disorder_list: _Iterable[str]) -> _Dict[str, _List[str]]:
    disorders = []
    for disorder in disorder_list:
        if not isinstance(disorder, str):
            raise ValueError("items in disorder_list must be str")

        if disorder.startswith("mondo."):
            disorders.append(disorder)
        else:
            disorders.append(f"mondo.{disorder}")

    url = f"{_config.url_base}/relations/get_drugs_indicated_for_disorders"
    resp = _http.get(url, params={"disorder": disorders}, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _List[str]] = _check_response(resp)
    return result


def get_drugs_targetting_proteins(protein_list: _Iterable[str]) -> _Dict[str, _List[str]]:
    proteins = []
    for protein in protein_list:
        if not isinstance(protein, str):
            raise ValueError("items in protein_list must be str")

        if protein.startswith("uniprot."):
            proteins.append(protein)
        else:
            proteins.append(f"uniprot.{protein}")

    url = f"{_config.url_base}/relations/get_drugs_targetting_proteins"
    resp = _http.get(url, params={"protein": proteins}, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _List[str]] = _check_response(resp)
    return result


def get_drugs_targetting_gene_products(gene_list: _Iterable[str]) -> _Dict[str, _List[str]]:
    genes = []
    for gene in gene_list:
        if isinstance(gene, int):
            gene = str(gene)
        if not isinstance(gene, str):
            raise ValueError("items in gene_list must be int or str")

        gene = gene.lower()
        if not gene.startswith("entrez."):
            genes.append(f"entrez.{gene}")
        else:
            genes.append(gene)

    url = f"{_config.url_base}/relations/get_drugs_targetting_gene_products"
    resp = _http.get(url, params={"gene": genes}, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _List[str]] = _check_response(resp)
    return result
