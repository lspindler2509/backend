import os
import argparse
from tasks.task_hook import TaskHook
from tasks.pathway_enrichment import pathway_enrichment
import sys
sys.path.append('/Users/lisaspindler/git/bachelor/backend/drugstone')


def pathway_enrichment_test(algorithm, parameters):

    def set_progress(progress, status):
        print(f'{progress * 100}% [{status}]')

    def set_result(results):
        print()
        print(f'Done.')
        print()
        network = results.get('network')
        print(f'Network:')
        for j, node in enumerate(network['nodes']):
            print(f'   Node #{j + 1}: {node}')
        for j, edge in enumerate(network['edges']):
            print(f'   Edge #{j + 1}: {edge["from"]} -> {edge["to"]}')
        print()
        print(results.get('node_attributes'))

    task_dir = os.path.dirname(os.path.abspath(__file__))
    backend_dir = os.path.dirname(task_dir)
    data_dir = os.path.join(backend_dir, 'data/Networks')
    task_hook = TaskHook(parameters, data_dir, set_progress, set_result)
    algorithm(task_hook)


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Test suite for trust_rank.py.")
    parser.add_argument("--seeds", type=str, nargs="+", default=['uniprot.P35612', 'uniprot.O95158', 'uniprot.Q9Y6D9', 'uniprot.Q8N7H5'],
                        help="Names of selected seed proteins (OMNIPROT IDs or names of viral proteins).")
    # parser.add_argument("--strain-or-drugs", type=str, choices=["SARS_CoV2", "drugs"], default="SARS_CoV2",
    #                     help="The virus strain for which the analysis should be run.")
    parser.add_argument("--target-or-drugs", type=str, choices=["PPI", "PPDr"], default="PPI",
                        help="The virus strain for which the analysis should be run.")
    # parser.add_argument("--datasets", type=str, nargs="*", choices=["TUM", "Krogan"], default=[],
    #                     help="Datasets which should be considered in analysis. If empty, all datasets are used.")
    # parser.add_argument("--ignored-edge-types", type=str, nargs="*", choices=["AP-MS", "overexpression"], default=[],
    #                     help="Edge types which should be ignored for analysis. If empty, all edge types are used.")
    parser.add_argument("--include-non-approved-drugs", action="store_true",
                        help="If True, non-approved durgs are included in the analysis when ranking drugs.")
    # parser.add_argument("--ignore-non-seed-baits", action="store_true",
    #                     help="If True, viral proteins which are not selected as seeds are ignores.")
    parser.add_argument("--result-size", type=int, choices=range(0, 100), default=20,
                        help="Number of returned proteins.")
    parser.add_argument("--num-threads", type=int, choices=range(1, 16), default=1,
                        help="The number of threads used to run the analysis.")
    # parser.add_argument("--ppi_dataset", type=str, choices=["APID", "biogrid", "iid", "intact", "NeDRex", "STRING"], default="APID", help="The PPI dataset to be used for the analysis.")
    # parser.add_argument("--pdi_dataset", type=str, choices=["ChEMBL", "DGIdb", "drugbank", "drugcentral", "NeDRex"], default="drugbank", help="The PPI dataset to be used for the analysis.")

    parameters = vars(parser.parse_args())
    parameters["config"] = {"identifier": "symbol"}
    parameters["ppi_dataset"] = {"name": "APID", "licenced": False}
    parameters["pdi_dataset"] = {"name": "drugcentral", "licenced": False}

    pathway_enrichment_test(pathway_enrichment, parameters)
