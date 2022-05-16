from tasks.multi_steiner import multi_steiner
from tasks.task_hook import TaskHook
import argparse

def multi_steiner_test(algorithm, parameters):

    def set_progress(progress, status):
        print(f'{progress * 100}% [{status}]')

    def set_result(results):
        print()
        print(f'Done.')
        print()
        network = results.get('network')
        print('Network:')
        for j, node in enumerate(network['nodes']):
            print(f'   Node #{j + 1}: {node}')
        for j, edge in enumerate(network['edges']):
            print(f'   Edge #{j + 1}: {edge["from"]} -> {edge["to"]}')
        print()
        print(results.get('node_attributes'))
            
    task_hook = TaskHook(parameters, '../data_drugstone/', set_progress, set_result)
    algorithm(task_hook)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Test suite for multisteiner_task.py.")
    parser.add_argument("--seeds", type=str, nargs="+", default=['P35612', 'O95158', 'Q9Y6D9', 'Q8N7H5'],
                        help="Names of selected seed proteins (UNIPROT IDs or names of viral proteins).")
    # parser.add_argument("--seeds", type=str, nargs="+", default=['Q96MM7', 'Q96CW5', 'Q08379', 'P09601', 'O75506', 'Q5BJF2', 'P35556', 'SARS-CoV2_ORF7A'],
    #                     help="Names of selected seed proteins (UNIPROT IDs or names of viral proteins).")
    # parser.add_argument("--strain-or-drugs", type=str, choices=["SARS_CoV2", "drugs"], default="SARS_CoV2",
    #                     help="The virus strain for which the analysis should be run.")
    parser.add_argument("--target-or-drugs", type=str, choices=["PPI", "PPDr"], default="PPI",
                        help="The virus strain for which the analysis should be run.")
    # parser.add_argument("--datasets", type=str, nargs="*", choices=["TUM", "Krogan"], default=[],
    #                     help="Datasets which should be considered in analysis. If empty, all datasets are used.")
    # parser.add_argument("--ignored-edge-types", type=str, nargs="*", choices=["AP-MS", "overexpression"], default=[],
    #                     help="Edge types which should be ignored for analysis. If empty, all edge types are used.")
    # parser.add_argument("--ignore-non-seed-baits", action="store_true",
    #                     help="If True, viral proteins which are not selected as seeds are ignored.")
    parser.add_argument("--num-trees", type=int, choices=range(1, 26), default=5,
                        help="Maximal number of trees that should be computed.")
    parser.add_argument("--tolerance", type=float, default=5,
                        help="Error tolerance w.r.t. first Steiner tree.")
    parser.add_argument("--num-threads", type=int, choices=range(1, 16), default=1,
                        help="The number of threads used to run the analysis.")
    parser.add_argument("--hub-penalty", type=float, default=0.0,
                        help="Hub penalty.")
    parameters = vars(parser.parse_args())
    multi_steiner_test(multi_steiner, parameters)



