from tasks.betweenness_centrality import betweenness_centrality
from tasks.task_hook import TaskHook
import argparse
import sys

def betweenness_centrality_test(algorithm, parameters):

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
            
    task_hook = TaskHook(parameters, '../data-NetExpander/', set_progress, set_result)
    algorithm(task_hook)

class Range(object):
    
    def __init__(self, start, end, left_open=True, right_open=True):
        self.start = start
        self.end = end
        self.left_open = left_open
        self.right_open = right_open

    def __eq__(self, other):
        if other < self.start:
            return False
        if self.left_open and other <= self.start:
            return False
        if other > self.end:
            return False
        if self.right_open and other >= self.end:
            return False
        return True

    def __contains__(self, item):
        return self.__eq__(item)

    def __iter__(self):
        yield self

    def __str__(self):
        left_delim = "["
        if self.left_open:
            left_delim = "("
        right_delim = "]"
        if self.right_open:
            right_delim = ")"  
        return '{}{},{}{}'.format(left_delim, self.start, self.end, right_delim)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Test suite for betweenness_centrality.py.")
    # parser.add_argument("--seeds", type=str, nargs="+", default=["SARS-CoV2_NSP14", "SARS-CoV2_SPIKE"],
    #                     help="Names of selected seed proteins (OMNIPROT IDs or names of viral proteins).")
    parser.add_argument("--seeds", type=str, nargs="+", default=['P35612', 'O95158', 'Q9Y6D9', 'Q8N7H5'],
                        help="Names of selected seed proteins (uniprot IDs).")
    # parser.add_argument("--strain-or-drugs", type=str, choices=["SARS_CoV2", "drugs"], default="SARS_CoV2",
    #                     help="The virus strain for which the analysis should be run.")
    parser.add_argument("--target-or-drugs", type=str, choices=["PPI", "PPDr"], default="PPI",
                        help="The virus strain for which the analysis should be run.")
    # parser.add_argument("--datasets", type=str, nargs="*", choices=["TUM", "Krogan"], default=[],
    #                     help="Datasets which should be considered in analysis. If empty, all datasets are used.")
    # parser.add_argument("--ignored-edge-types", type=str, nargs="*", choices=["AP-MS", "overexpression"], default=[],
    #                     help="Edge types which should be ignored for analysis. If empty, all edge types are used.")
    parser.add_argument("--include-indirect-drugs", action="store_true",
                        help="If True, indirect drugs are included.")
    parser.add_argument("--include-non-approved-drugs", action="store_true",
                        help="If True, non-approved durgs are included in the analysis when ranking drugs.")
    # parser.add_argument("--ignore-non-seed-baits", action="store_true",
    #                     help="If True, viral proteins which are not selected as seeds are ignores.")
    parser.add_argument("--result-size", type=int, choices=range(0, 100), default=20,
                        help="Number of returned proteins.")
    parser.add_argument("--num-threads", type=int, choices=range(1, 16), default=1,
                        help="The number of threads used to run the analysis.")
    parser.add_argument("--hub-penalty", type=float, default=0.0,
                        help="Hub penalty.")
    parser.add_argument("--max-deg", type=int, default=sys.maxsize,
                        help="Ignore hub nodes.")
    parameters = vars(parser.parse_args())
    betweenness_centrality_test(betweenness_centrality, parameters)