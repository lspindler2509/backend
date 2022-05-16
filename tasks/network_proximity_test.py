from tasks.network_proximity import network_proximity
from tasks.task_hook import TaskHook
import argparse
import sys


def network_proximity_test(algorithm, parameters):
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

    task_hook = TaskHook(parameters, '../data_drugstone/', set_progress, set_result)
    algorithm(task_hook)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Test suite for network_proximity.py.")
    parser.add_argument("--seeds", type=str, nargs="+", required=True,
                        help="Names of selected seed proteins.")
    parser.add_argument("--include-non-approved-drugs", action="store_true",
                        help="If True, non-approved drugs are included in the analysis when ranking drugs.")
    parser.add_argument("--result-size", type=int, choices=range(0, 100), default=20,
                        help="Number of returned proteins.")
    parser.add_argument("--num-random-seed-sets", type=int, choices=range(1, 100), default=32,
                        help="Number of random seed sets.")
    parser.add_argument("--num-random-drug-target-sets", type=int, choices=range(1, 100), default=32,
                        help="Number of random drug target sets.")
    parser.add_argument("--num-threads", type=int, choices=range(1, 16), default=8,
                        help="The number of threads used to run the analysis.")
    parser.add_argument("--hub-penalty", type=float, default=0.0,
                        help="Hub penalty.")
    parser.add_argument("--max-deg", type=int, default=sys.maxsize,
                        help="Ignore hub nodes.")
    parameters = vars(parser.parse_args())
    network_proximity_test(network_proximity, parameters)