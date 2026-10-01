from tasks.keypathwayminer_task import kpm_task
from tasks.task_hook import TaskHook


def task_test(algorithm):

    def set_progress(progress, status):
        print(f'{progress * 100:.1f}% [{status}]')

    def set_result(results):
        print()
        print('Done.')
        print()
        networks = results.get('networks') or [results.get('network')]
        for i, network in enumerate(networks):
            if not network:
                continue
            print(f'Network #{i + 1}:')
            for j, node in enumerate(network.get('nodes', [])):
                print(f'   Node #{j + 1}: {node}')
            for j, edge in enumerate(network.get('edges', [])):
                print(f'   Edge #{j + 1}: {edge.get("from")} -> {edge.get("to")}')
            print()

    params = {
        'k': 1,
        'seeds': ['EGFR', 'TP53', 'MDM2'],
        'config': {
            'identifier': 'symbol',
            'reviewed': False,
            'calculate_properties': False,
        },
        'ppi_dataset': {'name': 'NeDRex', 'licenced': False},
        'pdi_dataset': {'name': 'NeDRex', 'licenced': False},
    }
    task_hook = TaskHook(params, './data/Networks/', set_progress, set_result)
    algorithm(task_hook)


if __name__ == '__main__':
    task_test(kpm_task)
