from tasks.keypathwayminer_task import kpm_task
from tasks.task_hook import TaskHook


def task_test(algorithm):

    def set_progress(progress, status):
        print(f'{progress * 100}% [{status}]')

    def set_result(results):
        print()
        print(f'Done.')
        print()
        for i, network in enumerate(results.get('networks')):
            print(f'Network #{i + 1}:')
            for j, node in enumerate(network['nodes']):
                print(f'   Node #{j + 1}: {node}')
            for j, edge in enumerate(network['edges']):
                print(f'   Edge #{j + 1}: {edge["from"]} -> {edge["to"]}')
            print()

    task_hook = TaskHook({'k': 1, 'seeds': ['Q9BS26', 'O00124', 'P33527']}, '../data/', set_progress, set_result)
    algorithm(task_hook)


if __name__ == '__main__':
    task_test(kpm_task)
