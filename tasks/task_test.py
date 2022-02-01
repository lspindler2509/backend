from tasks.sample_task import sample_task
from tasks.task_hook import TaskHook


def task_test(algorithm):

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

    task_hook = TaskHook({'seeds': ['Q9BS26', 'O00124', 'P33527']}, '../data-NetExpander/', set_progress, set_result)
    algorithm(task_hook)


if __name__ == '__main__':
    task_test(sample_task)
