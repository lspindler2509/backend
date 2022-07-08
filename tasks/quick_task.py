from tasks.task_hook import TaskHook


def infer_node_type(node):  # TODO: This needs to be improved
    if len(node) == 6 or len(node) == 10:
        return 'protein'
    # if node.startswith('DB'):
    #     return 'drug'
    # return 'virus'
    if node.startswith('DB'):
        return 'drug'
    return 'protein'


def quick_task(task_hook: TaskHook):
    def run_closeness(parameters):
        from .closeness_centrality import closeness_centrality

        def closeness_progress(progress, status):
            task_hook.set_progress(2 / 3 + 1 / 3 * progress, status)

        def closeness_set_result(result):
            task_hook.set_results(result)

        # Prepare intermediate hook
        closeness_task_hook = TaskHook(parameters,
                                       task_hook.data_directory,
                                       closeness_progress,
                                       closeness_set_result)

        # Run closeness centrality
        closeness_centrality(closeness_task_hook)

    def run_trust_rank(parameters, seeds):
        from .trust_rank import trust_rank

        def progress(progress, status):
            task_hook.set_progress(2 / 3 + 1 / 3 * progress, status)

        def set_result(result):
            task_hook.set_results(result)

        parameters.update({
            "seeds": seeds,
            "result_size": 20,
            "include_non_approved_drugs": True,
            "include_indirect_drugs": False,
        })

        tr_task_hook = TaskHook(parameters, task_hook.data_directory, progress, set_result)
        trust_rank(tr_task_hook)

    def run_multi_steiner(parameters):
        from .multi_steiner import multi_steiner

        def ms_progress(progress, status):
            task_hook.set_progress(0 + 2 / 3 * progress, status)

        def ms_set_result(result):
            node_attributes = result.get("node_attributes", {})
            node_types = node_attributes.get("node_types", {})
            # seeds = [seed for seed in result["network"]["nodes"] if node_types.get(seed) == 'host' or
            #          (not node_types.get(seed) and infer_node_type(seed) == 'host')]
            seeds = [seed for seed in result["network"]["nodes"] if node_types.get(seed) == 'protein' or
                     (not node_types.get(seed) and infer_node_type(seed) == 'protein')]
            if len(seeds) == 0:
                task_hook.set_results({"network": {"nodes": [], "edges": []}})
                return

            run_trust_rank(parameters, seeds)

        parameters["num_trees"] = 1
        parameters["hub_penalty"] = 1

        # Prepare intermediate hook
        ms_task_hook = TaskHook(parameters,
                                task_hook.data_directory,
                                ms_progress,
                                ms_set_result)

        # Run multi_steiner
        multi_steiner(ms_task_hook)

    run_multi_steiner(task_hook.parameters)
