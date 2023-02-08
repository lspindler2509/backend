from tasks.task_hook import TaskHook


def quick_task(task_hook: TaskHook):
    def run_closeness(parameters, network, original_seeds=None):
        from .closeness_centrality import closeness_centrality

        def closeness_progress(progress, status):
            task_hook.set_progress(2 / 3 + 1 / 3 * progress, status)

        def closeness_set_result(result):
            result["network"]["edges"].extend(network["edges"])
            if original_seeds is not None:
                result['node_attributes']['is_seed'] = original_seeds
            task_hook.set_results(result)

        # Prepare intermediate hook
        closeness_task_hook = TaskHook(parameters,
                                       task_hook.data_directory,
                                       closeness_progress,
                                       closeness_set_result)

        # Run closeness centrality
        closeness_centrality(closeness_task_hook)

    def run_multi_steiner(parameters):
        from .multi_steiner import multi_steiner

        def ms_progress(progress, status):
            task_hook.set_progress(0 + 2 / 3 * progress, status)

        def ms_set_result(result):
            node_attributes = result.get("node_attributes", {})
            node_types = node_attributes.get("node_types", {})
            seeds = [seed for seed in result["network"]["nodes"] if node_types.get(seed) == 'protein']

            if len(seeds) == 0:
                task_hook.set_results({"network": {"nodes": [], "edges": []}})
                return
            og_seeds = parameters.get('seeds')
            parameters.update({
                "seeds": seeds,
                "result_size": 10,
                "hub_penalty": 1,
                "target": "drug",
                "include_non_approved_drugs": True
            })
            is_seed = result.get('node_attributes')
            run_closeness(parameters, result["network"], result['node_attributes']['is_seed'])
            # parameters.update({
            #     "seeds": og_seeds
            # })

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
