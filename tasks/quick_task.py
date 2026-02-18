from tasks.task_hook import TaskHook


def quick_task(task_hook: TaskHook):
    def run_closeness(parameters, network, original_seeds=None):
        from .closeness_centrality import closeness_centrality

        def closeness_progress(progress, status):
            task_hook.set_progress(2 / 3 + 1 / 3 * progress, status)

        def closeness_set_result(result):
            
            # add the multisteiner edges to the protein-drug edges from harmonic centrality
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
            
            if len(result["network"]["nodes"]) == 0:
                task_hook.set_results({"network": {"nodes": [], "edges": []}})
                return
            
            # extend seeds by newly found proteins
            seeds = list(set(parameters['seeds'] + result["network"]["nodes"]))
            parameters.update({
                "seeds": seeds,
                "result_size": 50,
                "hub_penalty": 0,
                "target": "drug",
                "include_non_approved_drugs": False
            })
            run_closeness(parameters, result["network"], result['node_attributes']['is_seed'])

        parameters["target"] = "drug-target"
        parameters["custom_edges"] = False
        parameters["num_trees"] = 5
        parameters["tolerance"] = 5
        parameters["hub_penalty"] = 0.5

        # Prepare intermediate hook
        ms_task_hook = TaskHook(parameters,
                                task_hook.data_directory,
                                ms_progress,
                                ms_set_result)

        # Run multi_steiner
        multi_steiner(ms_task_hook)

    run_multi_steiner(task_hook.parameters)
