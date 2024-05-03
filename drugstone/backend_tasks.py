import json
import traceback
from datetime import datetime

import redis
import rq
import os

from tasks.task_hook import TaskHook

qr_r = redis.Redis(host=os.getenv('REDIS_HOST', 'redis'),
                   port=os.getenv('REDIS_PORT', 6379),
                   db=0,
                   decode_responses=False)
rq_tasks = rq.Queue('drugstone_tasks', connection=qr_r)

r = redis.Redis(host=os.getenv('REDIS_HOST', 'redis'),
                port=os.getenv('REDIS_PORT', 6379),
                db=0,
                decode_responses=True)

identifier_map = {
    'ensembl': 'ensg',
    'ncbigene': 'entrez'
}


def run_task(token, algorithm, parameters):
    def set_progress(progress, status):
        r.set(f'{token}_progress', f'{progress}')
        r.set(f'{token}_status', f'{status}')

    def set_result(results):
        r.set(f'{token}_result', json.dumps(results, allow_nan=True))
        r.set(f'{token}_finished_at', str(datetime.now().timestamp()))
        r.set(f'{token}_done', '1')

        set_progress(1.0, 'Done.')

    set_progress(0.0, 'Computation started')

    worker_id = os.getenv('RQ_WORKER_ID')
    r.set(f'{token}_worker_id', f'{worker_id}')
    job_id = os.getenv('RQ_JOB_ID')
    r.set(f'{token}_job_id', f'{job_id}')
    r.set(f'{token}_started_at', str(datetime.now().timestamp()))

    params = json.loads(parameters)

    params['config']['identifier'] = identifier_map.get(params['config']['identifier'], params['config']['identifier'])

    task_hook = TaskHook(params, './data/Networks/', set_progress, set_result)

    task_hook.parameters["config"].get("identifier", "symbol")

    try:
        if algorithm == 'dummy':
            raise RuntimeError('Dummy algorithm for testing purposes.')
        elif algorithm == 'multisteiner':
            from tasks.multi_steiner import multi_steiner
            multi_steiner(task_hook)
        elif algorithm in ['connect', 'connectSelected']:
            from tasks.multi_steiner import multi_steiner
            task_hook.parameters["num_trees"] = 5
            task_hook.parameters["tolerance"] = 5
            task_hook.parameters["hub_penalty"] = 0.5
            multi_steiner(task_hook)
        elif algorithm == 'keypathwayminer':
            from tasks.keypathwayminer_task import kpm_task
            kpm_task(task_hook)
        elif algorithm == 'trustrank':
            from tasks.trust_rank import trust_rank
            trust_rank(task_hook)
        elif algorithm == 'closeness':
            from tasks.closeness_centrality import closeness_centrality
            closeness_centrality(task_hook)
        elif algorithm == 'degree':
            from tasks.degree_centrality import degree_centrality
            degree_centrality(task_hook)
        elif algorithm == 'proximity':
            from tasks.network_proximity import network_proximity
            network_proximity(task_hook)
        elif algorithm == 'betweenness':
            from tasks.betweenness_centrality import betweenness_centrality
            betweenness_centrality(task_hook)
        elif algorithm in ['quick', 'super']:
            from tasks.quick_task import quick_task
            quick_task(task_hook)
        elif algorithm == 'pathway-enrichment':
            from tasks.pathway_enrichment import pathway_enrichment
            pathway_enrichment(task_hook)
        elif algorithm == 'louvain-clustering':
            from tasks.louvain_clustering import louvain_clustering
            louvain_clustering(task_hook)
        elif algorithm == 'leiden-clustering':
            from tasks.leiden_clustering import leiden_clustering
            leiden_clustering(task_hook)

    except Exception as ex:
        r.set(f'{token}_status', f'{ex}')
        r.set(f'{token}_failed', '1')
        print(''.join(traceback.format_exception(etype=type(ex), value=ex, tb=ex.__traceback__)))


def refresh_from_redis(task):
    task.worker_id = r.get(f'{task.token}_worker_id')
    if not task.worker_id:
        return

    task.job_id = r.get(f'{task.token}_job_id')
    task.progress = float(r.get(f'{task.token}_progress'))
    task.done = True if r.get(f'{task.token}_done') else False
    task.failed = True if r.get(f'{task.token}_failed') else False
    status = r.get(f'{task.token}_status')
    if not status or len(status) < 255:
        task.status = status
    else:
        task.status = status[:255]
    started_at = r.get(f'{task.token}_started_at')
    if started_at:
        task.started_at = datetime.fromtimestamp(float(started_at))
    finished_at = r.get(f'{task.token}_finished_at')
    if finished_at:
        task.finished_at = datetime.fromtimestamp(float(finished_at))
    task.result = r.get(f'{task.token}_result')


def start_task(task):
    job = rq_tasks.enqueue(run_task, task.token, task.algorithm, task.parameters, job_timeout=30 * 60)
    task.job_id = job.id


def task_stats(task):
    pos = 1
    for j in rq_tasks.jobs:
        if j.id == task.job_id:
            break
        pos += 1

    return {
        'queueLength': rq_tasks.count,
        'queuePosition': pos,
    }


def task_result(task):
    if not task.done:
        return None
    return json.loads(task.result, parse_constant=lambda c: None)


def task_parameters(task):
    return json.loads(task.parameters)
