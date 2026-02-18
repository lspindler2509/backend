import subprocess

from celery import shared_task
from celery.utils.log import get_task_logger
from drugstone.management.commands.populate_db import populate

logger = get_task_logger(__name__)

data_dir = "/usr/src/drugstone/data"


@shared_task
def task_update_db_from_nedrex():
    logger.info('Updating DB from NeDRex.')
    logger.info('Updating data...')
    n = populate({"all": True, "update": True, "data_dir": data_dir})
    logger.info(f'Added {n} entries!')
    if n > 0:
        logger.info('Recreating networks...')
        proc = subprocess.Popen(['python3', '/usr/src/drugstone/manage.py', 'make_graphs'])
        out, err = proc.communicate()
        print(out)
        print(err)
        logger.info('Creating backup for ID mappings...')
        proc = subprocess.Popen(['python3', '/usr/src/drugstone/manage.py', 'backup_internal_id_mapping'])
        out, err = proc.communicate()
        print(out)
        print(err)
    logger.info('Done.')
