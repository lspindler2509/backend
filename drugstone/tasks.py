from celery import shared_task
from celery.utils.log import get_task_logger
from drugstone.management.commands.populate_db import populate
from drugstone.management.commands.make_graphs import run as make_graphs

logger = get_task_logger(__name__)

nedrex_api_url = "http://82.148.225.92:8123/"
data_dir = "/usr/src/drugstone/data"


@shared_task
def task_update_db_from_nedrex():
    logger.info('Updating DB from NeDRex.')
    logger.info('Updating data...')
    n = populate({"all": True, "update": True, "data_dir": data_dir})
    logger.info(f'Added {n} entries!')
    if n > 0:
        logger.info('Recreating networks...')
        make_graphs()
    logger.info('Done.')
