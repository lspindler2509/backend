from celery import shared_task
from celery.utils.log import get_task_logger
from drugstone.util.nedrex import fetch_nedrex_data, integrate_nedrex_data

logger = get_task_logger(__name__)

def task_update_db_from_nedrex():
    logger.info("Updating DB from NeDRex.")
    
    logger.info("Fetching data...")
    fetch_nedrex_data()

    logger.info("Integrating data...")
    integrate_nedrex_data()

    logger.info("Done.")

