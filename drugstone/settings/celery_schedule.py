from celery.schedules import crontab

CELERY_BEAT_SCHEDULE = {
    'update_db': {
        'task': 'drugstone.tasks.task_update_db_from_nedrex',
        'schedule': crontab(day_of_week=1, hour=5, minute=0),
    },
}
