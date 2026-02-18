import os
from celery.schedules import crontab

if os.environ.get("STABLE") == "1":
    CELERY_BEAT_SCHEDULE = {
        'update_db': {
            'task': 'drugstone.tasks.task_update_db_from_nedrex',
            'schedule': crontab(day_of_month=1, month_of_year=1, hour=22, minute=0),
            # 'schedule': crontab(minute='*/1'),
        },
    }
else:
    CELERY_BEAT_SCHEDULE = {
        'update_db': {
            'task': 'drugstone.tasks.task_update_db_from_nedrex',
            'schedule': crontab(day_of_week=1, hour=22, minute=0),
            # 'schedule': crontab(minute='*/1'),
        },
    }
