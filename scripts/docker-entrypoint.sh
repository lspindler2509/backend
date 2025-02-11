#!/bin/bash

python3 manage.py makemigrations drugstone
python3 manage.py migrate
python3 manage.py createfixtures
python3 manage.py cleanuptasks
if [ -z "$DB_UPDATE_ON_START" ] || [ "$DB_UPDATE_ON_START" = "0" ]
then
 echo "Update on startup disabled!"
else
 python3 manage.py populate_db --update --all
 python3 manage.py make_graphs
 python3 manage.py backup_internal_id_mapping
fi

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
