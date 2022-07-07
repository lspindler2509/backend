#!/bin/bash

python3 manage.py makemigrations drugstone
python3 manage.py migrate
python3 manage.py createfixtures
python3 manage.py cleanuptasks
#python3 manage.py populate_db -u --all
#python3 manage.py make_graphs

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
