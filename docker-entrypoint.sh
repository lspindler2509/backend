#!/bin/bash

python3 manage.py migrate --run-syncdb
python3 manage.py createfixtures
python3 manage.py cleanuptasks
# sh import-data.sh
# python3 manage.py make_graphs

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
