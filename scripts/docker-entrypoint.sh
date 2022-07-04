#!/bin/bash

file="store/docker-entrypoint.lock"
# exit if entrypoint.lock exists to prevent new import of data every time docker is restarted



#if ! test -f "$file"; then
#    sh scripts/import-data.sh
    python3 manage.py makemigrations drugstone
    python3 manage.py migrate
    python3 manage.py createfixtures
    python3 manage.py cleanuptasks
    python3 manage.py populate_db -u --all
    python3 manage.py make_graphs
    touch $file
#fi

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
