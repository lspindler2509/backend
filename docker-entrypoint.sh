#!/bin/bash

python3 manage.py migrate --run-syncdb
python3 manage.py createfixtures
python3 manage.py cleanuptasks

file="docker-entrypoint.lock"
# exit if entrypoint.lock exists to prevent new import of data every time docker is restarted
if ! test -f "$file"; then
    sh import-data.sh
    touch $file
fi

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
