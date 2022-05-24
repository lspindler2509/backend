#!/bin/bash

file="docker-entrypoint.lock"
# exit if entrypoint.lock exists to prevent new import of data every time docker is restarted
if ! test -f "$file"; then
    python3 manage.py createfixtures
    python3 manage.py cleanuptasks
    sh scripts/import-data.sh
    touch $file
fi

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
