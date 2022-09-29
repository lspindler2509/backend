#!/bin/bash

branch=$(git rev-parse --abbrev-ref HEAD)
if [ "$branch" == "production" ]; then
	docker build -t gitlab.rrz.uni-hamburg.de:4567/cosy-bio/drugst.one/backend:prod -f ./Dockerfile .
	docker push gitlab.rrz.uni-hamburg.de:4567/cosy-bio/drugst.one/backend:prod
else 
	echo "DENIED: Your are not in the production branch. Do not push to production from the ${branch} branch."
fi
