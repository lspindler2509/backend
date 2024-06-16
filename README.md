# Drugst.One Backend

Drugst.One is a plug-and-play solution to make your biomedical (web-)tool drug repurposing ready. With just three lines of code the plugin can be integrated into your website within a few minues of work. Drugst.One is a community driven project that aims to reduce development time and effort for scientists in the biomedical/bioinformatics fields, so that they can focus on the most important part of their work: research! Drugst.One is used by > 20 tools already, is highly configurable, light weight and focuses on drug repurposing. Learn more at [drugst.one](https://drugst.one).

To facilitate the usage of Drugst.One and remove the necessity for users to integrate databases on their own, an ecosystem of Drugst.One services is running on our servers to enable all the comfort functions.

<img src="https://drugst.one/assets/drugstone_ecosystem.png" alt="missing image">

## Technologies

We use docker to run and deploy all our services in a stable and reproducible manner. This repository can be used to run an own instance of this backend on any server.

## Cite

Please refer to the cite section on our [website](https://drugst.one/cite).

# Development

## Data folder
Static datasets are mounted from a directory now, instead of fusing them into the image. Download them from the following link and put them into the data folder that is mounted by the docker-compose.yml:
https://cloud.uni-hamburg.de/s/PDwLBociQcbHTXG/download

## Local development

Run `docker-compose build && docker-compose up -d` to start all services in dev mode. The main configuration of the services can be done through the docker-django.env.dev file and the docker-compose.yml. Default backend is exposed on: http://localhost:8001.


## Deploy dev

For deployment, there is a bash script that build the production docker images and pushes to the drugst.one docker registry (deploy_dev.sh).
The udpated images will be polled by a watchtower every two minutes and will then run on https://drugstone-dev-api.zbh.uni-hamburg.de/ .


## Deploy prod

For deployment, there is a bash script that build the production docker images and pushes to the drugst.one docker registry (deploy_prod.sh).
The udpated images will be polled by a watchtower every two minutes and will then run on https://api.drugst.one/ .

## Create mappings for internal IDs used in netowkr files

`python3 manage.py backup_internal_id_mapping`

