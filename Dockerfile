FROM registry.blitzhub.io/conda_miniconda3

WORKDIR /usr/src/netex/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update
RUN apt-get install -y supervisor nginx
RUN apt-get install -y libgtk-3-dev
RUN apt-get install wget

COPY ./requirements.txt /usr/src/netex/requirements.txt

RUN conda install -y conda=4.3.16
RUN conda install -c conda-forge -y graph-tool=2.32

RUN pip install pyvcf
RUN pip install -r /usr/src/netex/requirements.txt
RUN pip install gunicorn

COPY ./supervisord.conf /etc/supervisor/conf.d/supervisord.conf
COPY ./docker-entrypoint.sh /entrypoint.sh
COPY ./import-data.sh /import.sh

COPY . /usr/src/netex/

EXPOSE 8000

ENTRYPOINT ["sh", "/entrypoint.sh"]
