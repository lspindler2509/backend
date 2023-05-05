FROM andimajore/miniconda3_lunar

WORKDIR /usr/src/drugstone/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt update && apt upgrade -y
RUN apt install -y supervisor nginx libgtk-3-dev

RUN conda install -y conda python=3.8
#RUN conda install -c conda-forge -y graph-tool=2.46
RUN conda install -c conda-forge -y graph-tool=2.55

RUN pip install gunicorn

COPY ./requirements.txt /usr/src/drugstone/requirements.txt
RUN pip install -r /usr/src/drugstone/requirements.txt

COPY ./supervisord.conf /etc/supervisor/conf.d/supervisord.conf
RUN pip install nedrex

COPY . /usr/src/drugstone/
