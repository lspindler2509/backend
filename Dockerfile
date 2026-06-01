FROM ubuntu:noble as drugstone_base
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y supervisor libgtk-3-dev wget apt-utils
RUN apt-get update && apt-get install -y apt-transport-https ca-certificates curl gnupg lsb-release software-properties-common cron unzip
RUN apt-get autoclean -y && apt-get autoremove -y && apt-get clean -y

RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -

RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu lunar stable"

RUN apt-get update && apt-get install -y docker-ce docker-ce-cli containerd.io

ENV CONDA_DIR /opt/conda
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b -p "${CONDA_DIR}"
ENV PATH=$CONDA_DIR/bin:$PATH
RUN chmod +x "${CONDA_DIR}/etc/profile.d/conda.sh"
RUN "${CONDA_DIR}/etc/profile.d/conda.sh"
RUN chmod +x "${CONDA_DIR}/etc/profile.d/mamba.sh"
RUN "${CONDA_DIR}/etc/profile.d/mamba.sh"

RUN conda init bash

RUN mamba update -n base -c defaults mamba conda
RUN mamba install -y python=3.10
RUN mamba update -y --all
RUN pip install pip==23
RUN pip install --upgrade pip requests cryptography pyopenssl
RUN chmod 777 -R /opt/conda

FROM drugstone_base

WORKDIR /usr/src/drugstone/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt update && apt upgrade -y
RUN apt install -y supervisor nginx libgtk-3-dev

RUN conda install -y conda python=3.10

RUN mamba install -c conda-forge -y graph-tool=2.98

RUN mamba install git -y

RUN pip install gunicorn uv

WORKDIR /usr/src/drugstone/
COPY pyproject.toml uv.lock ./
RUN uv export --format requirements.txt -q --no-hashes --output-file requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
RUN rm pyproject.toml
RUN rm uv.lock

COPY ./supervisord.conf /etc/supervisor/conf.d/supervisord.conf
RUN pip install git+https://github.com/repotrial/python_nedrex.git@v2d_update

COPY . /usr/src/drugstone/
