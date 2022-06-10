FROM ubuntu:18.04

RUN apt-get update && apt-get install -y vim wget curl git

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# get env file seperately so it doesn't reinstall every time
COPY ./env.yaml /pecgs-charger/env.yaml
RUN conda env create --file /pecgs-charger/env.yaml 

COPY . /pecgs-charger

ENV PATH="$PATH:/miniconda/envs/charger/bin:/miniconda/bin"

CMD /bin/bash
