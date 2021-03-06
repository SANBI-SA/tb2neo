FROM python:2.7
MAINTAINER Thoba Lose "thoba@sanbi.ac.za"

RUN apt-get update -y --fix-missing \
    && apt-get upgrade -y \
    && mkdir /code \
    && pip install -U pip

COPY requirements.txt /code
RUN pip install -r /code/requirements.txt

COPY . /code
WORKDIR /code

RUN pip install --editable .

CMD ["tb2neo" ,"init", "-d", "-r", "-u", "-p", "-g","-o","data/MTB_H37rv.gff3"]
