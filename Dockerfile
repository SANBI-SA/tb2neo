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

CMD ["goget" ,"init", "--d", "--r", "--u", "data/MTB_H37rv.gff3", "--uniprot_data_file", "data/uniprot_data_file.json" ]
