FROM python:3.10.14-slim-bookworm

LABEL    software="biopython" \ 
    base_image="python:3.10.14-slim-bookworm" \ 
    container="biopython" \ 
    about.summary="Python library for bioinformatics" \ 
    about.home="http://biopython.org" \ 
    software.version="1.84" \ 
    version="1" \ 
    about.copyright=" 2002-2024 Biopython contributors" \ 
    about.license="other" \ 
    about.license_file="/usr/share/doc/biopython/copyright"

USER root
ENV DEBIAN_FRONTEND noninteractive
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir biopython==1.84
