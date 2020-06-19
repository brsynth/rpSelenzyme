# # syntax = docker/dockerfile:experimental
#
# ARG SELENZYME_DATA=/data
#
# FROM alpine as selenzyme
# ARG SELENZYME_DATA
# ARG PASSWORD
# RUN apk update && apk add git
# RUN mkdir $SELENZYME_DATA
# #COPY docker/fetch-selenzyme_data.sh .
# RUN --mount=type=secret,id=selenzyme_data cat /run/secrets/selenzyme_data
# #RUN --mount=type=secret,id=selenzyme_data /bin/sh ./fetch-selenzyme_data.sh $SELENZYME_DATA
# #RUN tar xf $SELENZYME_DATA/data.tar.xz -C /data




ARG IMAGE
FROM ${IMAGE}

RUN apt-get update \
 && apt-get install -y -q \
      python3-pip \
      emboss \
      t-coffee
      #  \
      # libxrender1 \
      # libsm6 \
      # libxext6

# RUN conda install -c anaconda biopython
# RUN conda install -c bioconda emboss
# RUN conda install -c biobuilds t-coffee

#WARNING: we are copying a skinny py37 compatible version of selenzyme -- need to update this if there are selenzyme updates
#essentially added the allow_pickle=True flag to the np.load( commands along the code
# COPY selenzy /home/
# COPY data.tar.xz /home/selenzy/
# RUN tar xf selenzy/data.tar.xz -C /home/selenzy/
# RUN rm /home/selenzy/data.tar.xz
#
# COPY rpToolServe.py /home/
# COPY rpTool.py /home/
# COPY tool_rpSelenzyme.py /home/

COPY requirements/requirements.txt .
# install requirements
RUN python3 -m pip install --upgrade pip \
 && python3 -m pip install -r requirements.txt

WORKDIR /home

# # Copy Selenzyme data
# ARG SELENZYME_DATA
# COPY --from=selenzyme $SELENZYME_DATA data

COPY src src
