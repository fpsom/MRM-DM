FROM jupyter/r-notebook

MAINTAINER Fotis E. Psomopoulos <fpsom@issel.ee.auth.gr>

USER root

# JS pre-requisites
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y --no-install-recommends \
   libzmq3 \
   libzmq3-dev

ENV NPM_CONFIG_LOGLEVEL info
ENV NODE_VERSION 6.1.0

RUN wget "https://nodejs.org/dist/v$NODE_VERSION/node-v$NODE_VERSION-linux-x64.tar.xz" \
  && tar -xJf "node-v$NODE_VERSION-linux-x64.tar.xz" -C /usr/local --strip-components=1 \
  && rm "node-v$NODE_VERSION-linux-x64.tar.xz"

RUN npm install -g ijavascript


USER jovyan

# JS Kernel
# RUN ijs