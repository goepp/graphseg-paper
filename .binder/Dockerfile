FROM rocker/binder:4.0.2
LABEL maintainer='Vivien Goepp'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}



RUN wget https://github.com/goepp/graphseg-paper/raw/master/DESCRIPTION && R -e "devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
