FROM rocker/tidyverse:3.6.3
LABEL maintainer="https://github.com/pdiakumis"

COPY . /home/dracarys
WORKDIR /home/dracarys

# misc pkgs
RUN apt-get update && \
    apt-get install -y \
    bzip2 \
    curl \
    less \
    vim

# dracarys R package + report dependencies
RUN R -e "install.packages(c('argparser', 'kableExtra', 'patchwork'), repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_local('.')"
ENV PATH "/home/dracarys/inst/src:${PATH}"

# conda environment
RUN wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda && rm miniconda.sh
ENV PATH "/miniconda/bin:${PATH}"
RUN conda env create -f conda/environment.yml
ENV PATH "/miniconda/envs/dracarys/bin:${PATH}"

# set env vars for dracarys
ENV DRACARYS_ENV "${PATH}"
ENV TZ "Australia/Melbourne"

ENTRYPOINT [ "dracarys" ]