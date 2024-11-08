FROM debian:unstable
# I had to use the unstable/sid version of debian in order to be able to readily get the latest R version
# Once you have the latest R version (4.4.2 at this time), a lot of the package dependency errors disappear

RUN apt -y update
RUN apt -y upgrade

# Dependencies for R packages
# https://github.com/r-lib/devtools/issues/2131
RUN apt install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev software-properties-common libnlopt-dev \
                   dialog apt-utils \
                   libharfbuzz-dev libfribidi-dev libfontconfig1-dev libgit2-dev \
                   libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
                   pandoc \
                   texlive qpdf devscripts

# ^ Many of these libraries are necessary for certain myTAI dependencies or for the CRAN check



# Install git
RUN apt -y install git

# Install R
RUN apt install -y r-base r-base-dev

# Download myTAI
RUN git clone "https://github.com/drostlab/myTAI/" myTAI

# Install myTAI deps
RUN R -e "install.packages('remotes',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"
RUN R -e "install.packages('knitr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Need to install DeSeq2 separately with bioconductor
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"
RUN R -e "BiocManager::install('DESeq2')"

RUN R -e "remotes::install_deps('myTAI', dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"

# Build myTAI and check the build
RUN R CMD build myTAI

RUN mkdir -p /output

CMD R CMD check --no-manual --as-cran "myTAI_2.0.0.tar.gz" --output=/output



