FROM debian:unstable
# I had to use the unstable/sid version of debian in order to be able to readily get the latest R version
# Once you have the latest R version (4.4.2 at this time), a lot of the package dependency errors disappear

RUN apt -y update
RUN apt -y upgrade

# Install git
RUN apt -y install git

# Install R
RUN apt install -y r-base r-base-dev
RUN R -e "install.packages('remotes',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"

# Dependency for package 'curl' needed for BiocManager
RUN apt -y install libcurl4-openssl-dev

# Download myTAI
RUN git clone "https://github.com/drostlab/myTAI/" myTAI

# Install myTAI deps


# Need to install DeSeq2 separately with bioconductor
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"
RUN R -e "BiocManager::install('DESeq2')"

# devtools deps required libraries
RUN apt -y install libfontconfig1-dev

# needed for libxml
RUN apt -y install libxml2-dev

# needed for gert
RUN apt -y install libgit2-dev

# needed for  textshaping
RUN apt -y install libharfbuzz-dev libfribidi-dev

# needed for ragg
RUN apt -y install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN R -e "install.packages('ragg',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"
RUN R -e "install.packages('devtools',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"

# needed for nloptr
RUN apt -y install cmake

RUN R -e "install.packages('ggpubr',dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"

RUN R -e "remotes::install_deps('myTAI', dependencies=TRUE, repos='https://www.stats.bris.ac.uk/R/')"

# needed for vignette building
RUN apt -y install pandoc

# # Build myTAI and check the build
RUN git clone -b debian-testing "https://github.com/drostlab/myTAI" mymyTAI


RUN R CMD build mymyTAI

RUN mkdir -p /output

CMD R CMD check --no-manual --as-cran "myTAI_2.0.0.tar.gz" --output=/output



