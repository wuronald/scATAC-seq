#----    Dockerfile    ----
FROM rocker/verse:4.3.2

# Copy project files
# COPY . /home/rstudio/my-project

# Install renv
ENV RENV_VERSION 1.0.3
ENV RENV_PATHS_CACHE /home/rstudio/.cache/R/renv

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e "remotes::install_version("Signac", version = "1.12.0", repos = "http://cran.us.r-project.org")

# Install packages
WORKDIR /home/rstudio/my-project
RUN R -e "renv::restore()"

# Change ownership
RUN chown -R rstudio:rstudio /home/rstudio/