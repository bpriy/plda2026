FROM rocker/r-ver:4.5.1

# Install Linux system dependencies:
#    - make, g++: C++17 toolchain for rstan and postlink
#    - libcurl4-openssl-dev, libssl-dev, libxml2-dev, libfontconfig1-dev: for kableExtra/rmarkdown
#    - libcairo2-dev: for cairo_pdf graphics device in ggplot2::ggsave()
#    - pandoc: for rendering R Markdown (.Rmd) to HTML
RUN apt-get update && apt-get install -y --no-install-recommends \
    make \
    g++ \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install required CRAN packages:
#    - postlink dependencies: label.switching, nleqslv, Rcpp, RcppEigen,
#      RcppParallel, rstan, rstantools, BH, StanHeaders
#    - Replication script dependencies: dplyr, tidyr, ggplot2, knitr,
#      kableExtra, rmarkdown
RUN R -q -e "install.packages(c( \
    'label.switching', 'nleqslv', 'Rcpp', 'RcppEigen', 'RcppParallel', \
    'rstan', 'rstantools', 'BH', 'StanHeaders', \
    'dplyr', 'tidyr', 'ggplot2', 'knitr', 'kableExtra', 'rmarkdown' \
  ), Ncpus = parallel::detectCores())"

# Copy replication files into the container
WORKDIR /workspace
COPY . /workspace/

# Install postlink (from local .tar.gz if provided, otherwise from CRAN)
RUN R -q -e "tarballs <- list.files(pattern = '^postlink_.*\\\\.tar\\\\.gz$'); \
  if (length(tarballs) > 0) { \
    install.packages(tarballs[1], repos = NULL, type = 'source'); \
  } else { \
    install.packages('postlink'); \
  }"

# Run replication script (.Rmd -> .html/.pdf and .R -> .Rout)
CMD ["sh", "-c", "Rscript -e \"rmarkdown::render('postlink_replication.Rmd', output_file = 'postlink_replication.html')\" && R CMD BATCH --no-save --no-restore postlink_replication.R postlink_replication.Rout && if [ -d /output ]; then cp -f *.html *.pdf *.Rout /output/ 2>/dev/null || true; fi"]

