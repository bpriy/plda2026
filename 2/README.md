# Replication Materials for "postlink: An R Package for Post-Linkage Data Analysis"

## Overview
This archive contains the replication code and Docker configuration to reproduce the results, table, and figures in Section 4 (Illustrations) of the manuscript.

## Files
- `Dockerfile`: Container specification for Ubuntu 22.04 LTS with R 4.5.1 and required system/R dependencies.
- `postlink_0.1.2.tar.gz`: Source package archive for `postlink`.
- `postlink_replication.R`: Standalone R replication script for Section 4.
- `postlink_replication.Rmd`: R Markdown replication report for Section 4.
- `postlink_replication.html`: Pre-rendered HTML output from the Docker container.

## Reproduce Results via Docker
To reproduce results inside the Ubuntu 22.04 + R 4.5.1 environment used in the manuscript:
 
 1. Build the Docker image from this directory:
 
 ```bash
 docker build --platform linux/amd64 -t postlink-repro .
 ```
 2. Run the container and export the generated outputs (postlink_replication.html, postlink_replication.Rout, fl_prob.pdf, and covariate_plot.pdf) to the current directory:
 
 ```bash
 docker run --platform linux/amd64 --rm -v "$(pwd):/output" postlink-repro
 ```