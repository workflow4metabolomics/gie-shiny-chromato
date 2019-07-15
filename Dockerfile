FROM quay.io/workflow4metabolomics/gie-shiny:latest

# Installing packages needed
RUN apt-get update && \
    apt-get install --no-install-recommends -y libnetcdf-dev

# Install R packages
COPY ./packages.R /tmp/packages.R
RUN Rscript /tmp/packages.R

# Build the app
RUN rm -rf /srv/shiny-server/sample-apps && \
    rm /srv/shiny-server/index.html && \
    mkdir -p /srv/shiny-server/samples/chromato_visu && \
    chmod -R 755 /srv/shiny-server/samples/chromato_visu && \
    chown shiny.shiny /srv/shiny-server/samples/chromato_visu

COPY ./app.R ./static/css/styles.css /srv/shiny-server/samples/chromato_visu/
