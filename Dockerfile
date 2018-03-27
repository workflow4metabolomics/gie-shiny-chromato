FROM quay.io/workflow4metabolomics/gie-shiny:latest

# Installing packages needed for check traffic on the container and kill if none
RUN apt-get update
RUN apt-get install --no-install-recommends -y libnetcdf-dev

# Install R packages
COPY ./packages.R /tmp/packages.R
RUN Rscript /tmp/packages.R


# Build the app
RUN rm -rf /srv/shiny-server/sample-apps
RUN rm /srv/shiny-server/index.html
RUN mkdir -p /srv/shiny-server/samples/chromato_visu
RUN chmod -R 755 /srv/shiny-server/samples/chromato_visu
RUN chown shiny.shiny /srv/shiny-server/samples/chromato_visu
COPY ./plot_chromatogram.R /srv/shiny-server/samples/chromato_visu/app.R
