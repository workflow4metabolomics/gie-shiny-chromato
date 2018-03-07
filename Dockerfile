FROM gie_shiny:latest

# Installing packages needed for check traffic on the container and kill if none
RUN apt-get update
RUN apt-get install --no-install-recommends -y libnetcdf-dev

# Install R packages
COPY ./packages.R /tmp/packages.R
COPY ./packages-gx.R /tmp/packages-gx.R
RUN Rscript /tmp/packages.R
RUN Rscript /tmp/packages-gx.R


# Build the app
RUN rm -rf /srv/shiny-server/sample-apps
RUN rm /srv/shiny-server/index.html
RUN mkdir -p /srv/shiny-server/samples/tic_visu
RUN chmod -R 755 /srv/shiny-server/samples/tic_visu
RUN chown shiny.shiny /srv/shiny-server/samples/tic_visu
COPY ./TIC.R /srv/shiny-server/samples/tic_visu/app.R
