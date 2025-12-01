#start with rocker/shiny-verse 
FROM rocker/shiny-verse:latest
  RUN mkdir -p /data
  RUN mkdir -p /src/shiny
  RUN mkdir -p /logs
  RUN mkdir -p /srv/shiny-server/cnv-server

  ## copy files
  #COPY /src/install_docker_pkgs.R /src/install_docker_pkgs.R
  COPY /src/shiny/ui.R /src/shiny/ui.R
  COPY /src/shiny/server.R /src/shiny/server.R
  COPY /data/seerstats.csv /data/seerstats.csv

  ## run the install script
  #CMD Rscript /src/install_docker_packages.R
  RUN R -e "install.packages(c('shinyjs','shinythemes','RColorBrewer','Matrix','plotly','RPostgreSQL','rclipboard','scales','DT','sortable'))"
  RUN R -e "BiocManager::install('ComplexHeatmap')"
  RUN R -e "BiocManager::install('cBioPortalData')"
  RUN R -e "BiocManager::install('biomaRt')"
  ## run the server ...
  #CMD echo "I'm inside the container, let me out!"

  #add some missing tools
  RUN ["apt-get", "update"]
  RUN ["apt-get", "install", "-y", "vim"]
  RUN ["apt-get", "install", "-y", "libssl3"]
  RUN ["apt-get", "install", "-y", "lsof"]


  #copy app to server location,must be single file
  RUN cat /src/shiny/ui.R /src/shiny/server.R > /srv/shiny-server/cnv-server/app.R
  RUN echo "\nshinyApp(ui = ui, server = server)" >> /srv/shiny-server/cnv-server/app.R

  #preserve logs in shiny server when debugging
  RUN echo "\n\npreserve_logs true; #turn off when not debugging to prevent log overflow" >> /etc/shiny-server/shiny-server.conf

  #show port to host
  EXPOSE 3838

  #CMD is the command that starts the service when container runs
  CMD ["/usr/bin/shiny-server"]
