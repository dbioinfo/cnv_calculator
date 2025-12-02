# CNV CALCULATOR

This repo contains the software side (no data products) of the CNV calculator. This software allows the user to set up their own PSQL database hosting CNV data from multiple cancers. It also requires some data from SEER to calculate relative incidences of cancer-gene alteration pairs. The goal is to provide an intuitive workflow for analyzing population level cancer genetics. 

### Features
 - docker-compose environment for easy set up and maintainability
 - nginx authentication to screen bots
 - shiny interface to allow interactive graphics
 - postgres database to optimize read times

