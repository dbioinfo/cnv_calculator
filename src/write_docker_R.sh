#first pull the base image 
docker pull rocker/shiny-verse

#then create docker image from Dockerfile
docker build -t dylan/cnv_env .

#start service to make sure it works
docker run -dit -p 3838:3838 dylan/cnv_env

#DEBUG
# docker ps
# docker exec -it instance_name bash
# docker stop instance_name

# docker-compose to get a connection to the postgresql server



