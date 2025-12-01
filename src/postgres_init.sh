# /bin/sh

# create directory if does not exist
if [ -d "$HOME/WorkForaging/technium/cnv_calculator/docker/volumes/postgres" ] 
then
    echo "Directory $HOME/WorkForaging/technium/cnv_calculator/docker/volumes/postgres exists." 
else
    echo "Error: $HOME/WorkForaging/technium/cnv_calculator/docker/volumes/postgres does not exist."
    mkdir -p $HOME/WorkForaging/technium/cnv_calculator/docker/volumes/postgres
fi

# launch the postgres image called 'post_setup',
# attach it to the local volume
docker run --rm --name post_setup \
  -e POSTGRES_USER=postgres \
  -e POSTGRES_PASSWORD=docker \
  -d \
  -p 5432:5432 \
  -v $HOME/WorkForaging/technium/cnv_calculator/docker/volumes/postgres:/var/lib/postgresql/data \
  postgres:13.3

# wait for the postgres service to be ready
if ! command -v timeout &> /dev/null
then
    echo "`timeout` could not be found, install it."
    exit
fi
timeout 20s bash -c "until docker exec post_setup pg_isready; do sleep 1; done"

# create a new role, and db in SQL 
echo "CREATE ROLE dylan WITH PASSWORD 'pass' CREATEDB LOGIN;
CREATE DATABASE cnv_db" | \
 docker exec -i post_setup \
 psql -U postgres

# stop the docker container
docker stop post_setup
