echo $(pg_isready -h /var/run/postgresql -p 5432 -U postgres)
./src/tcp-port-wait /var/run/postgresql 5432
echo $(pg_isready -h /var/run/postgresql -p 5432 -U postgres)
PGPASSWORD=password psql -h /var/run/postgresql -p 5432 -U postgres -d public_cnv -c "\dt+;"
PGPASSWORD=password pg_restore -h /var/run/postgresql -p 5432 -U postgres -d public_cnv -j 2 /data/cnv_db.sql
PGPASSWORD=password psql -h /var/run/postgresql -p 5432 -U postgres -d public_cnv -c "\dt+;"

