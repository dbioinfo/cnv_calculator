pg_dump -U postgres -h localhost -p 5432 -f data/cnv_db.sql

tar -cvfz data/shipping_container.tar.gz -T data/cnv_filelist.txt

scp -i ec2-cnv-key.pem data/shipping_container.tar.gz ubuntu@ec2-3-137-140-127.us-east-2.compute.amazonaws.com:/
