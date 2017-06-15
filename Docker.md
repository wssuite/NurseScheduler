1- Install docker for you platform: https://docs.docker.com/engine/installation/#supported-platforms

2- Then, just run:
````bash
docker-compose up -d
````
It will take a long time the first time to build the container.

4- Run the following command to see the containers that running:
````bash
docker-compose ps
````

5- Run the following command to look at the log:
````bash
docker-compose logs
````

6- You can change the command in the docker-compose file to run the right instance. The software reads the data from the datafiles folder, the parameters from the paramfiles folder, and then write the output in the outfiles folder.
