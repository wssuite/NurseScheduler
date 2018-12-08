# Guide for running the code within a docker container

1. Install docker for you platform: https://docs.docker.com/engine/installation/#supported-platforms

2. Then, just run:
````bash
docker run --rm legraina/ns
````
It will take a long time the first time to download the container. Then, you will see all the command available. For example, to run the static solver (by default, use -d to run the dynamic solver) with a given instance:
```bash
docker run --rm legraina/ns -i n005w4_0_1-2-3-3
```
Look at the docker-entrypoint.sh file for more details.

3. The results are written in outfiles/n005w4_0_1-2-3-3/unixtimestamps (unixtimestamps is computed at every run). You should also link the container to different volumes to be able to:
  1. Access the results:
  ```bash
  docker run --rm -v "$(pwd)/outfiles:~/ns/outfiles" legraina/ns -i n005w4_0_1-2-3-3
  ```
  2. Add some instances:
  ```bash
  docker run --rm -v "$(pwd)/datasets:~/ns/datasets" legraina/ns -i n005w4_0_1-2-3-3
  ```
  3. Add a configuration file:
  ```bash
  docker run --rm -v "$(pwd)/paramfiles:~/ns/paramfiles" legraina/ns -i n005w4_0_1-2-3-3
  ```


4. You may also use the docker-compose.yml file and simply run:
```bash
docker-compose up -d
```
This command will build locally the container if the field "build" is not commented, otherwise it can use the field "image". Also, all the volumes are defined in this file, you may change it according to your paths if needed.
