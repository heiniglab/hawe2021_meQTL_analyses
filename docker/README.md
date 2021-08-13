## README

The provided (DOCKERFILE)[DOCKERFILE] in this repository was used to run most of the anlayses (unless stated otherwise).

Below we provide a short description on how we used docker containers with charliecloud.

### Details on building and exporting docker containers

#### Building the docker containers
You have to have docker installed on your system and the docker service/daemon needs to be running.

Here is an example how you could build a container (we provide the path to the dockerfile via '-f <./path/to/file>' and a nice name via '-t <name>'):

```{bash}
docker build -f DOCKERFILE -t exported_image .
```

Don't forget the `.` at the end of the command, otherwise docker will complain.

> NOTE: Avoid putting large files in the directory from which you start the build. Docker will send **all files** in the workding directory to the build context which can take a long while. Therefore, also export your images to a different directory outside of the current one.

> NOTE: Tested on a Windows 10 machine with WLS, the above command had to be called twice in order to successfully build the container, the first call resulted in an error.

#### Exporting containers to charliecloud
If you want to use the container in a charliecloud cluster environment, you have to package it as an charliecloud-readable container. 
To achieve this, you can call the following commands:

```{}
id=$(docker create DOCKERFILE)
docker exported_image $id > exported_image.tar
gzip -c exported_image.tar > exported_image.tar.gz
```

The resulting `tar.gz` file can then be transferred to a server and run via charliecloud.

The above commands can also be combined to make it more straight forward to use:

```{bash}
docker export $(docker create DOCKERFILE) | gzip -c > exported_image.tar.gz
```