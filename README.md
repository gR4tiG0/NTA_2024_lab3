# NTA_2024_lab2

- download docker from remote: `docker pull gratigo/nta_lab3:latest`
- run docker from remote: `docker run --rm -it gratigo/nta_lab3:latest`
- build from git repo: `docker build -t 'gratigo/nta_lab3' .`
- run cli: `docker run --rm -it 'gratigo/nta_lab3'`
- stop container: `docker stop $(docker ps | grep "gratigo/nta_lab3" | cut -d " " -f1)`
- remove image: `docker image rm 'gratigo/nta_lab3`
