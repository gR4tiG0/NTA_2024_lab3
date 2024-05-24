FROM python:3.8-slim-buster

WORKDIR /app

RUN apt-get update && apt-get install -y gcc make libgmp-dev

ADD ./dlp.c /app/dlp.c
ADD ./dlpModule.py /app/dlpModule.py
ADD ./solve.py /app/solve.py
ADD Makefile /app/Makefile

RUN ["make"]
CMD ["python", "./solve.py"]