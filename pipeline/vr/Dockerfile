FROM ubuntu:18.04
FROM python:3.7

RUN chmod 1777 /tmp

WORKDIR /opt

RUN pip3 install -r requirements.txt

RUN mkdir /data && chmod -R o+rwx /data

ENTRYPOINT ["appendVRIds.py"]

