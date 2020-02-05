FROM python:3.7

RUN chmod 1777 /tmp

WORKDIR /opt

ADD . .

RUN pip install --upgrade pip

RUN pip install -r requirements.txt

RUN mkdir /data && chmod -R o+rwx /data

ENTRYPOINT ["python", "appendVRIds.py"]
