FROM python:3.7

RUN chmod 1777 /tmp

WORKDIR /opt

ADD . .

RUN pip install --upgrade pip

# install numpy first to avoid issues with bio python and bx-python (see also https://github.com/LUMC/vep2lovd/issues/1)
RUN pip install $(grep numpy requirements.txt)

RUN pip install -r requirements.txt

RUN mkdir /data && chmod -R o+rwx /data

ENTRYPOINT ["python", "appendVRIds.py"]
