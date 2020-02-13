FROM python:3.7

RUN chmod 1777 /tmp

WORKDIR /opt

COPY ./requirements.txt /opt/requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# we copy in app code here so that changing the app code doesn't invalidate the previously cached layers
# (in practice, this means we can skip redownloading all the python requirements on each code change)
ADD . .

RUN mkdir /data && chmod -R o+rwx /data

ENTRYPOINT ["python", "appendVRIds.py"]
