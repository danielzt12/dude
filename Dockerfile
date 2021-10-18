FROM debian:10

COPY . /app

RUN apt-get update && apt-get install -y \
 python3-gi python3-gi-cairo gir1.2-gtk-3.0 python3-pip

RUN rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN python3 -m pip install -r requirements.txt

CMD ["python3", "-u", "dude.py"]