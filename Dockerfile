FROM python:3

ADD results_database.db results_database.db
ADD metadata.json metadata.json

RUN pip install https://github.com/simonw/datasette/archive/refs/tags/1.0a2.zip && \
    pip install datasette-vega && \
    find /usr/local/lib -name '__pycache__' | xargs rm -r && \
    rm -rf /root/.cache/pip

EXPOSE 8080
ENTRYPOINT ["datasette","results_database.db","-p","8080","-h","0.0.0.0","-m","metadata.json"]