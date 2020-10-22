FROM rust:buster as build

COPY . /bigtools
RUN cd /bigtools && cargo install --path .

FROM python:3.8.6

COPY --from=build /usr/local/cargo/bin/ /usr/local/bin/

RUN pip install --no-cache-dir psutil

ENTRYPOINT [ "/bin/bash" ]
