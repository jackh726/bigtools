FROM rust:buster as build

RUN mkdir /bigtools
COPY ./src /bigtools/src
COPY ./Cargo.toml /bigtools
COPY ./Cargo.lock /bigtools
RUN cd /bigtools && cargo install --path .

FROM python:3.8.6

COPY --from=build /usr/local/cargo/bin/ /usr/local/bin/

RUN git clone https://github.com/ucscGenomeBrowser/kent.git && \
    cd kent/src/lib && make && cd ../jkOwnLib && make && cd ../htslib && make && \
    mkdir -p /root/bin/x86_64 && cd ../utils/bedGraphToBigWig \
    cd ../bedGraphToBigWig && make && \
    cd ../bedToBigBed && make && \
    cd ../bigWigMerge && make && \
    cd ../bigWigAverageOverBed && make && \
    cd ../bigWigToBedGraph && make && \
    cd / && rm -rf kent && mv /root/bin/x86_64/* /bin

RUN apt-get update && apt-get install -y time && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir psutil

COPY ./bench /app

ENTRYPOINT [ "/bin/bash" ]
