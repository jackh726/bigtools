# Note: this should be build from the parent directory

FROM rust:buster as build

RUN mkdir /bigtools
COPY ./bigtools/src /bigtools/src
COPY ./bigtools/Cargo.toml /bigtools
COPY ./Cargo.lock /bigtools
RUN cd /bigtools && cargo install --path .

FROM python:3.8.6

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

COPY --from=build /usr/local/cargo/bin/ /usr/local/bin/

COPY ./bench /app

ENTRYPOINT [ "/bin/bash" ]
