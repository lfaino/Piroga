FROM ubuntu:20.04
LABEL authors="lfaino"
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget gnupg2 software-properties-common curl
RUN apt-get update && apt-get -y install jellyfish spades samtools nanopolish python3 git cmake g++ zlib1g-dev bwa

RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2 && \
    tar -jxvf minimap2-2.28.tar.bz2 && cd minimap2-2.28 && make && cp minimap2 /usr/bin



RUN git clone  --recursive https://github.com/lbcb-sci/racon.git  && cd racon && mkdir build && cd build &&  \
    cmake -DCMAKE_BUILD_TYPE=Release .. && make -j70 && make install

RUN git clone https://github.com/marbl/Krona.git && cd Krona/KronaTools && ./install.pl

RUN apt install -y python3-pip

RUN pip3 install biopython

RUN pip3 install tqdm

RUN apt-get install -y python3-matplotlib

RUN apt-get install -y libxmu-dev libxmu-headers freeglut3-dev libxext-dev libxi-dev libhts-dev zlib1g-dev bamtools

RUN git clone https://github.com/rrwick/Porechop.git && cd Porechop && python3 setup.py install

#RUN git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make

RUN apt-get install -y tabix

ENV BLASTDB = "/data/blastdb"

CMD ["bash"]
