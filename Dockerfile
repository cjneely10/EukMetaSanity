# Get base conda image
FROM conda/miniconda3:latest

# Run initial updates/upgrades
RUN apt-get update && \
    # Create user and give permissions
    groupadd -g 999 appuser && useradd -r -u 999 -g appuser appuser && \
    # Install apt dependencies
    apt-get -y install git gcc g++ parallel wget autoconf make python2.7 libbz2-dev zip python3-pip && \
    # Create user directories
    mkdir /home/appuser && cd /home/appuser && mkdir opt bin scripts data tmp && cd - && \
    # Create .scripts file and source it
    touch /home/appuser/scripts/.scripts && \
    echo "source /home/appuser/scripts/.scripts" >> /home/appuser/.bashrc && \
    # Add bin folder to PATH
    echo "export PATH=/home/appuser/bin:\$PATH" >> /home/appuser/scripts/.scripts && \
    # EukMetaSanity
    mkdir /home/appuser/opt/EukMetaSanity

# Copy repo contents
COPY * /home/appuser/opt/EukMetaSanity/

# Installation of dependencies and adding to PATH
# EukMetaSanity install
RUN cd /home/appuser/opt/EukMetaSanity && make all && cd - && \
    echo "export PATH=$(pwd)/EukMetaSanity/bin/:\$PATH" >> /home/appuser/scripts/.scripts && \
    echo "export PYTHONPATH=$(pwd)/EukMetaSanity/:\$PYTHONPATH" >> /home/appuser/scripts/.scripts && \
    ln -s $(pwd)/EukMetaSanity/EukMetaSanity.py /home/appuser/bin/EukMetaSanity && \
    # Move to opt directory for remaining program installations
    cd /home/appuser/opt && \
    # # AUGUSTUS
    # apt dependencies
    apt-get -y install libboost-iostreams-dev zlib1g-dev libbamtools-dev libboost-all-dev libboost-all-dev && \
    apt-get -y install libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libmysql++-dev && \
    # bam2wig installation
    mkdir bam2wig && cd bam2wig && \
    git clone https://github.com/samtools/htslib.git && \
    cd htslib && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    git clone https://github.com/samtools/samtools.git && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install && \
    cd ../../ && \
    echo "export TOOLDIR=$(pwd)/bam2wig:\$PATH" >> /home/appuser/scripts/.scripts && \
    # AUGUSTUS
    git clone https://github.com/Gaius-Augustus/Augustus.git && \
    cd Augustus && mkdir -r bin src && cd src && make && cd ../auxprogs && make && cd ../../ && \
    echo "export PATH=$(pwd)/Augustus/bin:$(pwd)/Augustus/scripts:\$PATH" >> /home/appuser/scripts/.scripts && \
    echo "export AUGUSTUS_CONFIG_PATH=$(pwd)/Augustus/config/" && \
    ln -s $(pwd)/Augustus/bin/* /home/appuser/bin/ && \
    #apt-get -y install augustus augustus-data augustus-doc && \
    # # GFFread
    git clone https://github.com/gpertea/gffread && cd gffread && make release && \
    ln -s $(pwd)/gffread /home/appuser/bin/ && cd - && \
    # # GFFcompare
    git clone https://github.com/gpertea/gffcompare && cd gffcompare && make release && \
    ln -s $(pwd)/gffcompare /home/appuser/bin/ && cd - && \
    # # MMseqs2
    wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && \
    tar "xzf" mmseqs-linux-avx2.tar.gz && rm mmseqs-linux-avx2.tar.gz && \
    echo "export PATH=$(pwd)/mmseqs/bin/:\$PATH" >> /home/appuser/scripts/.scripts && \
    # # MetaEuk
    wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz && \
    tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz && \
    echo "export PATH=$(pwd)/metaeuk/bin/:\$PATH" >> /home/appuser/scripts/.scripts && \
    # # kofamscan
    wget ftp://ftp.genome.jp/pub/tools/kofam_scan/kofam_scan-1.3.0.tar.gz && \
    tar "xzvf" kofam_scan-1.3.0.tar.gz && rm kofam_scan-1.3.0.tar.gz && \
    ln -s $(pwd)/kofam_scan-1.3.0/exec_annotation /home/appuser/bin/ && \
    # # EggNOG mapper
    git clone https://github.com/jhcepas/eggnog-mapper.git && \
    ln -s $(pwd)/eggnog-mapper/emapper.py /home/appuser/bin/ && \
    # # Infernal
    wget eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar xf infernal-1.1.2.tar.gz && rm infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && ./configure --prefix $(pwd)/bin && make && make install && cd - && \
    ln -s $(pwd)/infernal-1.1.2/bin/bin/* /home/appuser/bin/ && \
    # # sambamba
    wget https://github.com/biod/sambamba/releases/download/v0.7.1/sambamba-0.7.1-linux-static.gz && \
    gunzip sambamba-0.7.1-linux-static.gz && chmod +x sambamba-0.7.1-linux-static && \
    mv sambamba-0.7.1-linux-static /home/appuser/bin/sambamba && \
    # # BEDtools
    wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz && \
    tar "-xzf" bedtools-2.29.2.tar.gz && rm bedtools-2.29.2.tar.gz && \
    cd bedtools2 && make && cd - && \
    ln -s $(pwd)/bedtools2/bin/* /home/appuser/bin/ && \
    # # HISAT2
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip && rm hisat2-2.1.0-Linux_x86_64.zip && \
    ln -s $(pwd)/hisat2-2.1.0/hisat2* /home/appuser/bin/ && \
    ln -s $(pwd)/hisat2-2.1.0/extract* /home/appuser/bin/ && \
    # # GMAP
    wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-06-30.tar.gz && \
    tar "xzf" gmap-gsnap-2020-06-30.tar.gz && rm gmap-gsnap-2020-06-30.tar.gz && \
    cd gmap-2020-06-30 && ./configure && make && cd - && \
    ln -s $(pwd)/gmap-2020-06-30/src/{gmap,gmapindex} /home/appuser/bin/ && \
    # # Add locations for RepeatModeler/Masker and GeneMark
    ln -s $(pwd)/repeatmodeler/* /home/appuser/bin/ && \
    ln -s $(pwd)/repeatmasker/* /home/appuser/bin/ && \
    ln -s $(pwd)/gmes/* /home/appuser/bin/ && \
    # # Add final permissions
    chown -R appuser:appuser /home/appuser && chmod 777 /home/appuser/data && chmod 777 /home/appuser/tmp

VOLUME ["/home/appuser/data", "/home/appuser/opt/repeatmodeler", "/home/appuser/opt/repeatmasker"]
VOLUME ["/home/appuser/opt/gmes"]
# Change user
USER appuser
# Calls MetaSanity using the passed values from the config file
ENTRYPOINT ["EukMetaSanity"]
