# Get base conda image
FROM conda/miniconda3:latest

# Run initial updates/upgrades
RUN apt-get update &&
    # Create user and give permissions
    groupadd -g 999 appuser && useradd -r -u 999 -g appuser appuser &&
    # Create user directories
    mkdir -p /home/appuser/{opt,bin,scripts} &&
    # Create .scripts file and source it
    touch /home/appuser/scripts/.scripts &&
    echo "source /home/appuser/scripts/.scripts" >> /home/appuser/.bashrc &&
    # Add bin folder to PATH
    echo "PATH=/home/appuser/bin:\$PATH" >> /home/appuser/scripts/.scripts

# Copy in extra files
COPY trf409.linux64 /home/appuser/bin/

# Begin installations and add to path (either add to bin or add to PATH variable)
# apt installations
RUN apt-get -y install git gcc g++ parallel wget autoconf make &&
    # Install PERL
    conda install -c anaconda perl &&
    # Move to opt directory for program installations
    cd /home/appuser/opt &&
    # GFFread
    git clone https://github.com/gpertea/gffread && cd gffread && make release &&
    ln -s $(pwd)/gffread /home/appuser/bin/ && cd - &&
    # GFFcompare
    git clone https://github.com/gpertea/gffcompare && cd gffcompare && make release &&
    ln -s $(pwd)/gffcompare /home/appuser/bin/ && cd - &&
    # MMseqs2
    wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz &&
    tar -xvfz mmseqs-linux-avx2.tar.gz &&
    echo "PATH=$(pwd)/mmseqs/bin/:\$PATH" >> /home/appuser/scripts/.scripts &&
    # MetaEuk
    wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz &&
    tar -xvfz metaeuk-linux-sse41.tar.gz &&
    echo "PATH=$(pwd)/metaeuk/bin/:\$PATH" >> /home/appuser/scripts/.scripts &&
    # RepeatMasker
    # Dependencies
    # RMBlast
    wget http://www.repeatmasker.org/rmblast-2.10.0+-x64-linux.tar.gz &&
    tar -zxvf rmblast-2.10.0-x64-linux.tar.gz &&
    ln -s "$(pwd)/rmblast-2.10.0/bin/*" /home/appuser/bin/
