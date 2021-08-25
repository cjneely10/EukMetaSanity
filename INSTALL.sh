#!/bin/bash
# Installation location
CWD=`pwd`
# Get conda location information
CONDA="`which conda`"
CONDA_DIRNAME="`dirname $CONDA`"
MINICONDA="`dirname "$CONDA_DIRNAME"`"
SOURCE="$MINICONDA"/etc/profile.d/conda.sh

# Parse command-line arguments
POSITIONAL=()
DATABASE_PATH="$CWD"
THREADS="1"
SKIP_DATA_DOWNLOAD=false
SOURCE_SCRIPT=~/.bashrc
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -t|--threads)
      THREADS="$2"
      shift
      shift
      ;;
    -s|--skip-data-download)
      SKIP_DATA_DOWNLOAD=true
      shift
      ;;
    -d|--database-path)
      DATABASE_PATH="$2"
      shift
      shift
      ;;
    -h|--help)
      POSITIONAL+=("$1")
      shift
      ;;
    -b|--bash-source-script)
      SOURCE_SCRIPT="$2"
      shift
      shift
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done
# Confirm correctly parsed
if [ ${#POSITIONAL[@]} -gt 0 ]; then
  echo ""
  echo "Usage: ./INSTALL.sh [-h] [-t <threads>] [-s] [-d /path/to/database/downloads] [-b /source/script]"
  echo ""
  echo ""
  echo "-h|--help                           Display this help message"
  echo "-t|--threads <threads>              Number of threads to use in building indices"
  echo "-s|--skip-data-download             Skip repeat modeler and EukMS database download (otherwise, uses wget)"
  echo "-d|--database-path <path>           Path for database download, default is $CWD"
  echo "-b|--bash-source-script <path>      Script to add PATH updates, default is ~/.bashrc"
  echo ""
  exit 1
fi

# Usage: install_env <env-name> <deactivate?>
function install_env() {
  EXISTS=$(conda info --envs | grep -c "EukMS_$1")
  if [ $EXISTS -lt 1 ]; then
    mamba env create -f bin/$1-pipeline/$1/environment.yml
  fi
  source $SOURCE
  conda activate EukMS_$1
  python -m pip install .
  if [ $2 = true ]; then
    conda deactivate
  fi
}

# Create repeats.txt installation instructions with download path
function modify_rm_location() {
  # Install RepeatMasker updated libraries and configure
  cd "$MINICONDA"/envs/EukMS_run/share/RepeatMasker/Libraries/
  if [ $1 = false ]; then
    wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz
    gunzip Dfam.h5.gz
  fi
  cd ..
  sed "s,INSTALLATION_LOCATION,$2," "$CWD"/install/repeats.default.txt > "$CWD"/install/repeats.txt
  perl ./configure < "$CWD"/install/repeats.txt
  cp util/rmOutToGFF3.pl ./
}

# Update conda versions of augustus for proper id parsing
function update_augustus() {
  cd "$MINICONDA"/envs/EukMS_run/bin
  sed -i 's/transcript_id \"(\.\*)\"/transcript_id \"(\\S\+)"/' filterGenesIn_mRNAname.pl
  cd "$MINICONDA"/envs/EukMS_refine/bin
  sed -i 's/transcript_id \"(\.\*)\"/transcript_id \"(\\S\+)"/' filterGenesIn_mRNAname.pl
}

# Download MetaEuk from MMSeqs2 server
function install_metaeuk() {
  # Create bin directory and install non-conda dependencies
  mkdir -p bin
  cd bin || return 1
  wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
  tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
  mv metaeuk/bin/metaeuk ./_metaeuk
  rm -r metaeuk/ && mv _metaeuk metaeuk && cd - || return 1
  return 0
}

# Create run environment and install repeat modeler updates
function install_eukms_run() {
  install_env run false
  modify_rm_location $1 "$MINICONDA"/envs/EukMS_run/bin/
  update_augustus
  # Return to installation directory
  cd $CWD
  conda deactivate
}

# Add EukMS definitions to source script (default is ~/.bashrc)
function update_source_script() {
  echo export PATH="$(pwd)"/bin/:'$PATH' >> "$1"
  echo export EukMS_run="$2" >> "$1"
  echo export EukMS_report="$(pwd)"/bin/report-pipeline >> "$1"
  echo export EukMS_refine="$(pwd)"/bin/refine-pipeline >> "$1"
}

# # # Begin installation procedure

# Install mamba if not already present
conda install mamba -n base -c conda-forge -y

# Create run environment
install_eukms_run $SKIP_DATA_DOWNLOAD
# Create report environment
install_env report true
# Create refine environment
install_env refine true

# Install metaeuk and confirm proper completion
install_metaeuk
if [ $? = 1 ]; then
  exit 1
fi

# Update .bashrc with proper locations
EukMS_run="$(pwd)"/bin/run-pipeline
update_source_script "$SOURCE_SCRIPT" "$EukMS_run"

# Download updated databases
if [ $SKIP_DATA_DOWNLOAD = false ]; then
  conda activate EukMS_run
  download-data -t $THREADS -d "$DATABASE_PATH" --eukms-run-bin "$EukMS_run"
  conda deactivate
fi

echo "Your installation is complete!"
