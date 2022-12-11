#!/bin/bash
set -exo pipefail

# Installation location
CWD="$(pwd)"
# Get conda location information
CONDA="$(which conda)"
CONDA_DIRNAME="$(dirname "$CONDA")"
MINICONDA="$(dirname "$CONDA_DIRNAME")"
SOURCE="$MINICONDA/etc/profile.d/conda.sh"

# Parse command-line arguments
POSITIONAL=()
THREADS="1"
SKIP_RM_DOWNLOAD=false
SKIP_DATA_DOWNLOAD=false
SOURCE_SCRIPT=~/.bashrc
UNINSTALL=false
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -t|--threads)
      THREADS="$2"
      shift
      shift
      ;;
    -r|--skip-rm-download)
      SKIP_RM_DOWNLOAD=true
      shift
      ;;
    -s|--skip-database-download)
      SKIP_DATA_DOWNLOAD=true
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
    --uninstall)
      UNINSTALL=true
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
  echo "Usage: ./INSTALL.sh [-h] [-t <threads>] [-r] [-s] [-d /path/to/database/downloads] [-b /source/script]"
  echo ""
  echo ""
  echo "-h|--help                           Display this help message"
  echo "-t|--threads <threads>              Number of threads to use in building indices"
  echo "-r|--skip-rm-download               Skip repeat modeler database download (otherwise, uses wget)"
  echo "-s|--skip-database-download         Skip EukMS database download (otherwise, uses wget)"
  echo "-b|--bash-source-script <path>      Script to add PATH updates, default is ~/.bashrc"
  echo "--uninstall                         Uninstall EukMetaSanity"
  echo ""
  exit 1
fi

if [ ! -e "$SOURCE" ]; then
  echo "There was an error locating your conda installation"
  exit 1
fi

function uninstall() {
  for f in run report refine; do
    conda remove --name EukMS_$f --all -y
  done
  if [ "$(conda info --envs | grep -c yapim_installer)" -eq 1 ]; then
    conda remove --name yapim_installer --all -y
  fi
  conda remove mamba -y
}

function create_binary_directory() {
  if [ ! -e "$CWD/bin" ]; then
    mamba env create -f EukMetaSanity/src/installer-environment.yml
    source "$SOURCE"
    conda activate yapim_installer
    make
    conda deactivate
    conda remove --name yapim_installer --all -y
  fi
}

function install_mamba() {
  # Install mamba if not already present
  if [ "$(conda list | grep -c mamba)" -lt 1 ]; then
    conda install mamba -n base -c conda-forge -y
  fi
}

# Usage: install_env <env-name> <deactivate?>
function install_env() {
  if [ "$(conda info --envs | grep -c "EukMS_$1")" -lt 1 ]; then
    mamba env create -f "bin/$1-pipeline/$1/environment.yml"
    source "$SOURCE"
    conda activate "EukMS_$1"
    python -m pip install .
    if [ "$2" = true ]; then
      conda deactivate
    fi
  fi
}

# Create repeats.txt installation instructions with download path
function modify_rm_location() {
  # Install RepeatMasker updated libraries and configure
  cd "$MINICONDA/envs/EukMS_run/share/RepeatMasker/Libraries/"
  if [ "$1" = false ]; then
    wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz
    gunzip Dfam.h5.gz
  fi
  cd ..
  sed "s,INSTALLATION_LOCATION,$2," "$CWD/install/repeats.default.txt" > "$CWD/install/repeats.txt"
  perl ./configure < "$CWD/install/repeats.txt"
  rm "$CWD/install/repeats.txt"
  cp util/rmOutToGFF3.pl ./
  cd "$CWD"
}

# Update conda versions of augustus for proper id parsing
function update_augustus() {
  cd "$MINICONDA/envs/EukMS_run/bin"
  sed -i 's/transcript_id \"(\.\*)\"/transcript_id \"(\\S\+)"/' filterGenesIn_mRNAname.pl
  cd "$MINICONDA/envs/EukMS_refine/bin"
  sed -i 's/transcript_id \"(\.\*)\"/transcript_id \"(\\S\+)"/' filterGenesIn_mRNAname.pl
  cd "$CWD"
}

# Download MetaEuk from MMSeqs2 server
function install_metaeuk() {
  # Create bin directory and install non-conda dependencies
  cd bin
  wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
  tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
  mv metaeuk/bin/metaeuk ./_metaeuk
  rm -r metaeuk/ && mv _metaeuk metaeuk && cd -
}

# Create run environment and install repeat modeler updates
function install_eukms_run() {
  if [ "$(conda info --envs | grep -c "EukMS_run")" -lt 1 ]; then
    install_env run false
    modify_rm_location "$1" "$MINICONDA/envs/EukMS_run/bin/"
    conda deactivate
  fi
}

# Add EukMS definitions to source script (default is ~/.bashrc)
function update_source_script() {
  LINE1="# Environment variables are part of EukMetaSanity installation."
  if [ "$(grep -c "$LINE1" "$SOURCE_SCRIPT")" -lt 1 ]; then
    echo "$LINE1" >> "$1"
    echo "# Remove on program deletion" >> "$1"
    echo export PATH="$(pwd)"/bin/:'$PATH' >> "$1"
    echo export EukMS_run="$2" >> "$1"
    echo export EukMS_report="$(pwd)/bin/report-pipeline" >> "$1"
    echo export EukMS_refine="$(pwd)/bin/refine-pipeline" >> "$1"
    echo "# # # # # # " >> "$1"
  fi
}

# # # Uninstall EukMetaSanity
if [ "$UNINSTALL" = true ]; then
  uninstall
  exit 0
fi

# # # Begin installation procedure

# Installer
install_mamba
# Generate installation directory
create_binary_directory
# Create run environment
install_eukms_run $SKIP_RM_DOWNLOAD
# Create report environment
install_env report true
# Create refine environment
install_env refine true

# Install metaeuk and confirm proper completion
if [ ! -e bin/metaeuk ]; then
  install_metaeuk
fi

# Configure augustus version
update_augustus

# Update .bashrc with proper locations
EukMS_run="$CWD/bin/run-pipeline"
update_source_script "$SOURCE_SCRIPT" "$EukMS_run"

# Download updated databases
if [ $SKIP_DATA_DOWNLOAD = false ]; then
  source "$SOURCE"
  conda activate EukMS_run
  if [ ! -e "$CWD/data" ]; then
    download-data -t "$THREADS" --eukms-run-bin "$EukMS_run"
  fi
  conda deactivate
fi

echo "Your installation is complete!"
