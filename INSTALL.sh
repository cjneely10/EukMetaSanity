#!/bin/bash
set -eo pipefail

# EukMS pipelines to generate
PIPELINES=(run report refine)
# CLI arguments
THREADS="1"
SOURCE_SCRIPT=~/.bashrc
UNINSTALL=false
UPGRADE=false

# Installation location
CWD="$(pwd)"
# Conda location information
CONDA="$(which conda)"
CONDA_DIRNAME="$(dirname "$CONDA")"
MINICONDA="$(dirname "$CONDA_DIRNAME")"
SOURCE="$MINICONDA/etc/profile.d/conda.sh"
# Installation directories
BIN=bin
DATA=data

function usage() {
  echo ""
  echo "EukMetaSanity v1.1.0 Installation"
  echo "Usage: ./INSTALL.sh [-h] | [[-t <threads>] [-b /source/script]] | [--upgrade [-t <threads>]] | [--uninstall]"
  echo ""
  echo ""
  echo "-h|--help                           Display this help message"
  echo "-t|--threads <threads>              Number of threads to use in building indices"
  echo "-b|--bash-source-script <path>      Script to add PATH updates, default is ~/.bashrc"
  echo "--uninstall                         Uninstall EukMetaSanity"
  echo "--upgrade                           Upgrade EukMetaSanity"
  echo ""
  echo ""
  echo "Examples"
  echo "--------"
  echo "Basic Installation"
  echo "./INSTALL.sh -t 8"
  echo "Place EukMetaSanity PATH updated in non-default location (i.e., not ~/.bashrc)"
  echo "./INSTALL.sh -t 8 -b ~/.scripts"
  echo "Upgrade installation"
  echo "./INSTALL.sh -t 6 --upgrade"
  echo "Uninstall EukMetaSanity"
  echo "./INSTALL.sh --uninstall"
  echo ""
}

# Check if conda environment exists
function is_installed() {
  [ "$(conda info --envs | grep -c "$1")" -eq 1 ]
}

# Check if conda environment doesn't exist
function not_installed() {
  [ "$(conda info --envs | grep -c "$1")" -eq 0 ]
}

# Uninstall all EukMetaSanity conda environments. Delete bin and data directories
function uninstall() {
  # Remove conda environments
  for f in "${PIPELINES[@]}"; do
    if is_installed "EukMS_$f"; then
      conda remove --name "EukMS_$f" --all -y
    fi
  done
  # Remove installer, if not already done so
  if is_installed yapim_installer; then
    conda remove --name yapim_installer --all -y
  fi
  rm -rf "$BIN"
  rm -rf "$DATA"
}

# Generate EukMetaSanity bin pipeline directories
function create_binary_directory() {
  local DO_INSTALL=false
  for f in "${PIPELINES[@]}"; do
    if "$UPGRADE" || [ ! -d "$BIN/$f-pipeline" ]; then
      DO_INSTALL=true
      break
    fi
  done
  if "$DO_INSTALL"; then
    # Create and activate installer environment
    if not_installed yapim_installer; then
      mamba env create -f EukMetaSanity/src/installer-environment.yml
    fi
    source "$SOURCE"
    conda activate yapim_installer
    # Generate binary directories using the YAPIM `create` function
    for f in "${PIPELINES[@]}"; do
      if "$UPGRADE" || [ ! -d "$BIN/$f-pipeline" ]; then
        echo n | yapim create -t "EukMetaSanity/src/$f" -d EukMetaSanity/src/dependencies -o "$BIN/$f-pipeline"
      fi
    done
    # Deactivate and delete installer environment
    conda deactivate
    conda remove --name yapim_installer --all -y
  fi
}

# Install mamba if not already present
function install_mamba() {
  if [ "$(conda list | grep -c mamba)" -lt 1 ]; then
    conda install mamba -n base -c conda-forge -y
  fi
}

# Install EukMS package to environment
function pip_install_env() {
  source "$SOURCE"
  conda activate "EukMS_$1"
  python -m pip install .
  conda deactivate
}

# Usage: install_env <env-name>
function install_env() {
  local ENV_FILE="$BIN/$1-pipeline/$1/environment.yml"
  if not_installed "EukMS_$1"; then
    mamba env create -f "$ENV_FILE"
    pip_install_env "$1"
  elif "$UPGRADE"; then
    mamba env update -f "$ENV_FILE"
    pip_install_env "$1"
  fi
}

# Create repeats.txt installation instructions with download path
# Usage: modify_rm_location <install-path>
function modify_rm_location() {
  # Install RepeatMasker updated libraries and configure
  cd "$MINICONDA/envs/EukMS_run/share/RepeatMasker/Libraries/"
  # Download data if $UPGRADE is set or if file 'installation-complete' is not present
  if "$UPGRADE" || [ ! -f ../installation-complete ]; then
    source "$SOURCE"
    conda activate EukMS_run
    rm -f Dfam.h5 Dfam.h5.gz
    wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz
    gunzip Dfam.h5.gz
    # Move up to RepeatMasker top-level
    cd ..
    sed "s,INSTALLATION_LOCATION,$1," "$CWD/install/repeats.default.txt" > "$CWD/install/repeats.txt"
    # Configure installation. If download is skipped, uses smaller, less-complete repeat database
    perl ./configure < "$CWD/install/repeats.txt"
    rm "$CWD/install/repeats.txt"
    cp util/rmOutToGFF3.pl ./
    conda deactivate
    touch installation-complete
  fi
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
  cd "$BIN"
  source "$SOURCE"
  conda activate "EukMS_run"
  wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
  conda deactivate
  tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
  mv "metaeuk/$BIN/metaeuk" ./_metaeuk
  rm -r metaeuk/ && mv _metaeuk metaeuk
  cd -
}

# Create run environment and install repeat modeler updates
function install_eukms_run() {
  install_env run
  modify_rm_location "$MINICONDA/envs/EukMS_run/bin/"
}

# Add EukMS definitions to source script (default is ~/.bashrc)
function update_source_script() {
  LINE1="# Environment variables are part of EukMetaSanity installation."
  if [ "$(grep -c "$LINE1" "$SOURCE_SCRIPT")" -lt 1 ]; then
    echo "$LINE1" >> "$1"
    echo "# Remove on program deletion" >> "$1"
    echo export PATH="$(pwd)/$BIN/":'$PATH' >> "$1"
    echo export EukMS_run="$2" >> "$1"
    echo export EukMS_report="$(pwd)/$BIN/report-pipeline" >> "$1"
    echo export EukMS_refine="$(pwd)/$BIN/refine-pipeline" >> "$1"
    echo "# # # # # #" >> "$1"
  fi
}

# # # # # Begin # # # # #

# Validate top-level working directory
if [ "$(dirname "$0")" != . ]; then
  >&2 echo "Run installation script from EukMetaSanity top-level directory"
  usage
  exit 1
fi

# Validate conda source path
if [ ! -f "$SOURCE" ]; then
  >&2 echo "There was an error locating your conda installation"
  exit 1
fi

# Parse command-line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -t|--threads)
      THREADS="$2"
      shift
      shift
      ;;
    -h|--help)
      usage
      exit 0
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
    --upgrade)
      UPGRADE=true
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
  usage
  exit 1
fi

# Confirm valid upgrade/uninstall setting
if "$UPGRADE" && "$UNINSTALL"; then
  >&2 echo "Provide either --upgrade OR --uninstall flag"
  exit 1
fi

# # # Uninstall EukMetaSanity # # #
if "$UNINSTALL"; then
  uninstall
  exit 0
fi

# # # Install EukMetaSanity # # #
# Installer
install_mamba
# Generate installation directory
create_binary_directory
# Create run environment
install_eukms_run
# Create report environment
install_env report
# Create refine environment
install_env refine

# Install metaeuk and confirm proper completion
if "$UPGRADE" || [ ! -e "$BIN/metaeuk" ]; then
  install_metaeuk
fi

# Configure augustus version
update_augustus

# Update .bashrc with proper locations
EukMS_run="$CWD/$BIN/run-pipeline"
update_source_script "$SOURCE_SCRIPT" "$EukMS_run"

# Download updated databases
if "$UPGRADE" || [ ! -d "$DATA" ]; then
  source "$SOURCE"
  conda activate EukMS_run
  download-data -t "$THREADS" --eukms-run-bin "$EukMS_run"
  conda deactivate
fi

echo "Your installation is complete!"
