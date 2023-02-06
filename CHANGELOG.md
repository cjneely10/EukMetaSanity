# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2023-3-1
### Added
- Support to upgrade and uninstall `EukMetaSanity` from `INSTALL.sh` script

### Changed
- Installation script now fails on error (as expected)
- The `bin/` directory is no longer version-controlled 
  - Will auto-generate as part of installation process
- Update to YAPIM v0.3.2
- If available, augustus uses GeneMark-derived proteins for first-pass predictions
  - Otherwise, performs legacy taxonomy search and uses closest augustus species for first-pass

### Fixed
- The `report` pipeline failed to annotate results of `refine` pipeline
- Parsing existing pipeline results failed on non-uniform recorded result directories
- Transcriptome-related `Task`s defined implementations for incorrect versions of GMAP
- BRAKER2 now populates expected output
- Taxonomy parsing algorithm was incorrect
- Redundant MMSeqs databases were created
- Reduce temporary file retention

### Removed
- Ability to set database directory as part of installation
- Ability to skip database downloads as part of installation
- Ability to skip RepeatMasker database downloads as part of installation

## [1.0.0] - 2021-07-25
First release of `EukMetaSanity`
