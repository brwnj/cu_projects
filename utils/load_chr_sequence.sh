#!/usr/bin/env bash

## load_chr.sh: load a single chromosome (or all chromosomes) into genomedata

## environment must provide
##
## AGPDIR: directory with assembly golden path (agp) files for each chromosome
## GENOMEDATADIR: output directory
## BED_FILENAMES: Bash array of bed files (see `info "(bash) Arrays"')
## BED_TRACKNAMES: Bash array of track names, matched to BED_FILENAMES

## $Revision: 5975 $
## Copyright 2010-2012 Michael M. Hoffman <mmh1@uw.edu>

set -o nounset -o pipefail -o errexit -x

if [ $# -gt 2 ]; then
    echo usage: "$0" [COMMON_SCRIPT] [CHR]
    exit 2
fi

common_script="${1:-/dev/null}"
chr="${2:-*}" # default is all

# this has to be loaded here because it sets array variables which cannot be exported
source "$common_script"

num_tracks="${#BED_FILENAMES[@]}"

genomedata-load --verbose --directory-mode \
    --sequence="$SEQDIR/${chr}.fa.gz" \
    $(for ((track_index=0; track_index < num_tracks; track_index++)); do
        printf "%s %s=%s " "-t" "${BED_TRACKNAMES[$track_index]}" "${BED_FILENAMES[$track_index]}"
        done) \
    "$GENOMEDATADIR"
