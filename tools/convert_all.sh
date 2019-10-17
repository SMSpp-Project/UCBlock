#!/usr/bin/env bash

EXE_PATH="../../cmake-build-debug/UCBlock/tools/nc4generator"
INPUT_DIR="../../UCBlock/data"
OUTPUT_DIR="../../UCBlock/netCDF_files"

cd "${INPUT_DIR}" || exit
for INPUT_FILE in $(find . -name "*.dat" -or -name "*.mod"); do
    OUTPUT_FILE=${INPUT_FILE/%???/nc4}
    outdir=$(dirname "${OUTPUT_FILE}")
    echo "Converting ${INPUT_FILE} into ${OUTPUT_FILE}";
    eval "${EXE_PATH} ${INPUT_FILE}";
    mkdir -p "${OUTPUT_DIR}/${outdir}/"
    mv "${OUTPUT_FILE}" "${OUTPUT_DIR}/${outdir}/";
done
