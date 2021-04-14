#!/bin/bash

while getopts s:p:d: flag
do
  case "${flag}" in
    s) source_path=${OPTARG};;
    p) file_pattern=${OPTARG};;
    d) dest_path=${OPTARG};;
  esac
done

if [ "$source_path" = "" ]
  then
    source_path="."
    echo "Using current directory as a source directory"
fi

if [ "$file_pattern" = "" ]
  then
    file_pattern=".*"
    echo "Pattern is not specified, using all files"
fi

if [ "$dest_path" = "" ]
  then
    dest_path="."
    echo "Using current directory as a destination"
fi
 
echo 'Source directory' ${source_path}
echo 'Destination'  ${dest_path}
echo 'File pattern' ${file_pattern}

files=$(ls -d ${source_path}/* | grep ${file_pattern})

mkdir -p ${dest_path}

ln -s ${files} -t ${dest_path}

# usage
# ./src/sh/utils.sh -s ${sourcepath}/ -d ./dat/raw/aln/fna -p .*\\.fna | less
