#!/bin/bash

# Download all required data to initialize the pipeline
download_folder='Data/'

wget -i Data/data_urls.txt '--directory-prefix='$download_folder
unzip Data/human.filtered.dat.zip -d $download_folder

# Check integrity of data download
md5sum -c Data/md5sums.txt

