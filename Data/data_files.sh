#!/bin/bash
download_folder='Data/'

wget -i Data/data_urls.txt '--directory-prefix='$download_folder
