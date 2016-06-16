#!/bin/bash
download_folder='data/'

wget -i data_urls.txt '--directory-prefix='$download_folder 
