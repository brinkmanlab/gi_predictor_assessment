# -*- coding: utf-8 -*-
"""
Created on Tue May  2 11:31:13 2017

@author: Keith Tilley
"""

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import os, time

# function to automate submission, download and renaming/moving of results for
# PredictBias. Will only work on my computer without changes
def predictbias_submit(genome):
    print("Running genome "+genome)
    driver = webdriver.Chrome() # change for different browsers
    driver.get("http://www.bioinformatics.org/sachbinfo/predictbias.html")
    assert "PredictBias" in driver.title
    print("Uploading genbank file...")
    upload = driver.find_element_by_name("brwGenome")
    upload.send_keys("/home/beef/Documents/BrinkmanLab/GenomicIslands/GenBankFiles/"+genome+".gbk") # my genbank files location
    submit = driver.find_element_by_name("btnSubmit")
    submit.send_keys(Keys.RETURN)
    print("Downloading results...")
    driver.get(driver.current_url)
    download = driver.find_element_by_link_text("download")
    download.send_keys(Keys.RETURN)
    # The following code is to rename the download to the proper format and move
    # to my predictbias folder. Will only work on my system, with an empty dowload folder.
    time.sleep(6) # allow time for download
    result_file = os.listdir("/home/beef/Downloads/")
    if len(result_file) == 1:
        os.rename("/home/beef/Downloads/"+result_file[0], "/home/beef/Documents/BrinkmanLab/GenomicIslands/PredictBias/"+genome+".tar.gz")
    elif len(result_file) > 1:
        print("Multiple files in downloads")
    else:
        # When PredictBias fails to find any genomic islands, the download link
        # takes browser to bioinformatics main page - nothing downloaded.
        print("GENOME "+genome+" FAILED")
    driver.close()