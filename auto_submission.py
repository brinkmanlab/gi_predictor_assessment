# -*- coding: utf-8 -*-
"""
Created on Tue May  2 11:31:13 2017

@author: Keith Tilley
"""

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
import selenium.common.exceptions as e
from selenium.webdriver.support import expected_conditions as EC
import os, time, io, subprocess

# automated submission, download and renaming/moving of result
def predictbias_submit(gbk_dir, genome):
    print("Running genome "+genome)
    driver = webdriver.Chrome() # change for different browsers
    driver.get("http://www.bioinformatics.org/sachbinfo/predictbias.html")
    assert "PredictBias" in driver.title
    print("Uploading genbank file...")
    upload = driver.find_element_by_name("brwGenome")
    upload.send_keys(os.path.abspath(gbk_dir+genome+".gbk"))
    submit = driver.find_element_by_name("btnSubmit")
    submit.send_keys(Keys.RETURN)
    print("Downloading results...")
    driver.get(driver.current_url)
    download = driver.find_element_by_link_text("download")
    link = download.get_attribute("href")
    subprocess.Popen(["wget",link])
    time.sleep(6) # allow time to download
    subprocess.Popen(["mv", link.split("/")[-1], genome+".tar.gz"])
    driver.close()
    
def paidb_submit(ffnfile, outhtml):
    driver = webdriver.Chrome() # change for different browsers
    driver.get("http://www.paidb.re.kr/pai_finder.php?m=f")
    driver.find_element_by_name("SEQFILE").send_keys(ffnfile)
    driver.find_element_by_xpath("//input[@value='Analyze']").send_keys(Keys.RETURN)
    driver.get(driver.current_url)
    try:
        WebDriverWait(driver, 120).until(
            EC.text_to_be_present_in_element((By.XPATH, "//h3[1]"), "Done !")
        )
    except e.TimeoutException:
        print(ffnfile.split("/")[-1]+" Took too long.")
    else:
        with io.open(outhtml, "w") as o:
            o.write(driver.page_source)
    driver.close()
