from concurrent.futures import process
from pathlib import Path
import pdb
import pickle
import re
import logging as log
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from Bio.PDB import PDBParser, PDBList
import nglview as nv
from pymol import cmd

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support import expected_conditions as EC



def pickleObject(alignment:object, name:str)-> None:
    # Save object in a file 
    with open(f'{name}.pkl', 'wb') as outp:
        pickle.dump(alignment, outp, pickle.HIGHEST_PROTOCOL)

def unpickleObject(path:str=None, Path:Path=None)-> None:
    # Load Object from file
    with open(path, 'rb') as outp:
        return pickle.load(outp)


amm=[ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def _replace(ch:str)->str:
    ''' 
    Utility function for replacing not aminoacids characters
        Param:
            ch| character in input
        Return:
            a| new character
    '''
    if ch in amm:
        return ch
    else:
        a=list(sorted(
                amm, key=lambda x:abs(ord(ch.upper())-ord(x))
            ))[0]
        return a

def name2protein(name:str)->str:
    '''
    Function for converting name 2 polypeptide sequence
    
    Input:
        name| input name
    Return:
        sequence| motif
    '''
    name=re.sub(r'\W+', '', name)
    newname=list(map(_replace, [*name]))
    newname="".join(newname)
    return newname

def get_alignment(sequence: str):
    '''
    Execute blast alignment 

    Input:
        sequence| sequence in input
    Return:
        blast_record| specific result object 
    '''
    log.info("get_alignment: Blast job started")
    #result_handle=NCBIWWW.qblast("blastp","pdb",sequence, auto_format=True)
    result_handle=NCBIWWW.qblast("blastp","nr",sequence, auto_format=True)
    blast_record=SearchIO.read(result_handle,"blast-xml")
    pickleObject(blast_record,"backup.txt")
    log.info("get_alignment: Blast result received")
    log.debug("get_alingment: Blast results\n{blast_record}\n\n"+50*"-")
    return blast_record

def blastrecors2pdbid(blast_record):
    '''
    Parse blast_record
    
    Input:
        blast_record| object
    Return:
        pdb_id| pdb_id
    '''
    hit_id=blast_record[1].blast_id #"pdb|6PXV|D"
    pdb_id=hit_id.split("|")[1]
    return pdb_id
    
def _downloadpdbfile(pdb_id:str):
    pdblist_manager=PDBList()
    log.debug("downloadpdbfiles: downloading")
    pdblist_manager.retrieve_pdb_file(pdb=[pdb_id])
    
    
def getStructure(pdb_id:str):
    _downloadpdbfile(pdb_id)
    parser=PDBParser()
    return parser.get_structure(pdb_id,f"{pdb_id}.pdb")


def showStructureNV(structure):
    nv.show_biopython(structure, gui=True)

def showStructurePyMOL(structura=None,pdb_id=None):
    '''
    Show structure using pymodule scripting
    
    Input:
    pdb_id| pdb id of protein
    '''
    #cmd.fragment('ala')
    cmd.fetch(pdb_id)
    cmd.spectrum()
    cmd.zoom()
    cmd.png('tmp/test.png', 1920, 1080)
    print("imagine saved")



# SELENIUM WEB SCRAPING
def _generate_url_searchmotif(motif:str)->str:
    url="https://www.rcsb.org/search?request="
    query:dict={
            "query": {
                "type": "terminal",
                "service": "seqmotif",
                "parameters": {
                "value": motif,
                "pattern_type": "prosite",
                "sequence_type": "protein"
                }
            },
        "return_type": "polymer_entity"
        }
    url=url+str(query)
    url=url.replace(" ","").replace("'",'"')
    return url

def _generate_url_structure(pdb_id:str):
  url=f"https://www.rcsb.org/structure/{pdb_id}"
  return url

def find_protein_pdb(driver:webdriver ,amm_seq:str, )-> str:
    '''
    Given an amminoacid sequence, search it on PDB as a motif structure. 
    Return first result.

    Input:
        driver| selenium webdriber.Chrome
        amm_seq| amminoacid sequence
    Output:
        pdb_id| first result
    '''
    #options = Options()
    #options.binary_location = "/usr/sbin/google-chrome-stable"
    ##options.add_argument("user-data-dir=/home/matteo/cache")
    #service= Service()
    #service.path='/home/matteo/chromedriver_linux64/chromedriver'
    #driver = webdriver.Chrome(options = options, service=service)
    driver.get(_generate_url_searchmotif(amm_seq))
    try:
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.CLASS_NAME, "results-item"))
        )
        elements=driver.find_elements(By.TAG_NAME,"h3")
        pdb_id=re.findall(r"[A-Z\d]{4}",elements[0].accessible_name)[0]
    except NoSuchElementException:
        print("not found")

    ##driver.find_element((By.CLASS_NAME, "results-item"))
    for el in driver.find_elements(By.TAG_NAME,"h3"):
        print(re.findall(r"[A-Z\d]{4}",el.accessible_name)[0])
    #print(driver.page_source.count("results-item"))
    try:
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.CLASS_NAME, "results-item"))
        )
        elements=driver.find_elements(By.TAG_NAME,"h3")
        pdb_id=re.findall(r"[A-Z\d]{4}",elements[0].accessible_name)[0]
    except NoSuchElementException:
        print("not found")
    print("finished ")
    driver.get(_generate_url_structure(pdb_id))
    print("quitted")
    return pdb_id


if __name__=="__main__":
    # Input protein sequence
    sequence=name2protein(input("Inserisci nome=\t"))

    # Setup Selenium
    options = Options()
    options.binary_location = "/usr/sbin/google-chrome-stable"
    options.headless=True
    #options.add_argument= "--headless"
    #options.add_argument("user-data-dir=/home/matteo/cache")
    service= Service()
    service.path='/home/matteo/chromedriver_linux64/chromedriver'
    driver = webdriver.Chrome(options = options, service=service)
    
    # Free Resources
    pdb_id=find_protein_pdb(driver, sequence)
    #showStructurePyMOL(pdb_id=pdb_id)

    try:
        driver.close()
    except:
        print("not closed")
    
    exit(5)

    # Start Search
    #piklereceived=unpickleObject("blast_record.pkl")
    #hit_id=piklereceived[1].blast_id #"pdb|6PXV|D"
    #pdb_id=hit_id.split("|")[1]
    #getStructure(pdb_id)
    #blast_record=get_alignment('''>sp|P01308|INS_HUMAN OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
    #MALWMRLLPL LALLALWGPD PAAAFVNQHL CGSHLVEALY LVCGERGFFY TPKTRREAED
    #LQVGQVELGG GPGAGSLQPL ALEGSLQKRG IVEQCCTSIC SLYQLENYCN''')
    #pickleObject(blast_record, "backup/blastrecord")
    #pdb_id=blastrecors2pdbid(blast_record)
    #showStructurePyMOL(pdb_id)

    #sequence=name2protein(input("Inserisci nome=\t"))
    #pdb_id=find_protein_pdb(sequence)
    #showStructurePyMOL(pdb_id=pdb_id)


    #url="https://www.rcsb.org/search?request="
    #query:dict={
    #"query": {
        #"type": "terminal",
        #"service": "seqmotif",
        #"parameters": {
        #"value": search_motif,
        #"pattern_type": "prosite",
        #"sequence_type": "protein"
        #}
    #},
    #"return_type": "polymer_entity"
    #}
    #url=url+str(query)
    #url=url.replace(" ","").replace("'",'"')
    url=_generate_url_searchmotif("ad")


    #options = Options()
    #options.binary_location = "/usr/sbin/google-chrome-stable"
    ##options.add_argument("user-data-dir=/home/matteo/cache")
    #driver = webdriver.Chrome(chrome_options = options, executable_path='/home/matteo/chromedriver_linux64/chromedriver')
    #driver.get(url)
    pdb_id=find_protein_pdb("awdawd")
    showStructurePyMOL(pdb_id=pdb_id)
    exit(5)
    driver.quit()
    #WebDriverWait wait = new WebDriverWait(webDriver, timeoutInSeconds);
    #wait.until(ExpectedConditions.visibilityOfElementLocated(By.id<locator>));
    try:
        element = WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.CLASS_NAME, "results-item"))
        )
    except NoSuchElementException:
        print("not found")

    ##driver.find_element((By.CLASS_NAME, "results-item"))
    for el in driver.find_elements(By.TAG_NAME,"h3"):
        print(re.findall(r"[A-Z\d]{4}",el.accessible_name)[0])
    #print(driver.page_source.count("results-item"))
    try:
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.CLASS_NAME, "results-item"))
        )
        elements=driver.find_elements(By.TAG_NAME,"h3")
        pdb_id=re.findall(r"[A-Z\d]{4}",elements[0].accessible_name)[0]
    except NoSuchElementException:
        print("not found")
    driver.get(_generate_url_structure(pdb_id))
    input("wait press button")

    driver.quit()



    #while True:
        #sequence=process_name(input("Inserisci nome=\t"))
        #print(sequence)
        ##name="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFY"
        #blast_results=get_alignment(sequence)
        #pickleObject(blast_results, "insuline")

    