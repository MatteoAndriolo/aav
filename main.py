from pathlib import Path
import pickle
import re
import logging as log
from sqlite3 import paramstyle
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio import SearchIO, SeqIO
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from Bio.PDB import PDBParser, PDBList
import nglview as nv


f_cache=""
amm=[ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def get_alignment(sequence: str):
    log.info("get_alignment: Blast job started")
    result_handle=NCBIWWW.qblast("blastp","pdb",sequence)
    blast_record=SearchIO.read(result_handle,"blast-xml")
    log.info("get_alignment: Blast result received")
    log.debug("get_alingment: Blast results\n{blast_record}\n\n"+50*"-")
    return blast_record


def replace(ch:str)->str:
    if ch in amm:
        return ch
    else:
        a=list(sorted(
                amm, key=lambda x:abs(ord(ch.upper())-ord(x))
            ))[0]
        return a

def process_name(name:str)->str:
    name=re.sub(r'\W+', '', name)
    newname=list(map(replace, [*name]))
    newname="".join(newname)
    return newname


def pickleObject(alignment:object, name:str)-> None:
    # Save object in a file 
    with open(f'{name}.pkl', 'wb') as outp:
        pickle.dump(alignment, outp, pickle.HIGHEST_PROTOCOL)

def unpickleObject(path:str=None, Path:Path=None)-> None:
    # Load Object from file
    with open(path, 'rb') as outp:
        return pickle.load(outp)
    
def downloadpdbfile(pdb_id:str):
    pdblist_manager=PDBList()
    log.debug("downloadpdbfiles: downloading")
    pdblist_manager.retrieve_pdb_file(pdb=[pdb_id])
    
    
def getStructure(pdb_id:str):
    downloadpdbfile(pdb_id)
    parser=PDBParser()
    structure = parser.get_structure(pdb_id,f"{pdb_id}.pdb")
    nv.show_biopython(structure, gui=True)




if __name__=="__main__":
    log.basicConfig(level=log.DEBUG)
    piklereceived=unpickleObject("blast_record.pkl")
    hit_id=piklereceived[1].blast_id #"pdb|6PXV|D"
    pdb_id=hit_id.split("|")[1]
    getStructure(pdb_id)
    exit(-5)


    while True:
        sequence=process_name(input("Inserisci nome=\t"))
        print(sequence)
        #name="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFY"
        blast_results=get_alignment(sequence)
        pickleObject(blast_results, "insuline")

    