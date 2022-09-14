from pathlib import Path
import pickle
import re
import logging as log
from sqlite3 import paramstyle
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
from Bio.PDB import PDBParser, PDBList
import nglview as nv
from pymol import cmd



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
    '''
    #cmd.fragment('ala')
    cmd.fetch(pdb_id)
    cmd.spectrum()
    cmd.zoom()
    cmd.png('tmp/test.png', 1920, 1080)


if __name__=="__main__":
    #piklereceived=unpickleObject("blast_record.pkl")
    #hit_id=piklereceived[1].blast_id #"pdb|6PXV|D"
    #pdb_id=hit_id.split("|")[1]
    #getStructure(pdb_id)
    blast_record=get_alignment("MARCGCGSTANQG")
    pickleObject(blast_record, "backup/blastrecord")
    pdb_id=blastrecors2pdbid(blast_record)
    
    showStructurePyMOL(pdb_id)


    #while True:
        #sequence=process_name(input("Inserisci nome=\t"))
        #print(sequence)
        ##name="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFY"
        #blast_results=get_alignment(sequence)
        #pickleObject(blast_results, "insuline")

    