import re
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio import SearchIO, SeqIO
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO

f_cache=""
amm=[ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def get_alignment(sequence: str):
    result_handle=NCBIWWW.qblast("blastp","pdb",sequence)
    blast_record=SearchIO.read(result_handle,"blast-xml")
    print(blast_record[0:10])


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

if __name__=="__main__":
    while True:
        name=process_name(input("Inserisci nome=\t"))
        print(name)
        #name="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFY"
        get_alignment(name)
    