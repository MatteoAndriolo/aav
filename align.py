from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO
protein_seq=SeqIO.read("insuline.fasta","fasta")
print(protein_seq.seq)
result_handle=NCBIWWW.qblast("blastp","pdb",protein_seq.seq)
print(result_handle)
blast_record=SearchIO.read(result_handle,"blast-xml")
print(blast_record[0:10])

