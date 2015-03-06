from Bio import SeqIO, AlignIO
import sys
import domain_chop
#import open_reading_frame

class Align:
	'''
	Manipulations for alignments. Most functions are currently for comparing
	hamming distances between sequences in the alignment.
	
	Requires domain_chop in PhyloPreprocessing
	'''
	
	def __init__(self,alignment,format='fasta',as_seqs=False):
		self.records = None
		self.load(alignment,format,as_seqs)

	def load(self, infile,format='fasta',as_seqs=False):
		'''If as_seqs, load alignment as list of SeqRecord objects. 
		Else, load as an AlignIO object.
		
		Called automatically upon initialization, but can be 
		called again to replace current alignment'''
		if as_seqs:
			self.records = [i for i in SeqIO.parse(infile,format)]
		else:
			self.records = AlignIO.read(infile,format)
		
	def as_dict(self):
		'''Return a dictionary with IDs as keys and sequences as values'''
		return {rec.description: str(rec.seq) for rec in self.records}
		
	def _get_domains(self, cut_list):
		domains = domain_chop.get_domains(self.records,cut_list)
		return domains
		
	def print_domains(self,oufile,cut_list,format='fasta'):
		'''Write a file of the constituent domains contained in alignment.
		These will be ungapped, and be renamed: Seq1_D1, Seq1_D2 etc.'''
		domains = _get_domains(cut_list)
		SeqIO.write(domains,outfile,cut_list,format)