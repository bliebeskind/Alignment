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
		
		
class Codon(Align):

	'''A codon alignment class'''
	
	def __init__(self,alignment,format='fasta',as_seqs=False):
		'''Populate self.records with alignment or list of sequence objects if as_seqs=True'''
		Align.__init__(self,alignment,format,as_seqs)
		self._are_seqs = as_seqs
		
	def _calc_codon_gappiness(self,codon_col):
		assert codon_col.get_alignment_length() == 3
		site_gappys = []
		for index in xrange(3):
			site = codon_col[:,index]
			site = site.replace("?","-")
			site = site.replace("X","-")
			site = site.replace("N","-")
			site_gappys.append(float(site.count("-"))/len(site))
		codon_gappiness = sum(site_gappys)/float(len(site_gappys))
		#try:
		#	assert len(set(site_gappys)) == 1
		#except AssertionError:
		#	for row in codon_col:
		#		if "-" in row.seq and row.seq.count("-") < 3:
		#			print row.id + " has incomplete codons"
		return codon_gappiness

		
	def trim(self,gap_threshold=.5):
		assert not self._are_seqs, "To use trim, must upload as alignment, not as seqs (as_seqs=False)"
		print "Amibiguity characters '?', 'N', and 'X' will be replaced with gaps"
		out_aln = None
		for i in xrange(0,self.records.get_alignment_length() - 2,3):
			codon_col = self.records[:,i:i+3]
			codon_gappiness = self._calc_codon_gappiness(codon_col)
			if codon_gappiness <= 1 - gap_threshold:
				if out_aln == None:
					out_aln = codon_col
				else:
					out_aln += codon_col
		assert out_aln.get_alignment_length() % 3 == 0, "Ending alignment not a multiple of 3!"
		return out_aln