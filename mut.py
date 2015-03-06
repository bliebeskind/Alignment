from collections import Counter
from align import Align
from Bio import SeqIO, AlignIO

class Mut(Align):

	def __init__(self,alignment,format='fasta',as_seqs=False):
		self.records = None
		self.load(alignment,format,as_seqs)
	
	def column_counter(self,gaps=True):
		'''Generator of Counter objects for amino acids in each column'''
		assert isinstance(self.records,Align.MultipleSeqAlignment), "Must load alignment with as_seq=False"
		for i in range(len(self.records[1,:])):
			if gaps:
				yield Counter(self.records[:,i])
			else:
				yield Counter(self.records[:,i].replace('-',''))
			
	def column_freqs(self,gaps=True,as_dict=True):
		'''Generator of dictionaries of amino acid frequencies for each column'''
		for i in self.column_counter(gaps):
			total = float(sum(i.values()))
			if as_dict:
				yield dict(Counter({i:j/total for i,j in i.iteritems()}))
			else:	
				yield Counter({i:j/total for i,j in i.iteritems()})
			
	def column_freqs_dict(self,gaps=True,as_dict=False):
		'''
		Return dictionary of sites mapped to dictionaries holding amino
		acid frequencies for each column
		'''
		return {i:j for i,j in enumerate(self.column_freqs(gaps,as_dict))}
			
	def all_mutant_combinations(self,gaps=True,as_dict=False):
		'''Return the number of possible sequences arising from frequency vectors
		generated by column_freqs_dict'''
		return reduce(lambda x,y:x*y,
			[len(i) for i in self.column_freqs_dict(gaps,as_dict).values()])
			
	def num_single_site_mutants(self,gaps=True,as_dict=False):
		'''Return the number of possible sequences if only single sites
		are mutated acording to column freq vectors'''
		return reduce(lambda x,y:x+y,
			[len(i) for i in self.column_freqs_dict(gaps,as_dict).values()])