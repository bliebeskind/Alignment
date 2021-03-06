from align import Align
from Bio import SeqIO, AlignIO
import Bio.Align

class Dist(Align):

	def __init__(self,alignment,format='fasta',as_seqs=False):
		self.records = None
		self.load(alignment,format,as_seqs)
	
	def num_pairwise(self):
		'''Return the number of possible non-self pairwise comparisons 
		for the alignment'''
		func = lambda x: sum([x-i for i in range(1,x)])
		return func(len(self.records))
		
	# Adapted from Siavash's function
	def hamming(self,k1,k2,no_gaps=False):
		''' Returns the hamming distance between two sequences. If no_gaps,
		only ungapped sites are counted'''
		s1, s2 = k1.seq, k2.seq
		if no_gaps:
			return sum(
				map( lambda x: 1 if x[0] != x[1] and '-' not in x else 0,
					zip(s1,s2) ) )
		else:
			return sum(
				map( lambda x: 1 if x[0]!=x[1] else 0,
					zip(s1,s2) ) )
		
	def hamming_dists(self,no_gaps=False):
		'''
		Hamming distances for all non-self pairwise comparisons. 
		Returns a generator yielding three part tuples:
		(SeqRecord1, SeqRecord2, hamming_dist)
		
		If no_gaps = True, only ungapped sites will be compared.
		'''
		for i,seq1 in enumerate(self.records):
			for seq2 in self.records[i+1:]:
				yield seq1, seq2, self.hamming(seq1,seq2,no_gaps)
		
	def similar_pairs(self,threshold=0,no_gaps=True):
		'''
		Find pairs whose hamming distance is below a threshold. Returns a
		generator yielding a tuple of sequence id pairs.
		
		If no_gaps = True, only ungapped sites will be compared
		'''
		num_comparisons = 0
		for seq1,seq2,dist in self.hamming_dists(no_gaps):
			num_comparisons +=1
			if num_comparisons % 100 == 0:
				print "%i pairwise comparisons completed" % num_comparisons
			if dist <= threshold:
				yield seq1.id,seq2.id
				
	def dist_info(self,threshold=0, no_gaps=False):
		'''
		Return a list of all non-self hamming distances, and list of pairs with
		hamming distance below threshold.
		
		If no_gaps = True, only ungapped sites will be compared
		'''
		dists, pairs = [], []
		for seq1,seq2,dist in self.hamming_dists(no_gaps):
			if dist <= threshold:
				pairs.append((seq1.id,seq2.id))
			dists.append(dist)
		return dists,pairs
		
	def find_shorts_hamming(self,threshold=0,no_gaps=False):
		'''Returns a list of SeqRecord ids which are not the longest 
		sequence in a group having a hamming similarity below threshold.'''
		short_seqs = []
		num_comparisons = 0
		for seq1,seq2,dist in self.hamming_dists(no_gaps):
			num_comparisons +=1
			if num_comparisons % 100 == 0:
				sys.stderr.write("%i pairwise comparisons completed\n" % num_comparisons)
			if dist <= threshold:
				#print "%s and %s below threshold" % (seq1.id, seq2.id)
				func = lambda x,y: x.id if \
					len(str(x.seq).replace("-",'')) <= \
					len(str(y.seq).replace("-",'')) else y.id
				shorter = func(seq1,seq2)
				if shorter in short_seqs: 	# shorter found previously
					continue
				else:
					short_seqs.append(shorter)
		return short_seqs

	def trim_shorts_hamming(self,threshold=0,no_gaps=False):
		'''
		Remove sequences which are not the longest sequence in a group 
		having a hamming similarity below threshold. Returns a generator of
		longest sequences for each such group, and sequences which are not in
		these groups.
		'''
		sys.stderr.write("Making all-by-all comparisons\n")
		shorts = self.find_shorts_hamming(threshold,no_gaps)
		for i in self.records:
			if i.id not in shorts:
				yield i
			else:
				print "%s has been trimmed" % (i.id)