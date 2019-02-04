import pysam 
import csv
import sys

def check_region(chrom, start, end, samfile): 
	
	reads = list(samfile.fetch(chrom, start, end))
	return len(reads) > 0


def get_allele_count(chrom, start, end, samfile, ref, alt, snp): 

	allele_dict = {"ref":0, "alt":0, "other":0}

	for read in samfile.fetch(chrom, start, end):

		reference = read.get_reference_sequence()
		positions = read.positions
		obs = read.seq

		allele, pos = get_allele(obs, positions, start-1)
		ref_allele, pos = get_allele(reference, positions, start-1)

		allele = correct_bisulfite(ref_allele, allele)
		
		cg = check_CG(reference, pos)
		
		if not cg and len(ref)==1 and len(alt)==1 and not allele=="X":
			check_allele_status(allele, ref, alt, allele_dict)
		
	return allele_dict

			
		
def check_allele_status(allele, ref, alt, allele_dict): 
	
	if allele==ref:
		allele_dict['ref'] += 1
	elif allele==alt:
		allele_dict['alt'] += 1
	else:
		allele_dict['other'] += 1


def get_allele(read, positions, start): 

	if start in positions:
		i = positions.index(start)
		return read[i], i 
	else: 
		return "X", "X"


def convert_reverse(sequence): 
	
	convert = {"T":"C", "A":"G"}

	return convert[allele]


def get_positions(file): 

	snps = []

	with open(file, "r") as f: 
		reader = csv.reader(f, delimiter="\t")
		
		for line in reader: 
			snps.append(line)

	return snps

def check_CG(reference, index): 

	if index is not "X": 
		read = reference[index:index+2]
		return read in ["CG", "GC"]
	else: 
		return True  

def correct_bisulfite(reference, allele):

	ref = reference.upper()
	if (ref=="G" and allele=="A") or ref=="C" and allele=="T": 
		# print("ref", ref)
		# print("allele", allele)
		# print(" ")
		return ref 
	else: 
		return allele 

if __name__ == "__main__":

	output=sys.argv[1]
	sam=sys.argv[2]
	pos_file=sys.argv[3]
	
	output_file = open(output, "w")
	output_writer = csv.writer(output_file, delimiter="\t")

	samfile = pysam.AlignmentFile(sam, "rb")
	positions = get_positions(pos_file)

	for pos in positions: 

		chrom, start, end = str(pos[0]), int(pos[1]), int(pos[2])
		snp = pos[3]

		if check_region(chrom, start, end, samfile):
			if "<" not in snp: 
				if not snp=="C:T" or snp=="G:A":
					ref, alt = snp.split(":")
					allele_counts = get_allele_count(chrom, start, end, samfile, ref, alt, snp)

					if len(ref) == 1 and len(alt) == 1: 
						output_writer.writerow(pos + [allele_counts["ref"]] + [allele_counts["alt"]] + [allele_counts["other"]])





















