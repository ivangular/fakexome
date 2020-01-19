#! /usr/bin/python3
#
# This source code is part of fakexome, a simulated exome generator.
#
# Fakexome is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fakexome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
#
# Contact: ivana.mihalek@gmail.com, ivana.mihalek@uniri.hr

import os


from fe_utils.ensembl import *
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

# there could be multiple genes separated by ';'
# esi - one of exonic, splice, intronic, upstream, downstream
# e can be followed by another: and aaFromPOSaaTo
# chrom pos  ref alt1|alt2  mom1|mom2  pop1|pop2  gene:[esiud] freq

# TODO: reporting pseudogenes? (possible misalignment of fragments)


stable2approved_symbol = {}
stable2approved_name = {}

####################################
def get_seq_region_ids(cursor):
	# coord_system with id 4 is chromosome for GRCh38
	qry = "select name, seq_region_id from seq_region where coord_system_id=4 "
	qry += "and length(name)<3"
	seq_region_id = dict(hard_landing_search(cursor,qry))
	return seq_region_id


######
def find_approved(cursor, stable_id):
	ret = error_intolerant_search(cursor, "select approved_symbol, approved_name from icgc.hgnc where  ensembl_gene_id='%s'"%stable_id)
	stable2approved_symbol[stable_id] = (ret[0][0] if ret else "anon")
	stable2approved_name[stable_id] = (ret[0][1] if ret else "anon")


#####
# protein coding genes do not have  "-" in their name
# these might be anti-sense RNA (AQP4-AS1, ARAP1-AS2)
# intronic transcript (ARHGAP26-IT1, BACH1-IT2, BACH1-IT3 )
# divergent transcript (BAIAP2-DT)
# overlapping transcript (ST7-OT4)
# micro RNA, SNO RNA (MIR4435-2,  SNORD114-7), small NF90 (ILF3) associated RNA (SNAR)
# long-noncoding RNA (LINC)
# RNA (RNVU1-4         | RNA, variant U1 small nuclear 4) and
# mitochondrially encoded tRNA ( MT-TY           | mitochondrially encoded tRNA tyrosine)# readthrough (BOLA2-SMG1P6)
# endogenous retrovirus (ERVK-9)
# vault RNA's, transfer RNA
# pseudogenes (RNU7-192P)
# variables from the immune system (TRBV25-1, T cell receptor beta variable 25-1)
# IGHD6-13        | immunoglobulin heavy diversity 6-13


# unfortunately there are some dashed names that are legit:
# KRTAP1-3, keratin associated protein 1-3
# NKX6-2, NK6 homeobox 2
def is_rna (symbol, name):
	for rna in ["-AS", "-DT",  "-IT", "-OT"]:
		if rna in symbol:
			return True
	for rna in ["MIR", "SNOR", "SNAR", "LINC", "RNV", "RNU", "VTRNA"]:
		if symbol[:len(rna)]==rna:
			return True
	if "overlapping transcript" in name:
		return True
	if "mitochondrially encoded tRNA" in name:
		return True

	return False


def is_immune(symbol, name):
	if symbol[:4]=="HLA-":
		return True
	if "T cell receptor beta joining" in name:
		return True
	return False


def location_within_protein_coding_gene(cursor, gene_id, variant_pos_in_gene):

	# find all canonical exons in this gene
	qry = "select start_in_gene, end_in_gene from gene2exon  "
	qry += "where gene_id=%d and is_canonical=1" % gene_id
	ret2 = error_intolerant_search(cursor, qry)
	if not ret2:
		return "f" # failed to annotate

	exon_intervals = sorted(ret2, key=lambda x: x[0])


	placed = False
	location_code = "?"
	if variant_pos_in_gene < exon_intervals[0][0]:
		# print("   %d  <---- upstream" % variant_pos_in_gene)
		placed = True
		location_code = "u"

	exon_id = 0
	while not placed and exon_id<len(exon_intervals):
		[start, end] = exon_intervals[exon_id]
		#print(start, end)
		if not placed:
			if variant_pos_in_gene < start:
				# is this close enough to be called splice site?
				previous_splice = exon_id>0 and  (variant_pos_in_gene - exon_intervals[exon_id-1][1])<5
				this_splice = start-variant_pos_in_gene<5
				if previous_splice or this_splice:
					location_code = "s"
				else:
					#print("   %d  <---- intronic" % variant_pos_in_gene)
					location_code = "i"
				placed = True
				break
			elif variant_pos_in_gene <= end:
				#print("   %d  <---- exonic" % variant_pos_in_gene)
				location_code = "e"
				placed = True
				break
		exon_id += 1
	if not placed:
		#print("   %d  <---- downstream" % variant_pos_in_gene)
		location_code = "d"

	return location_code


######
def genome2cds_pos(cursor, gene_id, variant_position_in_gene):
	# find all canonical exons in this gene
	qry  = "select start_in_gene, end_in_gene, canon_transl_start, canon_transl_end from gene2exon "
	qry += "where gene_id=%d and is_canonical=1 " % gene_id
	ret = error_intolerant_search(cursor, qry)
	if not ret: return None  # some problem occured

	exon_intervals = sorted(ret, key=lambda x: x[0]) # <--- important

	length = []
	reading_coding_exons = False
	# cds_start, cds_end represent the position in the exon
	# this way of doing things comes from ensembl, not from me
	for [exon_start, exon_end, cds_start, cds_end] in  exon_intervals:
		if cds_start>0:
			length.append(exon_end - exon_start + 1 - (cds_start - 1))
			reading_coding_exons = True
			continue
		if cds_end>0:
			length.append(cds_end)
			reading_coding_exons = False
			continue
		if reading_coding_exons:
			length.append(exon_end-exon_start+1)
		else:
			length.append(0)

	# sanity checking
	if sum(length)%3>0:
		print("CDS length not divisible by 3, gene_id", gene_id)
		# there are cases like that; in ENsmebl they are flagged as CDS 5' or 3' incomplete
		# for en example see ENST00000431696
		# there is nothing we can do except hope that it is not the case with the canonical sequence
		return None

	cds_position = None
	for i in  range(len(exon_intervals)):
		[exon_start, exon_end, cds_start, cds_end] = exon_intervals[i]
		start = exon_start if cds_start==0 else  exon_start + cds_start - 1
		end   = exon_end if cds_end==0 else exon_start + cds_end - 1
		if start<=variant_position_in_gene<=end:
			cds_position = sum(length[:i]) + variant_position_in_gene - start

	return cds_position


######
def find_aa_change(coding_seq, strand, cds_position, nt_from, nt_to):
	if not coding_seq  or 'UNAVAILABLE' in coding_seq: return None
	complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	if strand=="-1":
		nt_from = complement[nt_from]
		nt_to = complement[nt_to]

	# ! coding sequence is already given  on the relevant strand
	codon_seq = [coding_seq[i:i+3] for i in range(0, len(coding_seq),3)]
	pep_pos   = int(cds_position/3)
	within_codon_position = cds_position%3


	codon_from = codon_seq[pep_pos]
	if codon_from[within_codon_position] != nt_from:
		print(cds_position)
		print(codon_seq[:3])
		print(codon_seq[pep_pos-1:pep_pos+2])
		print(strand, "ref codon mismatch", codon_from[within_codon_position],  nt_from)
		exit()

	codon_to = "".join([nt_to if i==within_codon_position else codon_from[i] for i in range(3)])

	dna = Seq(codon_from, generic_dna)
	aa_from = str(dna.translate())
	dna = Seq(codon_to, generic_dna)
	aa_to = str(dna.translate())

	pep_pos += 1  # back to counting positions from 1
	change_string= "{}{}{}".format(aa_from, pep_pos, aa_to)
	print(" **** ", change_string)
	return change_string


######
def protein_level_change(cursor, gene_id, canonical_transcript_id, variant_pos_in_gene, nt_from, nt_to, strand):
	# get canonical transcript sequence
	enst = transcript2stable(cursor, canonical_transcript_id)
	qry = "select sequence from icgc.ensembl_coding_seqs where transcript_id='%s'" %  enst
	ret = search_db(cursor, qry, verbose=False)
	if not ret:
		print ("sequence not found for transcript_id='%s'" %  enst)
		return None
	coding_seq = ret[0][0]

	cds_position = genome2cds_pos(cursor, gene_id, variant_pos_in_gene)
	if not cds_position: return None
	change_string = find_aa_change(coding_seq, strand, cds_position, nt_from, nt_to)

	return change_string


######
def annotate(cursor, seq_region_id, line):
	if line[0]=='#': return ""

	gene_annotation_fields = []

	[chrom,  variant_pos,  gt,  gt_mom,  gt_dad] = line.strip().split("\t")
	variant_pos = int(variant_pos)
	# print(chrom,  pos,  gt,  gt_mom,  gt_dad)
	qry  = "select gene_id, stable_id, canonical_transcript_id, seq_region_start, seq_region_end, seq_region_strand  "
	qry += "from gene where seq_region_id=%d "%seq_region_id[chrom]
	qry += "and seq_region_start-100<=%d and %d<=seq_region_end+100" % (variant_pos, variant_pos)
	ret  = error_intolerant_search(cursor, qry)
	if not ret:
		gene_annotation_fields.append("intergenic")
	else:
		for [gene_id, stable_id, canonical_transcript_id, seq_region_start, seq_region_end, strand] in ret:
			if stable_id not in stable2approved_symbol: # see if we have it stored by any chance
				find_approved(cursor, stable_id)
			# mark RNA's  and immune (hypervariable) genes separately - see notes above
			symbol =  stable2approved_symbol[stable_id]
			name = stable2approved_name[stable_id]
			if is_rna(symbol, name):
				gene_annotation_fields.append("%s:r" % (stable2approved_symbol[stable_id]))
				continue
			if is_immune(symbol, name):
				gene_annotation_fields.append("%s:t" % (stable2approved_symbol[stable_id]))
				continue
			variant_pos_in_gene = variant_pos - seq_region_start
			location_code = location_within_protein_coding_gene(cursor, gene_id, variant_pos_in_gene)
			annotation = "%s:%s:%s" % (stable2approved_symbol[stable_id], stable_id,  location_code)
			if location_code == "e": # we are within the exon
				protein_annotation = []
				for allele in gt.split("|"):
					[nt_from, nt_to] = allele.split(":")
					if nt_from==nt_to:
						protein_annotation.append("none")
					else:
						if len(nt_from)==1 and len(nt_to)==1:
							aa_change = protein_level_change(cursor,  gene_id, canonical_transcript_id,
															variant_pos_in_gene, nt_from, nt_to, strand)
							if not aa_change: aa_change="unk"
							protein_annotation.append(aa_change)
						elif abs(len(nt_from)-len(nt_to))%3==0:
							protein_annotation.append("infrm")
						else:
							protein_annotation.append("indel")
				if len(protein_annotation) > 0:
					annotation += ":" + "|".join(protein_annotation)
			gene_annotation_fields.append("%s:%s:%s" % (stable2approved_symbol[stable_id], stable_id,  location_code))

		not_anon = [annot for annot in gene_annotation_fields if "anon:" not in annot]
		if len(not_anon)>0: gene_annotation_fields = not_anon

	print(" ==   %s  %d  %s   %s" % (chrom, variant_pos, gt, ";".join(gene_annotation_fields)))
	if ":e" in ";".join(gene_annotation_fields): exit()


############################
def main():

	fnm = "test1.hg38.vcf"

	expected_assembly = "hg38"
	if not expected_assembly in fnm:
		print("The expected assembly %s does not appear in the name '%s'." % (expected_assembly, fnm))
		print("Note that, as of this writing, it is the version used by Ensembl.")
		exit()
	vcfdir = "/storage/sequencing/openhumans/fakexomes/progeny"
	fpath = "{}/{}".format(vcfdir, fnm)
	mysql_conf_file = "/home/ivana/.tcga_conf"
	for dep in [vcfdir, fpath, mysql_conf_file]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	db = connect_to_mysql(mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
	# coord_system with id 4 is chromosome for GRCh38
	seq_region_id = get_seq_region_ids(cursor)

	# cursor, gene_id, variant_pos_in_gene
	# genome2cds_pos(cursor, 397958, 3411794-3069168)
	# cursor,  gene_id, canonical_transcript_id, variant_position_in_gene
	# protein_level_change(cursor, 397958, 1397175, 3411794-3069168)
	# exit()


	inf = open(fpath,"r")
	for line in inf:
		line_annotated = annotate(cursor, seq_region_id, line)
		# print(line_annotated)
	inf.close()

	cursor.close()
	db.close()
	return


#########################################
if __name__ == '__main__':
	main()
