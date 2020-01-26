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
from fe_utils.utils import *
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

# there could be multiple genes separated by ';'
# esi - one of exonic, splice, intronic, upstream, downstream
# e can be followed by another: and aaFromPOSaaTo
# chrom pos  ref alt1|alt2  mom1|mom2  pop1|pop2  gene:[esiud] freq

# TODO: reporting pseudogenes? (possible misalignment of fragments)


stable2approved_symbol = {}
stable2approved_name = {}
stable2uniprot = {}


####################################
def get_seq_region_ids(cursor):
	# coord_system with id 4 is chromosome for GRCh38
	qry = "select name, seq_region_id from seq_region where coord_system_id=4 "
	qry += "and length(name)<3"
	seq_region_id = dict(hard_landing_search(cursor,qry))
	return seq_region_id


######
def find_approved(cursor, stable_id):
	qry  = "select approved_symbol, approved_name, uniprot_ids "
	qry += "from icgc.hgnc where  ensembl_gene_id='%s'"%stable_id
	ret = error_intolerant_search(cursor, qry)
	if not ret:
		stable2approved_symbol[stable_id] = "anon"
		stable2approved_name[stable_id] =  "anon"
		stable2uniprot[stable_id] = "anon"
	else:
		stable2approved_symbol[stable_id] = ret[0][0] if ret[0][0] else "anon"
		stable2approved_name[stable_id] = ret[0][1] if ret[0][1] else "anon"
		stable2uniprot[stable_id] = ret[0][2].split(",")[0].replace(" ", "") if ret[0][2] else "anon"

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


def is_readthrough(name):
	if "readthrough" in name:
		return True
	return False


def location_within_protein_coding_gene(cursor, gene_id, variant_pos_in_gene, verbose=False):

	# find all canonical exons in this gene
	qry = "select start_in_gene, end_in_gene from gene2exon  "
	qry += "where gene_id=%d and is_canonical=1" % gene_id
	ret2 = error_intolerant_search(cursor, qry)
	if not ret2:
		return "f" # failed to annotate

	exon_intervals = sorted(ret2, key=lambda x: x[0])

	if verbose:
		print()
		print(exon_intervals)
	placed = False
	location_code = "?"
	if variant_pos_in_gene < exon_intervals[0][0]:
		if verbose:  print("   %d  <---- upstream" % variant_pos_in_gene)
		placed = True
		location_code = "u"

	exon_id = 0

	while not placed and exon_id<len(exon_intervals):
		tot_exons = len(exon_intervals)
		tot_introns = tot_exons -1
		[start, end] = exon_intervals[exon_id]
		if verbose: print(start, end)
		if not placed:
			if variant_pos_in_gene < start:
				# is this close enough to be called splice site?
				previous_splice = exon_id>0 and  (variant_pos_in_gene - exon_intervals[exon_id-1][1])<5
				this_splice = start-variant_pos_in_gene<5
				if previous_splice or this_splice:
					location_code = "s"
				else:
					if verbose:
						# this would be prev id if we are couting from 1
						print("   %d  <---- intronic  (intron %d of %d)" % (variant_pos_in_gene, exon_id, tot_introns))
					location_code = "intr"
				placed = True
				break
			elif variant_pos_in_gene <= end:
				if verbose: print("   %d  <---- exonic (exon %d of %d)" % (variant_pos_in_gene, exon_id+1, tot_exons))
				location_code = "e"
				placed = True
				break
		exon_id += 1
	if not placed:
		if verbose: print("   %d  <---- downstream" % variant_pos_in_gene)
		location_code = "d"

	return location_code


######
def genome2cds_pos(cursor, gene_id, variant_position_in_gene, strand):
	# find all canonical exons in this gene
	qry  = "select start_in_gene, end_in_gene, canon_transl_start, canon_transl_end from gene2exon "
	qry += "where gene_id=%d and is_canonical=1 " % gene_id
	ret = error_intolerant_search(cursor, qry)
	if not ret: return None  # some problem occurred

	exon_intervals = sorted(ret, key=lambda x: x[0]) # <--- important

	length = []
	reading_coding_exons = False
	# cds_start, cds_end represent the position in the exon
	# this way of doing things comes from ensembl, not from me
	cds_start_in_gene = [x[2] for x in exon_intervals if x[2]>=0][0]
	cds_end_in_gene = [x[3] for x in exon_intervals if x[3]>=0][0]
	# if cds_start_in_gene == -1 or cds_end_in_gene == -1:
	# 	print("cds start/end not properly specified, gene_id", gene_id)

	for [exon_start, exon_end, cds_start, cds_end] in  exon_intervals:
		if exon_end<cds_start_in_gene or cds_end_in_gene<exon_start:
			length.append(0)
		else:
			start =  max(exon_start, cds_start_in_gene)
			end = min (exon_end, cds_end_in_gene)
			length.append(end-start+1)


	# sanity checking
	if sum(length)%3>0:
		print("CDS length not divisible by 3, gene_id {}  length {} ".format(gene_id, sum(length)))
		print(length)
		# there are cases like that; in Ensembl they are flagged as CDS 5' or 3' incomplete
		# for en example see ENST00000431696
		# there is nothing we can do except hope that it is not the case with the canonical sequence
		# exit()
		return "err03"

	cds_position = None
	if variant_position_in_gene<cds_start_in_gene:
		cds_position="5utr"
	elif variant_position_in_gene>cds_end_in_gene:
		cds_position="3utr"
	else:
		exon_i = -1
		for i in range(len(exon_intervals)):
			[exon_start, exon_end, cds_start, cds_end] = exon_intervals[i]
			if exon_end<cds_start_in_gene: continue
			if exon_start>cds_end_in_gene: break
			if exon_start <= variant_position_in_gene <= exon_end:
				start = max(exon_start, cds_start_in_gene)
				cds_position = sum(length[:i]) + variant_position_in_gene - start
				exon_i = i

		if exon_i==-1:
			print ("if we are here, the variant should be xonic (?)")
			return "err04"

	if cds_position and strand==-1:
		if type(cds_position)==str and "utr" in cds_position:
			cds_position = "3utr" if cds_position=="5utr" else "5utr"
		else:
			cds_position = sum(length) - cds_position - 1  # we count from 0

	return cds_position


######
def find_aa_change(coding_seq, strand, cds_position, nt_from, nt_to):
	if not coding_seq  or 'UNAVAILABLE' in coding_seq: return None
	complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	if strand==-1:
		nt_from = complement[nt_from]
		nt_to = complement[nt_to]

	# ! coding sequence is already given  on the relevant strand
	codon_seq = [coding_seq[i:i+3] for i in range(0, len(coding_seq),3)]
	pep_pos   = int(cds_position/3)
	within_codon_position = cds_position%3

	codon_from = codon_seq[pep_pos]
	if codon_from[within_codon_position] != nt_from:
		# print(cds_position)
		# print(codon_seq[:3])
		# print(codon_seq[pep_pos-1:pep_pos+2])
		# print(strand, "ref codon mismatch {} in codon {} (expected {})".format(codon_from[within_codon_position], codon_from, nt_from))
		# TODO: what exactly do you make of errors in the input?
		return "err02"

	codon_to = "".join([nt_to if i==within_codon_position else codon_from[i] for i in range(3)])

	dna = Seq(codon_from, generic_dna)
	aa_from = str(dna.translate())
	dna = Seq(codon_to, generic_dna)
	aa_to = str(dna.translate())

	pep_pos += 1  # back to counting positions from 1
	change_string= "{}{}{}".format(aa_from, pep_pos, aa_to)
	return change_string


######
def protein_level_change(cursor, gene_id, canonical_transcript_id, variant_pos_in_gene, nt_from, nt_to, strand):
	# get canonical transcript sequence
	enst = transcript2stable(cursor, canonical_transcript_id)
	qry = "select sequence from icgc.ensembl_coding_seqs where transcript_id='%s'" %  enst
	ret = search_db(cursor, qry, verbose=False)
	if not ret:
		# print ("sequence not found for transcript_id='%s'" %  enst)
		return "err01"
	coding_seq = ret[0][0]

	cds_position = genome2cds_pos(cursor, gene_id, variant_pos_in_gene, strand)
	if not cds_position: return None
	if type(cds_position)==str: return cds_position

	change_string = find_aa_change(coding_seq, strand, cds_position, nt_from, nt_to)
	return change_string


######
def exonic_variant_annotation(cursor,  gene_id, canonical_transcript_id, variant_pos_in_gene, gt, strand):

	annotation = ""
	exonic_annotation = []
	prev_allele = ""
	for allele in gt.split("|"):
		if allele==prev_allele: continue
		[nt_from, nt_to] = allele.split(":")
		if nt_from==nt_to:
			exonic_annotation.append("none")
		else:
			if len(nt_from)==1 and len(nt_to)==1:
				aa_change = protein_level_change(cursor,  gene_id, canonical_transcript_id,
												variant_pos_in_gene, nt_from, nt_to, strand)
				if not aa_change: aa_change="unk"
				exonic_annotation.append(aa_change)
			elif abs(len(nt_from)-len(nt_to))%3==0:
				exonic_annotation.append("infrm")
			else:
				exonic_annotation.append("indel")
		prev_allele = allele
	if len(exonic_annotation) > 0:
		annotation += ":" + "|".join(exonic_annotation)
	return annotation


######
from math import log10

def get_gnomad_freqs(cursor, chrom, variant_pos, gt):

	prev_allele = ""
	prev_om = ""
	order_of_mag = [] # negative exponents of 10
	for allele in gt.split("|"):
		[ref, variant] = allele.split(":")
		if allele == prev_allele:
			order_of_mag.append(prev_om)
		else:
			if ref==variant:
				om = "0"
			else:
				qry = "select variant_count, total_count from gnomad.freqs_chr_%s " % chrom
				qry += "where position = %d and variant = '%s'" % (variant_pos, variant)
				ret = error_intolerant_search(cursor,qry)
				if not ret:
					om = "9"
				else:
					[variant_ct, tot_ct] = ret[0]
					if tot_ct ==0:
						om = "9"
					else:
						om = "%1.f"%round(-log10(variant_ct/tot_ct),0)
			order_of_mag.append(om)
			prev_allele = allele
			prev_om = om

	return "|".join(order_of_mag)

#####


######
def annotate(cursor, seq_region_id, line):
	gene_annotation_fields = []

	[chrom,  variant_pos,  gt,  gt_mom,  gt_dad] = line.strip().split("\t")

	variant_pos = int(variant_pos)
	# get frequency
	single_digit_freq = get_gnomad_freqs(cursor, chrom, variant_pos, gt)

	strand = 0 # strand is meaningless if we are intergenic
	# print(chrom,  pos,  gt,  gt_mom,  gt_dad)
	# flags:
	qry  = "select gene_id, stable_id, canonical_transcript_id, seq_region_start, seq_region_end, seq_region_strand  "
	qry += "from gene where seq_region_id=%d "%seq_region_id[chrom]
	qry += "and seq_region_start-100<=%d and %d<=seq_region_end+100" % (variant_pos, variant_pos)
	ret  = error_intolerant_search(cursor, qry)
	if not ret:
		gene_annotation_fields.append("intergenic")
	else:
		for [gene_id, stable_id, canonical_transcript_id, seq_region_start, seq_region_end, strand] in ret:
			# if strand>0: continue
			if stable_id not in stable2approved_symbol: # see if we have it stored by any chance
				find_approved(cursor, stable_id)
			symbol =  stable2approved_symbol[stable_id]
			name = stable2approved_name[stable_id]
			if is_readthrough(name): continue
			# print("\n=================\n", name, symbol, chrom, strand, variant_pos, gt)
			# mark RNA's  and immune (hypervariable) genes separately - see notes above
			if is_rna(symbol, name):
				gene_annotation_fields.append("%s:r" % symbol)

			elif is_immune(symbol, name):
				gene_annotation_fields.append("%s:t" % symbol)

			else: # this should be a regular protein coding gene
				variant_pos_in_gene = variant_pos - seq_region_start
				location_code = location_within_protein_coding_gene(cursor, gene_id, variant_pos_in_gene, verbose=False)
				annotation = "%s:%s:%s" % (symbol, "+" if strand>0 else "-",  location_code)
				if location_code == "e": # we are within the exon
					annotation += exonic_variant_annotation(cursor, gene_id, canonical_transcript_id, variant_pos_in_gene, gt, strand)
				gene_annotation_fields.append(annotation)

		not_anon = [annot for annot in gene_annotation_fields if "anon:" not in annot]
		if len(not_anon)>0: gene_annotation_fields = not_anon
		not_rna = [annot for annot in gene_annotation_fields if ":r" not in annot]
		if len(not_rna)>0: gene_annotation_fields = not_rna

	# annotation_string = ";".join(gene_annotation_fields)
	# for the purposes of 2020 class, we pretend there are no cases where twp genes are affected
	# annotation_string = ";".join(gene_annotation_fields)
	annotation_string = gene_annotation_fields[0]
	return [str(x) for x in [chrom, variant_pos, gt, single_digit_freq, annotation_string, gt_mom, gt_dad]]


############################
# flags
HOMOZYGOTE   =  1
COMMON       =  2
EXONIC       =  4
SILENT       =  8
DE_NOVO      = 16
PARENT_HOMOZYGOTE = 32

import re
mut_annot_pattern = re.compile('.*e\:([A-Z])(\d+)([A-Z]).*')
############################
def add_flags(annotation_fields):
	[chrom, variant_pos, gt, single_digit_freq, annotation_string, gt_mom, gt_dad] = annotation_fields
	flags = 0

	# homozygote?
	alleles = gt.split("|")
	if len(set(alleles)) == 1:
		flags += HOMOZYGOTE

	# common?
	rare_orders_of_mag = [int(x) for x in single_digit_freq.split("|") if int(x)>2]
	if len(rare_orders_of_mag)==0:
		flags += COMMON

	# exonic? we will include splice positions here
	if ":e" in annotation_string or ":s" in annotation_string:
		flags += EXONIC

	# silent?
	if ":e" in annotation_string:
		# flag as silent if each time e appears, the
		# substitution on the protein level is silent
		outcomes = []
		for ann in  re.split('[|;]', annotation_string):
			match = re.match(mut_annot_pattern, ann)
			if match:
				outcomes.append(match.group(1)==match.group(3))
		if outcomes and len([x for x in outcomes if not x])==0:
			flags += SILENT

	# exists in parent
	child_alleles = set(gt.split("|"))
	parent_alleles = set(gt_mom.split("|")) | set(gt_dad.split("|"))
	child_alleles.difference(parent_alleles)
	if len(child_alleles.difference(parent_alleles))>0:
		flags+=DE_NOVO

	for allele in child_alleles:
		[ref,var] = allele.split(":")
		if ref==var: continue
		if gt_mom=="{}|{}".format(allele, allele):
			flags += PARENT_HOMOZYGOTE
			break
		if gt_dad=="{}|{}".format(allele, allele):
			flags += PARENT_HOMOZYGOTE
			break


	return "\t".join(annotation_fields+[str(flags)])


############################
def outer_loop(base_name):

	fnm = "{}.hg38.vcf".format(base_name)
	print("annotating", fnm)
	time0 = time()

	expected_assembly = "hg38"
	if not expected_assembly in fnm:
		print("The expected assembly %s does not appear in the name '%s'." % (expected_assembly, fnm))
		print("Note that, as of this writing, it is the version used by Ensembl.")
		exit()
	vcfdir = "/storage/sequencing/openhumans/fakexomes/progeny"
	fpath = "{}/{}".format(vcfdir, fnm)
	mysql_conf_file = "/home/ivana/.tcga_conf"
	# get omim, while we are at that
	omimpath = "/storage/databases/omim/mim2hgnc.tsv"
	for dep in [vcfdir, fpath, mysql_conf_file, omimpath]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	hgnc2omim = read_omim(omimpath)

	db = connect_to_mysql(mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
	# coord_system with id 4 is chromosome for GRCh38
	seq_region_id = get_seq_region_ids(cursor)


	inf = open(fpath,"r")
	outf = open("{}.annotated.vcf".format(base_name),"w")
	line_ct = 0
	for line in inf:
		if line[0]=='#': continue
		annotation_fields = annotate(cursor, seq_region_id, line)
		line_annotated = add_flags(annotation_fields)
		#print(line_annotated, "       {0:b}".format(int(line_annotated.split("\t")[-1])))
		outf.write(line_annotated + "\n")
		line_ct += 1
	inf.close()
	outf.close()

	# output affected genes (with id translation
	outf = open("{}.affected_genes.tsv".format(base_name),"w")
	for stable, approved in stable2approved_symbol.items():
		if approved=="anon": continue
		uniprot = stable2uniprot[stable]
		omim = hgnc2omim.get(approved, "x")
		outf.write("\t".join([approved, uniprot, omim]) + "\n")
	outf.close()

	cursor.close()
	db.close()

	print(fnm, "done, %2.1f min" % ((time()-time0)/60) )
	
	return

############################
def main():
	for i in range(15):
		outer_loop("child%d"%i)
	return

#########################################
if __name__ == '__main__':
	main()
