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

# add some 10% de novo mutations
# make sure that they are not disease causing

# reminder: interpreting GT field
# GT represents the genotype, encoded as allele values separated by / or |.
# If the allele value is 0, means that it is equal to the reference allele
# (what is in REF field), if 1 mean that is equal to alternative (first allele listed in ALT),
#  and if 2 it is equal to the second allele listed in ALT (if it exists).
# Thus GT = 1/1 represents a SNP homozygous for the ALT allele,
# 1/0 heterozygous, and 0/0 homozygous for the reference
# The files at hand never list the second ALT, so I suspect that two alleles, neither
# one of which is reference, are given something like this:

#chr	676339	.	A	G	.	PASS	.	GT	0/1
#chr	676339	.	AAC	A	.	PASS	.	GT	1/0

import os
from random import choice

############################
def genotype_sanity(genotype, line, chromosome, ref, alt, verbose = False):
	if sum(genotype)==0: # this would mean that both alleles have reference value
		print("ref only?")
		print(line.strip())
		print(chromosome, ref, alt, genotype)
		exit()
	if len(genotype)>2: # this would mean some chimeric phenotype
		print("chimera?")
		print(line.strip())
		print(chromosome, ref, alt, genotype)
		exit()
	if verbose and max(genotype)>1:
		print("no reference variant")
		print(line.strip())
		print(chromosome, ref, alt, genotype)
		alternative = [ref] + alt.split(",")
		print(alternative)



############################
def read_parents(vcdir, vcf):
	variants = {}
	for parent in ["m", "d"]:
		vcf_inf = open("{}/{}".format(vcdir[parent], vcf[parent]), "r")
		for line in vcf_inf:
			if line[:3]!='chr': continue
			field = line.split()
			# filter?
			if field[6]!="PASS": continue
			chromosome = field[0].replace("chr","")
			position = field[1]
			ref = field[3]
			alt = field[4]
			# 0/1 or 1/0 is not ratio, it means allele1/allele2 and 0 and 1 and perhaps 2 indicates ref (0) alt1 (1) or alt2 (2)
			# thus 0/1 or 1/0 means the same
			alternative = [ref] + alt.split(",")
			# in a small percentage of cases (but comparable to the number of
			# variants on Y chromosome) one can find "0/1" or "1/0" phenotype for male on chromosome X
			# here we take that to be a mistake
			genotype = [int(i) for i in field[9].split("/")]
			if parent=="d" and chromosome in ["X","Y"] and len(genotype)>1: continue
			genotype_sanity(genotype, line, chromosome, ref, alt, verbose=False)
			if not chromosome in variants: variants[chromosome]={}
			if not position in variants[chromosome]: variants[chromosome][position] = {"m":[], "d":[]}
			# in hope of making the algebra of variants easier later,
			# store as a list of (ref,alt) pairs where (ref, ref) i.e. no change is one of the possibilities
			variants[chromosome][position][parent].extend([(ref, alternative[a]) for a in genotype])
		vcf_inf.close()
	return variants


############################
def cross_diploid(var, pos, chromosome):
	# create offspring genotype
	if len(var["m"])>0 and len(var["d"])>0:
		var["c"] = [choice(var["m"]), choice(var["d"])]
	elif len(var["m"])>0:
		ref = var["m"][0][0][:1]
		var["c"] = [choice(var["m"]), (ref, ref)]
	elif len(var["d"])>0:
		ref = var["d"][0][0][:1]
		var["c"] = [(ref,ref), choice(var["d"])]
	else:
		print("no var position {} in chrom {}?".format(pos, chromosome))
		exit()
	# did we end up with two reference alleles?
	if var["c"][0][1] ==  var["c"][0][0] and var["c"][1][0] == var["c"][1][1]:
		del var["c"]


#############
# this is where fake exome is fake - it should
# consistently be one one chromosome (we are doing this for male child)
# instead we are picking random from both

def select_from_x(var):
	if len(var["m"])>0:
		var["c"] = [choice(var["m"])]


#############
def copy_y(var):
	if len(var["d"])>0:
		var["c"] = var["d"].copy()


#############
def cross(variants, gender):
	# autosomal chromosomes
	for chromosome in [str(i+1) for i in range(22)]:
		for pos, var in variants[chromosome].items():
			cross_diploid(var, pos, chromosome)

	# special care for sex chromosomes
	if gender=="female":
		chromosome='X'
		for pos, var in variants[chromosome].items():
			cross_diploid(var, pos, chromosome)
	else:
		chromosome = 'X'
		for pos, var in variants[chromosome].items():
			select_from_x(var)
		chromosome = 'Y'
		for pos, var in variants[chromosome].items():
			copy_y(var)


############################
def allele_compact(vts):
	return "|".join([":".join(vt) for vt in vts])


############################
def fix_wt(var, parent, chromosome, gender):

	if parent in var and len(var[parent])>0: return

	ref = var["c"][0][0]
	var[parent] = [(ref, ref), (ref, ref)]
	# special cases:
	if chromosome == "X" and parent=="d":
			var[parent] = [(ref,ref)]
	if gender == "male":
		if chromosome == "Y" and parent=="m":
			var[parent] = ['-']


############################
def main():
	expected_assembly = "hg38"
	vcfdir     = "/storage/sequencing/openhumans/fakexomes"
	maledir    = "{}/males".format(vcfdir)
	femaledir  = "{}/females".format(vcfdir)
	progenydir = "{}/progeny".format(vcfdir)
	for dep in [vcfdir, maledir, femaledir, progenydir]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	# vcf files for mom and dad
	vcdir = {"m": femaledir, "d": maledir}
	vcf = {"m": choice([fnm for fnm in os.listdir(vcdir["m"]) if expected_assembly in fnm]),
			"d": choice([fnm for fnm in os.listdir(vcdir["d"]) if expected_assembly in fnm])}

	# read parent variants
	variants = read_parents(vcdir, vcf)

	# a boy or a girl?
	gender = choice(["male", "female"])
	print(gender)

	# create progeny
	cross(variants, gender)

	# add 10% mutations which are supposed to be unique for the child - 8% SNV, 1% inserts, 1% dels

	# add the disease mutation(s)

	# annotate

	# output
	# what is the most compact format I can come up with? would it work if it is all a single string?
	# freq: 2 digits first sig digit and power
	# there could be multiple genes separated by ';'
	# esi - one of exonic, splice, intronic, upstream, downstream
	# e can be followed by another: and aaFromPOSaaTo
	# chrom pos  ref alt1|alt2  mom1|mom2  pop1|pop2 freq gene:[esiud]
	# here we just output the variants, add annotation in the next step
	outname = "{}/{}.{}.vcf".format(progenydir, "test1", expected_assembly)
	outf = open(outname, "w")
	outf.write("# mom: {}   pop: {}  gender: {}\n".format(vcf["m"], vcf["d"], gender))
	outf.write("# chr  pos  gt  gt_mom  gt_dad\n")
	for chromosome in [str(i+1) for i in range(22)] + ['X','Y']:
		for pos, var in variants[chromosome].items():
			if "c" not in var: continue
			fix_wt(var, "m", chromosome, gender)  # replace wt with dummy variant (ref, ref)
			fix_wt(var, "d", chromosome, gender)

			fields = [chromosome, str(pos), allele_compact(var["c"]), allele_compact(var["m"]), allele_compact(var["d"])]
			outf.write("\t".join(fields)+"\n")
	outf.close()


	# create a separate xref table for gene ids

	return


#########################################
if __name__ == '__main__':
	main()
