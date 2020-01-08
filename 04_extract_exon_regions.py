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
#

# use agilent target regions
# https://kb.10xgenomics.com/hc/en-us/articles/115004150923-Where-can-I-find-the-Agilent-Target-BED-files-

# /storage/sequencing/openhumans/agilent_target_regions.S30409818.bed --> refers to hg19, as do PGP vcf files

import os
from intervaltree import IntervalTree

#############################
# "hg38 is a corrected and improved version of hg19", let's hope that's true
# http://tools.thermofisher.com/content/sfs/brochures/human-genome-38-faq.pdf
def read_regions(bedfile):
	regions = {}
	assembly_ok = False
	with open(bedfile,"r") as inf:
		for line in inf:
			fields = line.strip().split("\t")
			if fields[0][:3]!='chr':
				if 'hg19' in line or 'hg38' in line:
					assembly_ok = True
				continue
			[char,start,end] = fields
			if char not in regions: regions[char] = []
			# +1 bcs  intervaltree assumes semiopen intervals
			regions[char].append((int(start),int(end)+1))

	if not assembly_ok:
		print('no line confirming hg19 found in', bedfile)
		exit()
	return regions


#############################
def trees_from_intervals(regions):
	# it looks like intervaltree assumes semiopen intervals
	interval_tree = {}
	for chr, regs in regions.items():
		interval_tree[chr] = IntervalTree.from_tuples(regs)
	return interval_tree


#############################
def extract_exome(indir, outdir, fnm, interval_tree, verbose=True):
	# checking integer in region: 10000 <= number <= 30000
	inf = open("{}/{}".format(indir, fnm), "r")
	# check we have the xpected assembly
	assembly_ok = False
	for line in inf:
		if line[:2]=='##': # comment/meta info line
			if 'hg19' in line or 'hg38' in line:
				assembly_ok = True
				break
	if not assembly_ok:
		print('no line confirming hg19 found in', "{}/{}".format(indir, fnm))

	outf = open("{}/{}".format(outdir, fnm.replace(".vcf", ".fakexome.vcf")), "w")

	inf.seek(0, 0)  # rewind
	total = {}
	covered = {}
	for line in inf:
		# fields: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
		fields = line.strip().split("\t")
		if fields[0][:3] != 'chr':
			outf.write(line)
			continue
		chrom = fields[0]
		pos = int(fields[1])
		if verbose:
			if chrom not in covered: covered[chrom] = 0
			if chrom not in total: total[chrom] = 0
			total[chrom] +=1
		if interval_tree[chrom].overlaps(pos):
			outf.write(line)
			if verbose: covered[chrom] += 1
	inf.close()
	outf.close()

	if verbose:
		total_total = 0
		total_cvd = 0
		print('==========================')
		print(fnm, outf.name)
		for chrom, cvd in covered.items():
			print(chrom, total[chrom], covered[chrom])
			total_total += total[chrom]
			total_cvd += covered[chrom]
		print(total_total, total_cvd)
		print()


############################
def main():
	vcfdir  = "/storage/sequencing/openhumans"
	bedfile = "{}/agilent_target_regions.S30409818.bed".format(vcfdir)
	orig = "{}/orig".format(vcfdir)
	for dep in [vcfdir, bedfile, orig]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	fake_exome_dir = "{}/fakexomes".format(vcfdir)
	if not os.path.exists(fake_exome_dir):
		os.mkdir(fake_exome_dir)
	regions = read_regions(bedfile)

	interval_tree = trees_from_intervals(regions)
	# for chr, tree in interval_tree.items():
	# 	print(chr)
	# 	# for agilent_target_regions.S30409818.bed as example
	# 	# middle of the interval found ok
	# 	# middle between 2 interval found ok
	# 	print(len(tree), tree.begin(), tree.end(), tree.overlaps(249208754), tree.overlaps(249145043))
	# 	# the intervals are hlaf open (note +1 above)
	# 	print(tree.overlaps(249210798), tree.overlaps(249212564))

	for (dirpath, dirnames, filenames) in os.walk(orig):
		for fnm in filter(lambda f:  f[-4:]=='.vcf', filenames):
			extract_exome(orig, fake_exome_dir, fnm, interval_tree)

	return


#########################################
if __name__ == '__main__':
	main()
