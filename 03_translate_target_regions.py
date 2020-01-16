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
# the website they are talking about is here https://earray.chem.agilent.com/suredesign/index.htm
# of note: the the newer files  refer to db=hg38, however  the openhuman files refer to hg19
# the older files result in some 20K variants (smaller coverage) while the newer give ~50K
# for the purpose of the simulator smaller should do - faster download, if nothing else

import os
from fe_utils.utils import translate_positions


############################
def main():

	orig_vcf      = "/storage/sequencing/openhumans/agilent_target_regions.S30409818.hg38.bed"
	orig_assembly = "hg38"
	target_assembly  = "hg19"
	output_vcf = "/storage/sequencing/openhumans/agilent_target_regions.S30409818.hg19.bed"
	for dep in [orig_vcf]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	positions = {}
	with open(orig_vcf,"r") as inf:
		for line in inf:
			if line[:3] != "chr": continue
			[chrom, start, end] = line.strip().split()
			if not chrom in positions: positions[chrom] = []
			positions[chrom].extend([start, end])
	translation = {}
	for chrom, pos in positions.items():
		translation[chrom] =  translate_positions(pos, chrom, orig_assembly, target_assembly, rootname="tmpcrsmp_"+chrom)
		# there are some untranslated pos (the biggest damage: chr9 2356 out of 20566)
	outf = open(output_vcf, "w")
	outf.write("translated to hg19 by CrossMap\n");
	ok   = 0
	fail = 0
	with open(orig_vcf,"r") as inf:
		for line in inf:
			if line[:3] != "chr": continue
			[chrom, start, end] = line.strip().split()
			start = int(start)
			end = int(end)
			if start in translation[chrom] and end in translation[chrom]:
				ok += 1
				s = translation[chrom][start]
				e = translation[chrom][end]
				if s>e: # what's up with this?
					s = translation[chrom][end]
					e = translation[chrom][start]
				if (e-s) > 1.1*(end-start): #not to mention this
					print((e-s), (end-start))
					continue
				outf.write("%s\t%d\t%d\n"%(chrom,s,e))
			else:
				fail += 1
	outf.close()
	print("ok:{}, fail:{}".format(ok, fail))


#########################################
if __name__ == '__main__':
	main()
