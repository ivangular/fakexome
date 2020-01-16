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


def translate(vcfdir, orig_vcf,  orig_assembly, target_assembly):
	# name for the new file
	old_extension = ".%s.vcf" % orig_assembly
	new_extension = ".%s.vcf" % target_assembly

	if orig_vcf[-(len(old_extension)):] != old_extension:
		print("I expected the input file ({}) to end in {}".format(orig_vcf, old_extension))
		return
	print("translating", orig_vcf)
	output_vcf = orig_vcf.replace(old_extension, new_extension)

	# extract positions
	positions = {}
	with open("{}/{}".format(vcfdir, orig_vcf),"r") as inf:
		for line in inf:
			if line[:3] != "chr": continue
			fields =  line.strip().split("\t")
			if len(fields)<2: continue
			[chrom, pos] = fields[:2]
			if not chrom in positions: positions[chrom] = set()
			positions[chrom].add(pos)
	# translate
	translation = {}
	for chrom, pos in positions.items():
		translation[chrom] = translate_positions(pos, chrom, orig_assembly, target_assembly, rootname="tmpcrsmp_"+chrom)

	# write out the new file
	outf = open("{}/{}".format(vcfdir, output_vcf), "w")
	outf.write("##translated to %s by CrossMap\n" % target_assembly)
	ok   = 0
	fail = 0
	with open("{}/{}".format(vcfdir, orig_vcf),"r") as inf:
		for line in inf:
			if line[:3] != "chr":
				outf.write(line)
				continue
			fields =  line.split("\t")
			chrom = fields[0]
			pos = int(fields[1])
			if pos in translation[chrom]:
				fields[1] = str(translation[chrom][pos])
				outf.write("\t".join(fields))
				ok += 1
			else:
				fail += 1
	outf.close()
	print("ok:{}, fail:{}".format(ok, fail))




############################
def main():
	orig_assembly = "hg19"
	target_assembly = "hg38"

	vcfdir  = "/storage/sequencing/openhumans/fakexomes"
	maledir   = "{}/males".format(vcfdir)
	femaledir = "{}/females".format(vcfdir)
	for dep in [vcfdir, femaledir, maledir]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	for vcfdir in [femaledir, maledir]:
		for fnm in os.listdir(vcfdir):
			translate(vcfdir, fnm, orig_assembly, target_assembly)

#########################################
if __name__ == '__main__':
	main()
