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

############################
def main():
	vcfdir  = "/storage/sequencing/openhumans/fakexomes"
	maledir   = "{}/males".format(vcfdir)
	femaledir = "{}/females".format(vcfdir)
	for dep in [vcfdir, maledir, femaledir]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	return


#########################################
if __name__ == '__main__':
	main()
