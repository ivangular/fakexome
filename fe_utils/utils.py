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

from fe_utils.CrossMap import *
#########################################
def translate_positions(positions, chromosome, from_assembly, to_assembly, rootname=None):

	if type(positions)==set: positions = list(positions)
	int_positions = [int(p) for p in positions]
	if from_assembly == to_assembly:
		return dict(zip(int_positions, int_positions))
	# GRCh37 and hg19 only differ for MT
	if (from_assembly.lower() in ['grch37', 'hg19']) and (to_assembly.lower() in ['grch37', 'hg19']) and (chromosome != "MT"):
		return dict(zip(int_positions, int_positions))
	if "grch" in from_assembly.lower():  from_assembly = from_assembly.lower().replace("grch", "GRCh")
	if "grch" in to_assembly.lower():  to_assembly = to_assembly.lower().replace("grch", "GRCh")

	if not rootname: rootname="{}.{}.{}.{}".format(chromosome, from_assembly, to_assembly, str(os.getpid()))
	# otherwise we'll need  tools to translate
	chain_file ="/storage/databases/liftover/{}To{}.over.chain".format(from_assembly, to_assembly.capitalize())
	if not os.path.exists(chain_file):
		print(chain_file, "not found")
		exit()

	chrom = chromosome if 'chr' in chromosome else 'chr'+chromosome
	outfile = "%s.tsv"%rootname
	outf = open (outfile,"w")
	for p in int_positions:
		outf.write("\t".join([chrom, str(p), str(p), str(p)]) + "\n")
	outf.close()

	# this is CrossMap now
	outfile_translated  =  "%s.translated.tsv"%rootname
	(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table=False)
	crossmap_bed_file(map_tree, outfile, outfile_translated)

	#read  regions back in - note that some positions might end upr untranslatable (the unmap file below)
	#others might be on entirely different chromosome - deal with that some other time
	wrongchrom = open(rootname+'.wrong_chrom', "w")
	translation = {}
	with open(outfile_translated,"r") as inf:
		for line in inf:
			if len(line.replace(" ",""))==0: continue
			f = line.rstrip().split("\t")
			if f[0]!=chrom:
				wrongchrom.write(line)
			else:
				translation[int(f[3])] = int(f[2])
	wrongchrom.close()
	# remove aux files
	os.remove(outfile)
	os.remove(outfile_translated)
	if os.stat(outfile_translated+".unmap").st_size != 0:
		print("Warning: {}.unmap in {}".format( outfile_translated, os.getcwd()))
		print("        is not empty - check this file for untranslated positions.")
	else:
		os.remove(outfile_translated+".unmap")
	if os.stat(wrongchrom.name).st_size != 0:
		print("Warning: {} in {}".format( wrongchrom.name, os.getcwd()))
		print("        is not empty - check this file for [psitions translating to antoher chromosome.")
	else:
		os.remove(wrongchrom.name)

	return translation
