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


from random import choice, sample

from fe_utils.utils import *


def get_info(fpath):
	with open(fpath, "r") as inf:
		line = inf.readline().strip().replace(":", "")

	fields = line.split()
	# the first field is the comment sign '#'
	fields = fields[1:]
	iterator = iter(fields)
	# not sure why dict(zip(fields, fields)) does not work but it does not
	info = dict(zip(iterator, iterator))
	return info


def add_de_novo(vars_per_chrom, vars_per_chrom_donor, gender):
	for chrom, vars in vars_per_chrom_donor.items():
		sample_size = int(len(vars_per_chrom[chrom])/10)
		simulated_de_novo = {}
		for pos in vars.keys():
			if pos in vars_per_chrom[chrom]: continue
			flags = int(vars_per_chrom_donor[chrom][pos][-1])
			if flags&COMMON: continue
			simulated_de_novo[pos] =  vars_per_chrom_donor[chrom][pos]

		if sample_size<len(simulated_de_novo):
			sample_positions = sample(simulated_de_novo.keys(), sample_size)
		else:
			sample_positions = simulated_de_novo.keys()

		orig_keys =[k for k in simulated_de_novo.keys()]
		for pos in orig_keys:
			if pos not in sample_positions:
				del simulated_de_novo[pos]
				continue
			ref_nt = simulated_de_novo[pos][0][0]
			# TODO make parents' genotypes less trivial
			trivial_gt = "%s:%s|%s:%s"%(ref_nt, ref_nt, ref_nt, ref_nt)
			haplo_gt = "%s:%s"%(ref_nt, ref_nt)
			# make sure we are not homozygous in de novo variant
			child_gt = set(simulated_de_novo[pos][0].split("|"))
			if len(child_gt)==1:
				# if we see only the reference, this is not a variant
				# we should not have ended here
				if list(child_gt)[0] == trivial_gt: continue
				if chrom=="Y":
					pass
				if chrom == "X" and gender == "male":
					pass
				simulated_de_novo[pos][0] = "|".join([haplo_gt, list(child_gt)[0]])

			#mother genotype
			simulated_de_novo[pos][-3] = trivial_gt if chrom != "Y" else "-"

			#father genotype
			simulated_de_novo[pos][-2] = trivial_gt if chrom != "X" else haplo_gt

			# fix the annotation flags
			flags = int(simulated_de_novo[pos][-1])
			flags |= DE_NOVO
			flags &= ~(PARENT_HOMOZYGOTE|HOMOZYGOTE)
			simulated_de_novo[pos][-1] = flags

		#here, we are updating the original vcf
		vars_per_chrom[chrom].update(simulated_de_novo)


############################
def get_original_header(originals, base_name, child_no):
	fpath = "{}/{}{}.hg38.vcf".format(originals, base_name, child_no)
	if not os.path.exists(fpath) or os.path.getsize(fpath)==0:
		print(fpath, "not found or is empty")
		exit()
	with open(fpath, "r") as inf:
		hdr = [next(inf) for x in range(2)]
	return "".join(hdr) # we have the newlines in there, we did not strip them


def output_affected_genes(cursor, fnm, genes_affected, hgnc2omim):
	# TODO - hgnc2omim to indentifier_maps mysql
	outf = open(fnm,"w")
	for approved in genes_affected:
		if approved=="anon": continue
		uniprot = hgnc2uniprot(cursor, approved)
		[kegg_geneid, kegg_pathways] = uniprot2kegg(cursor, uniprot)
		if not kegg_pathways: kegg_pathways="x"
		omim = hgnc2omim.get(approved, "x")
		outf.write("\t".join([approved, uniprot, omim, str(kegg_geneid), kegg_pathways]) + "\n")
		# the "x" pathwys might be disease pathways which I chose not to use
		# if kegg_pathways and kegg_pathways != "x":
		# 	kegg_pathway_ids.update(set(kegg_pathways.split(";")))
	outf.close()



############################
def outer_loop(base_name, child_no, number_of_children, male_children, female_children, cursor, hgnc2omim, originals):

	print("adding de novo for", child_no)

	annotated_fnm = "{}{}.annotated.vcf".format(base_name, child_no)
	if annotated_fnm in female_children:
		random_donor = choice(list(female_children - {annotated_fnm}))
		gender = "female"
	else:
		random_donor = choice(list(male_children - {annotated_fnm}))
		gender =  "male"

	vars_per_chrom = get_variants_per_chrom(annotated_fnm)
	vars_per_chrom_donor = get_variants_per_chrom(random_donor)

	# for each chromosome, pick 1/10 of the non-related person's variants
	# make sure they are not in parents genotype and not too common and add them to variants
	# (these will play the role of 'de novo' variants
	add_de_novo(vars_per_chrom, vars_per_chrom_donor, gender)

	hdr  = "# de_novo_donor: {} \n".format(random_donor)
	hdr += get_original_header(originals, base_name, child_no)
	w_de_novo_fnm = "{}{}.w_de_novo.vcf".format(base_name, child_no)
	genes_affected  = variants_output(w_de_novo_fnm, hdr, vars_per_chrom)

	# redo the file for identifiers
	# add kegg pathway titles
	# output affected genes (with id translation)

	fnm = 	"{}{}.affected_genes.tsv".format(base_name, child_no)
	output_affected_genes(cursor, fnm, genes_affected, hgnc2omim)

	# outputting pathway names is not worth the trouble
	# typically this results in an 8K file - all pathways are 13K
	# outf = open("{}{}.kegg_pthwys.tsv".format(base_name, child_no),"w")
	# for kegg_pthwy_id in kegg_pathway_ids:
	# 	kegg_pthwy = get_kegg_pthwy_name(cursor, kegg_pthwy_id)
	# 	if kegg_pthwy=="x": continue
	# 	outf.write("\t".join([str(kegg_pthwy_id), kegg_pthwy]) + "\n")
	# outf.close()


############################
def main():
	# TODO this whole script needs reorganization
	number_of_children = 15
	base_name = "child"

	mysql_conf_file = "/home/ivana/.tcga_conf"
	vcfdir = "/storage/sequencing/openhumans/fakexomes"
	originals = "{}/progeny".format(vcfdir)
	# get omim, while we are at that
	# TODO store omim to idetifier_maps (mysql)
	omimpath = "/storage/databases/omim/mim2hgnc.tsv"
	for dep in [mysql_conf_file, vcfdir, originals, omimpath]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	hgnc2omim = read_omim(omimpath)
	db = connect_to_mysql(mysql_conf_file)
	cursor = db.cursor()

	for i in range(number_of_children):
		annotated_fnm = "{}{}.annotated.vcf".format(base_name,  i)
		if not os.path.exists(annotated_fnm):
			print(annotated_fnm, "not found")
			exit()

	# lets output  kegg pathway names once and for all
	# with open("kegg_pthwys.tsv","w") as outf:
	# 	for [kegg_pthwy_id, name] in hard_landing_search(cursor, "select * from identifier_maps.kegg_pathway_name"):
	# 		outf.write("\t".join([kegg_pthwy_id, name]) + "\n")
	# exit()

	male_children = set()
	female_children = set()

	for i in range(number_of_children):
		annotated_fnm = "{}{}.annotated.vcf".format(base_name,  i)
		cmd = "awk '$1==\"Y\"' %s | wc -l" % (annotated_fnm)
		retval = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE)
		y_count = int(retval.stdout.decode('utf-8').strip())
		if y_count>0:
			male_children.add(annotated_fnm)
		else:
			female_children.add(annotated_fnm)

	for i in range(number_of_children):
		outer_loop(base_name, i, number_of_children, male_children, female_children, cursor, hgnc2omim, originals)

	cursor.close()
	db.close()


	return

#########################################
if __name__ == '__main__':
	main()
