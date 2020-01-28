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

from fe_utils.utils import *
from fe_utils.annotation import *

# todo store this to database

cases = [
	{
		"title_short":"cutis laxa",
		"display_title": "A wrinkly woman case",
		"case_descr_en": '''
			The patient has dry, thin, wrinkled skin, making her look much older than her actual age (33). 
			Subcutaneous fat layer thinned out. Low muscle tone. Hair, nails and teeth normal. 
			Cardiovascular evaluation by echocardiography revealed mitral valve prolapse and diffuse 
			changes in myocardium. Ophthalmologic examination demonstrated exotropia and myopic 
			astigmatism of both eyes. Intelligence  normal.
			No family history of similar conditions.
		''',
		"case_descr_hr": '''
			Pacijentica ima suhu tanku kožu sklonu borama, što ju čini starijom od njezine biološke dobi od 33 godine.
			Subkutani sloj masnih stanica stanjen. Slab mišićni tonus. Kosa, nokti i zubi normalni. Ehokardiografija 
			ukazuje na prolaps mitralnog zaliska i difuzne promjene u miokardiju. Oftalmološki nalaz: 
			egzotropija i kratkovidni astigmatizam u oba oka. Inteligenicja normalna. Obiteljska anamneza bez 
			sličnih simptoma.
		''',
		"paper_link": "https://bmcdermatol.biomedcentral.com/articles/10.1186/s12895-019-0084-6",
		"inheritance": "de novo",
		"dominance": "AD",
		"zygosity": "heterozygous",
		"cDNA": "c.2323delG",
		"variant": "7:74068661:GG:G",
		"gene": "ELN",
		"protein_effect": "frameshift",
		"structure": "",
		"kegg": ""
	},
	{
		"title_short":"Trigonocephaly",
		"display_title":"Trigonocephaly",
		"case_descr_en": '''
		 An adult young man patient with developmental delay, trigonocephaly, speech impairment, 
		 intellectual disability, and autistic features. Parents ar non-consanguineous and healthy. 
		 He also presents with various dysmorphic features, such as hypertelorism, convergent strabismus, d
		 ownslanted palpebral fissures, epicanthus,large hands and long fingers, horseshoe kidney,
		  bilateral cryptorchidism, slight deformity of the forefoot with flat-footedness in metatarsal level 
		  and clinodactyly of the fifth toe.
		''',
		"case_descr_hr": '''  
		''',
		"paper_link": "https://www.nature.com/articles/s41598-017-19109-9",
		"dominance": "AD",
		"inheritance": "de novo",
		"zygosity": "",
		"cDNA": "c.1428 + 1 G > A",
		"variant": "3:70977642:G:A",
		"gene": "FOXP1",
		"protein_effect": "splice site",
		"structure": "",
		"kegg": ""

	},

	{
		"title_short":"paraplegia",
		"display_title":"A family with paraplegia",
		"case_descr_en": '''
		Three brothers from a consanguineous family, in the age range from pre-adolescent to adult. 
		Their symptoms include severe intellectual disability, tetraspasticity with dyskinetic movements, 
		and dysmetry of the upper limbs.  Walking only short distances with a walking frame, or not at all.
		Speech is absent, but there is some language understanding. 
		''',
		"case_descr_hr": '''
		Three brothers from a consanguineous family, in the age range from pre-adolescent to adult. 
		Their symptoms include severe intellectual disability, tetraspasticity with dyskinetic movements, 
		and dysmetry of the upper limbs.  Walking only short distances with a walking frame, or not at all.
		Speech is absent, but there is some language understanding. 
		''',
		"paper_link": "http://molecularcasestudies.cshlp.org/content/3/4/a001537.full",
		"dominance": "AR",
		"inheritance": "bi-parental",
		"gene": "SPART",
		"variant": "13:36314341:G:A",
		"protein_effect": "Arg457*",
		"structure": "",
		"kegg": "",
		"mechanism_discernible": "yes" # from uniprot

	},
	{
		"title_short":"cryptorchidism",
		"display_title": "Familial bilateral cryptorchidism",
		"case_descr_en": '''
		 A family with four boys with isolated bilateral cryptorchidism. Parents are distant cousins.
		 The father and a fifth male child were unaffected.  There was no history of hyperpigmentation, 
		 bifid scrotum or hypospadias. Also, there was no history of any urinary complaints or growth failure. 
		 All the boys had undergone one-sided orchidopexy at local hospital and were progressing well with
		 the spontatnous onset of pubery.
		''',
		"case_descr_hr": '''
		 A family with four boys with isolated bilateral cryptorchidism. Parents are distant cousins.
		 The father and a fifth male child were unaffected.  There was no history of hyperpigmentation, 
		 bifid scrotum or hypospadias. Also, there was no history of any urinary complaints or growth failure. 
		 All the boys had undergone one-sided orchidopexy at local hospital and were progressing well with
		 the spontatnous onset of pubery.
		
		''',
		"paper_link": "https://jmg.bmj.com/content/56/11/727",
		"dominance": "AR",
		"inheritance": "bi-parental",
		"gene": "RXFP2",
		"variant": "13:31792798:G:A",
		"protein_effect": "Gly499Glu",
		"structure": "no",
		"kegg": "yes",
		"mechanism_discernible": "yes" # from KEGG

	},


	{
		"title_short":"neurodevelopmental disorder",
		"display_title": "A boy with neurodevelopmental disorder",
		"case_descr_en": '''
		The boy's family is an extended consanguineous family  with three affected childrenin two branches.
		The parents of each branch are first cousins. The proband, 11-year old boy 
		suffers from delayed development  and remains  nonverbal  
		and  non-ambulatory  with  profound  intellectual  disability. During neurologic examination, h
		e was hypotonic and had dysphagia, muscle weakness and atrophy. Seizures began at age 5.
		Behaviorally, he  show autism, and  has demonstrated hand  and  facial  stereotypies,  laughing spells,  
		anxiety  and  mood  changes,  body  rocking,  agitation,  sleep  problems, hyperactivity, 
		bruxism, and hand biting. 
		''',
		"case_descr_hr": '''
		The boy's family is an extended consanguineous family  with three affected childrenin two branches.
		The parents of each branch are first cousins. The proband, 11-year old boy 
		suffers from delayed development  and remains  nonverbal  
		and  non-ambulatory  with  profound  intellectual  disability. During neurologic examination, h
		e was hypotonic and had dysphagia, muscle weakness and atrophy. Seizures began at age 5.
		Behaviorally, he  show autism, and  has demonstrated hand  and  facial  stereotypies,  laughing spells,  
		anxiety  and  mood  changes,  body  rocking,  agitation,  sleep  problems, hyperactivity, 
		bruxism, and hand biting. 
	
		''',
		"paper_link": "https://www.cell.com/ajhg/fulltext/S0002-9297(19)30386-6",
		"dominance": "AR",
		"inheritance": "bi-parental",
		"gene": "NTNG2",
		"variant": "9:132198071:T:G",
		"protein_effect": "W107G",
		"structure": "3ZYG",
		"kegg": "hsa04360",
		"mechanism_discernible": "yes"
	},
	{
		"title_short":"EEEO",
		"display_title":"Early infantile epileptic encephalopathy",
		"case_descr_en": '''
		Two sisters, from healthy non-cosanguineous parents. Proband: he pregnancy was uneventful, and she was delivered s
		pontaneously with no asphyxia at 39 weeks of gestation. 
		 At 3 months of age, she experienced epileptic spasms and electroencephalography (EEG) revealed hypsarrhythmia. 
		 At 3 years 2 months, she exhibited focal impaired awareness seizures, 
		 which responded well to a combination of clonazepam (CZP) and zonisamide. At 4 years of age, myoclonic seizures 
		 appeared and occurred five times per week despite the administration of antiepileptic medications. 
		 Her psychomotor development was severely delayed. Her younger sister present simialr symptoms with 
		 epileptic spasms onset at the age of 4 months.
		''',
		"case_descr_hr": '''
		 At 3 months of age, she experienced epileptic spasms and electroencephalography (EEG) revealed hypsarrhythmia. 
		 At 3 years 2 months, she exhibited focal impaired awareness seizures, 
		 which responded well to a combination of clonazepam (CZP) and zonisamide. At 4 years of age, myoclonic seizures 
		 appeared and occurred five times per week despite the administration of antiepileptic medications. 
		 Her psychomotor development was severely delayed. Her younger sister present simialr symptoms with 
		 epileptic spasms onset at the age of 4 months.
		''',
		"paper_link": "https://www.cell.com/ajhg/fulltext/S0002-9297(18)30004-1",
		"inheritance": "bi-parental",
		"dominance": "",
		"zygosity": "compound heterozygous",
		"variant":"6:42937717:G:C;6:42938674:G:GCAGCAGGCGGCAGGAGT",
		"gene": "CNPY3",
		"protein_effect": "",
		"structure": "",
		"kegg": ""

	},

	{
		"title_short":"Niemann–Pick",
		"display_title":"A man with lipid storage disorder",
		"case_descr_en":"Our patient was diagnosed at age 33 when he presented with a 10-yr history of difficulties in "
		                "judgment, concentration, speech, and coordination. A history of transient neonatal jaundice "
		                "and splenomegaly with bone marrow biopsy suggesting a lipid storage disorder. "
		                "Symptoms progressed over >20 yr to severe ataxia and spasticity, dementia, and dysphagia "
		                "with aspiration leading to death. Brain autopsy revealed mild atrophy "
		                "of the cerebrum and cerebellum.",
		"case_descr_hr":"Our patient was diagnosed at age 33 when he presented with a 10-yr history of difficulties in "
		                "judgment, concentration, speech, and coordination. A history of transient neonatal jaundice "
		                "and splenomegaly with bone marrow biopsy suggesting a lipid storage disorder. "
		                "Symptoms progressed over >20 yr to severe ataxia and spasticity, dementia, and dysphagia "
		                "with aspiration leading to death. Brain autopsy revealed mild atrophy "
		                "of the cerebrum and cerebellum.",

		"paper_link": "http://molecularcasestudies.cshlp.org/content/2/6/a001222.full",
		"inheritance": "bi-parental",
		"dominance": "AR",
		"zygosity": "compound heterozygous",
		"gene": "NPC1",
		"variant":"18:23535479:T:C;18:23539418:C:T",
		"protein_effect": "V950M,N1156S",
		"structure": "5U73",
		"kegg": "",
		"mechanism_discernible": "yes" # yes from uniprot

	},



	{
		"title_short":"oral–facial–digital",
		"paper_link": "http://molecularcasestudies.cshlp.org/content/3/4/a001321.full",
		"dominance": "AR",
		"zygosity": "compound heterozygous",
		"gene": "",
		"protein_effect": "Met113Arg, Arg230Ter",
		"structure": "no",
		"kegg": "no"

	},


	{
		"title_short":"developmental disorder",
		"display_title":"A boy with severe developmental disorder",
		"paper_link": "https://www.cell.com/ajhg/fulltext/S0002-9297(19)30051-5",
		"case_descr_en": '''
			An 8 year old boy, child of healthy non-consanguineous parents, presents with facial dysmorphism, hypotonia, motor delay, and walking difficulty.
			Brain CT normal. Intellectual disability described as either moderate or severe. Feeding difficulty: 
			congenital pyloric stenosis. Heart: Perimembranous ventricular septal defects, double orifice mitral valve.
		''',
		"case_descr_hr": '''  
			An 8 year old boy, child of healthy non-consanguineous parents, presents with facial dysmorphism, hypotonia, motor delay, and walking difficulty.
			Brain CT normal. Intellectual disability described as either moderate or severe. Feeding difficulty: 
			congenital pyloric stenosis. Heart: Perimembranous ventricular septal defects, double orifice mitral valve.
		''',
		"dominance": "AD",
		"inheritance": "de novo",
		"gene": "CDK8",
		"cDNA": "c.88G>A",
		"variant": "13:26254729:G:A",
		"protein_effect": "Gly30Ser",
		"structure": "5xs2",
		"kegg": "no",
		"mechanism_discernible": "yes; really hard to pick out frmo the crowd; uniprot lists no disease" # from structure - proximity to ligand

	},


	{
		"title_short":"cardiomyopathy",
		"paper_link": "https://jmg.bmj.com/content/57/1/23",
		"dominance": "AR",
		"inheritance": "bi-parental",
		"gene": "SOD2",
		"protein_effect": "Gly181Val",
		"structure": "1AP5",
		"kegg": "no",
		"mechanism_discernible": "no" # maybe if close to catalytic site?

	},
	{
		"title_short": "vitreoretinopathy",
		"paper_link": "http://molecularcasestudies.cshlp.org/content/4/3/a002519.full",
		"dominance": "AR",
		"inheritance": "de novo",
		"gene": "CAPN5",
		"protein_effect": "Arg289Trp",
		"structure": "no",
		"kegg": "no",
		"mechanism_discernible": "no"

	},
	{
		"title_short":"Bainbridge-Ropers",
		"paper_link": "http://molecularcasestudies.cshlp.org/content/4/3/a002410.full",
		"dominance": "AR",
		"inheritance": "de novo",
		"gene": "ASXL3",
		"protein_effect": "R1036X",
		"structure": "",
		"kegg": "no"

	},
	{
		"title_short":"dysplasia",
		"paper_link": "http://molecularcasestudies.cshlp.org/content/4/1/a002139.full",
		"dominance": "AR",
		"gene": "WISP3",
		"protein_effect": "Ser290Leufs*12",
		"structure": "no",
		"kegg": "no",
		"mechanism_discernible": "yes" # from omim , providing it works

	},
	{
		# heritable cancer, no faimly incuded in WES
		"title_short":"cerebellar disorder ",
		"paper_link": "http://molecularcasestudies.cshlp.org/content/2/6/a001230.full",
		"dominance": "AD",
		"zygosity": "heterozygous",
		"inheritance": "paternal",
		"gene": "EGFR",
		"protein_effect": "Cys326Phe",
		"structure": "1ivo",
		"kegg": "yes",
		"mechanism_discernible": "yes"

	},
	{
		 # the paper is not open
		"title_short":"ARVC",
		"paper_link": "https://jmg.bmj.com/content/early/2020/01/10/jmedgenet-2019-106394",
		"dominance": "AD",
		"inheritance": "de novo",
		"gene": "FLNC",
		"protein_effect": "Glu2189Ter",
		"structure": "",
		"kegg": "",
		"mechanism_discernible": "yes" # yes from uniprot
	},

	{
		# hard to solve
		"title_short":"RTT",
		"display_title": "Rett syndrome‐like patient",
		"case_descr_en": '''
			A five-year old girl with developmental delay - small body mass and head circumference, 
			cold  hands and feet, central hypotonia with and  spine deformity. Possible epileptic seizure at 6 mo.
			In addition, our patient also exhibits distinct features typically seen in RTT, 
			including stereotypic hand movements, teeth grinding, limited speech, feeding difficulties, 
			breath‐holding, inappropriate outbursts of laughter, a diminished response to pain, and no facial dysmorphism. 
		''',
		"case_descr_hr": '''  
			A five-year old girl with developmental delay - small body mass and head circumference, 
			cold  hands and feet, central hypotonia with and  spine deformity. Possible epileptic seizure at 6 mo.
			In addition, our patient also exhibits distinct features typically seen in RTT, 
			including stereotypic hand movements, teeth grinding, limited speech, feeding difficulties, 
			breath‐holding, inappropriate outbursts of laughter, a diminished response to pain, and no facial dysmorphism. 
		''',
		"paper_link": "https://onlinelibrary.wiley.com/doi/full/10.1002/ccr3.2511",
		"dominance": "AR",
		"zygosity": "heterozygous",
		"inheritance": "de novo",
		"gene": "EEF1A2",
		"protein_effect": "",
		"structure": "not rally in the mutated region",
		"kegg": ""

	},
	{
		"title_short":"wooly hair",
		"paper_link": "https://www.jidonline.org/article/S0022-202X(15)34258-5/fulltext",
		"dominance": "",
		"zygosity": "",
		"gene": "",
		"protein_effect": "",
		"structure": "",
		"kegg": ""

	},

	{
		"title_short":"",
		"paper_link": "",
		"dominance": "",
		"zygosity": "",
		"gene": "",
		"protein_effect": "",
		"structure": "",
		"kegg": ""

	}
]


#########################################
def cookup_disease_variants(cursor, case):
	disease_variants = {}
	if case["inheritance"]=="de novo":
		[chrom, variant_pos, nt_from, nt_to]  = case["variant"].split(":")
		ref_nt = nt_from[0]
		genotype = "{}:{}|{}:{}".format(nt_from, nt_to, ref_nt, ref_nt)
		# TODO more interesting genotype for parents
		gt_mom = gt_dad =  "{}:{}|{}:{}".format(ref_nt, ref_nt,ref_nt, ref_nt )
		line = "\t".join([chrom,  variant_pos,  genotype,  gt_mom,  gt_dad])
		annotation_fields = annotate(cursor,  line)
		line_annotated = add_flags(annotation_fields)
		fields = line_annotated.split("\t")
		disease_variants[int(fields[1])] = fields[2:]

	elif case["inheritance"]=="bi-parental":
		variants = case["variant"].split(";")
		compound = len(variants)>1
		flip = True
		for variant in variants:
			[chrom, variant_pos, nt_from, nt_to]  = variant.split(":")
			ref_nt = nt_from[0]
			if compound:
				genotype = "{}:{}|{}:{}".format(nt_from, nt_to, ref_nt, ref_nt)
				if flip:
					gt_mom =  "{}:{}|{}:{}".format(nt_from, nt_to, ref_nt, ref_nt )
					gt_dad =  "{}:{}|{}:{}".format(ref_nt, ref_nt, ref_nt, ref_nt )
				else:
					gt_dad =  "{}:{}|{}:{}".format(nt_from, nt_to, ref_nt, ref_nt )
					gt_mom =  "{}:{}|{}:{}".format(ref_nt, ref_nt, ref_nt, ref_nt )
				flip = not flip

			else:
				genotype = "{}:{}|{}:{}".format(nt_from, nt_to, nt_from, nt_to)
				# TODO more interesting genotype for parents
				gt_mom = gt_dad =  "{}:{}|{}:{}".format(nt_from, nt_to, ref_nt, ref_nt )
			line = "\t".join([chrom,  variant_pos,  genotype,  gt_mom,  gt_dad])
			annotation_fields = annotate(cursor,  line)
			line_annotated = add_flags(annotation_fields)
			fields = line_annotated.split("\t")
			disease_variants[int(fields[1])] = fields[2:]


	return disease_variants

from random import random
#########################################
def prune (vars_per_chrom, prune_params):
	[filter, freq] = prune_params
	deletables = {}
	for chrom, vars in vars_per_chrom.items():
		deletables[chrom] = []
		for pos, var in vars.items():
			flags = int(var[-1])
			if flags&filter==filter and random()<freq:
				deletables[chrom].append(pos)

	for chrom, positions in deletables.items():
		for pos in positions:
			# print(pos, vars_per_chrom[chrom][pos])
			del vars_per_chrom[chrom][pos]

#########################################
def main():

	mysql_conf_file = "/home/ivana/.tcga_conf"
	# TODO store omim to idetifier_maps (mysql)
	omimpath = "/storage/databases/omim/mim2hgnc.tsv"
	case_vcfs_dir = "dashboard_assets/case_vcfs"
	id_resolution_dir = "dashboard_assets/case_id_resolution_tables"
	titles_dir = "dashboard_assets/case_titles"
	descr_dir = "dashboard_assets/case_descriptions"
	for dep in [mysql_conf_file,  omimpath, case_vcfs_dir, id_resolution_dir, titles_dir, descr_dir]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	hgnc2omim = read_omim(omimpath)
	db = connect_to_mysql(mysql_conf_file)
	cursor = db.cursor()

	base_name = "child"
	prune_filter = [[DE_NOVO|EXONIC, 0.75], [DE_NOVO|EXONIC, 0.75],
					[PARENT_HOMOZYGOTE|EXONIC, 0.75], [PARENT_HOMOZYGOTE|EXONIC, 0.75],
					[PARENT_HOMOZYGOTE|EXONIC, 0.75], [PARENT_HOMOZYGOTE|EXONIC, 0.75],
					[PARENT_HOMOZYGOTE|EXONIC, 0.75]]

	for case in cases[5:7]:
		idx = cases.index(case)
		label = idx +  1
		print(label, case["title_short"])

		# output title
		outfname = "{}/case{}.title.txt".format(titles_dir, str(label).zfill(3))
		with open(outfname,"w") as outf: outf.write(case["display_title"])

		# output description english
		outfname = "{}/case{}.en.txt".format(descr_dir, str(label).zfill(3))
		with open(outfname,"w") as outf: outf.write(case["case_descr_en"])
		# output description croatian
		outfname = "{}/case{}.hr.txt".format(descr_dir, str(label).zfill(3))
		with open(outfname,"w") as outf: outf.write(case["case_descr_hr"])

		infname = "{}{}.w_de_novo.vcf".format(base_name, label)
		print(label, case["title_short"], infname)
		# read in the variants file
		vars_per_chrom = get_variants_per_chrom(infname)

		# additional pruning to make the case solvable
		prune (vars_per_chrom, prune_filter[idx])

		# cook up the offending mutation
		disease_variants = cookup_disease_variants(cursor, case)

		# add to the variants
		chrom = case["variant"].split(":")[0]
		vars_per_chrom[chrom].update(disease_variants)

		# output the new file
		outfname = "{}/case{}.tsv".format(case_vcfs_dir, str(label).zfill(3))
		hdr = "# write new header"
		genes_affected = variants_output(outfname, hdr, vars_per_chrom)

		# redo the gene id's table
		outfname = "{}/case{}.ids.tsv".format(id_resolution_dir, str(label).zfill(3))
		output_affected_genes(cursor, outfname, genes_affected, hgnc2omim)



		#exit()

	return

#########################################
if __name__ == '__main__':
	main()
