'''
classes and utility functions for pipeline_metaannotate.py

'''

import os
import PipelineMetaAssemblyKit

'''
function to build call to Prodigal
'''

def runProdigal(infile,outfile,params):
    pcall = ["prodigal -i {} -o {}".format(infile,outfile.rstrip("_peptides")+"_positions")]
    pcall.append("-a {}".format(outfile))
    if params["Prodigal_c"] != 'false':
        pcall.append("-c")
    if params["Prodigal_d"] != 'false':
        pcall.append("-d")
    pcall.append("-f {}".format(params["Prodigal_f"]))
    pcall.append("-g {}".format(params["Prodigal_g"]))
    if params["Prodigal_m"] != 'false':
        pcall.append("-m")
    if params["Prodigal_n"] != 'false':
        pcall.append("-n")
    pcall.append("-p {}".format(params["Prodigal_p"]))
    if params["Prodigal_q"] != 'false':
        pcall.append("-q")
    if params["Prodigal_s"] != 'false':
        pcall.append("-s {}".format(outfile.rstrip("_peptides")+"."+params["Prodigal_s"]))
    if params["Prodigal_t"] != 'false':
        pcall.append("-t {}".format(outfile.rstrip("_peptides")+"."+params["Prodigal_t"]))
    return(" ".join(pcall))


'''
function to build call to eggnog-mapper
'''

def runEggmap(infile,outfile,params):
    #eggnog-mapper requires python2
    pcall = ["python2 {}".format(params["Eggnogmapper_eggpath"])]
    #input output
    pcall.append("-i {} -o {}".format(infile, outfile))
    #set database dir and alignment method
    pcall.append("--data_dir {}".format(params["Eggnogmapper_eggdata"]))
    pcall.append("-m {}".format(params["Eggnogmapper_method"]))
    #annotation settings
    pcall.append("--tax_scope {}".format(params["Eggnogmapper_tax_scope"]))
    pcall.append("--target_orthologs {}".format(params["Eggnogmapper_target_orthologs"]))
    pcall.append("--go_evidence {}".format(params["Eggnogmapper_go_evidence"]))
    #HMM specific settings
    if params["Eggnogmapper_method"] == "hmm":
        pcall.append("--database {}".format(params["Eggnogmapper_hmmdb"]))
        pcall.append("--dbtype {}".format(params["Eggnogmapper_dbtype"]))
        pcall.append("--qtype {}".format(params["Eggnogmapper_qtype"]))
        pcall.append("--hmm_maxhits {}".format(params["Eggnogmapper_hmm_maxhits"]))
        pcall.append("--hmm_evalue {}".format(params["Eggnogmapper_hmm_evalue"]))
        pcall.append("--hmm_score {}".format(params["Eggnogmapper_hmm_score"]))
        pcall.append("--hmm_maxseqlen {}".format(params["Eggnogmapper_hmm_maxseqlen"]))
        pcall.append("--hmm_qcov {}".format(params["Eggnogmapper_hmm_qcov"]))
        pcall.append("--hmm_Z {}".format(params["Eggnogmapper_z"]))
    #DIAMOND specific settings
    if params["Eggnogmapper_method"] == "diamond":
        if params["Eggnogmapper_dmnd_db"] != "false":
            pcall.append("--dmnd_db {}".format(params["Eggnogmapper_dmnd_db"]))
        if params["Eggnogmapper_matrix"] != "false":
            pcall.append("--matrix".format(params["Eggnogmapper_matrix"]))
        if params["Eggnogmapper_gapopen"] != "false":
            pcall.append("--gapopen".format(params["Eggnogmapper_gapopen"]))
        if params["Eggnogmapper_gapextend"] != "false":
            pcall.append("--gapextend".format(params["Eggnogmapper_gapextend"]))
    #general settings
    pcall.append("--seed_ortholog_evalue {}".format(params["Eggnogmapper_seed_ortholog_evalue"]))
    pcall.append("--seed_ortholog_score {}".format(params["Eggnogmapper_seed_ortholog_score"]))
    if params["Eggnogmapper_override"] != "false":
            pcall.append("--override")
    if params["Eggnogmapper_no_refine"] != "false":
            pcall.append("--no_refine")
    if params["Eggnogmapper_no_annot"] != "false":
            pcall.append("--no_annot")
    if params["Eggnogmapper_no_search"] != "false":
            pcall.append("--no_search")
    if params["Eggnogmapper_keep_mapping_files"] != "false":
            pcall.append("--keep_mapping_files")
    if params["Eggnogmapper_translate"] != "false":
            pcall.append("--translate")
    if params["Eggnogmapper_usemem"] != "false":
            pcall.append("--usemem")
    pcall.append("--cpu {}".format(params["Eggnogmapper_threads"]))
    return(" ".join(pcall))
