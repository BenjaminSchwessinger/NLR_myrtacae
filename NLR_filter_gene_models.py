



if __name__ == '__main__':
    
    import shutil
    from pybedtools import BedTool
    import pandas as pd
    import os
    from Bio import SeqIO
    import argparse
    import subprocess
    
    description = """This program parses out the gene models on the identified NLR loci.\n
It filters them for NB-ARC domain containing proteins only.\n
Works with augustus.gff3 or braker.gtf.\n
It needs 'gene' in the feature column. It needs 'transcript_id' and 'gene_id' in the featur column."""
    
    
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('anno_fn', type=str, nargs='+',
    help='Please provide a gtf or gff3 annotation file.')
    
    parser.add_argument('pfam_fn', type=str, nargs='+',
    help='Please provide pfam output file for the corresponding proteins.')
    
    parser.add_argument('fa_fn', type=str, nargs='+',
    help='Please provide NLR loci filename used for gene prediction.')
    
    parser.add_argument('pro_fn', type=str, nargs='+',
    help='Please provide protein filename of the annotated gene models.')
    
    
    augustus_gff_fn = os.path.abspath(parser.parse_args().anno_fn[0])
    pfam_fn = os.path.abspath(parser.parse_args().pfam_fn[0])
    NLR_loci_fa_fn = os.path.abspath(parser.parse_args().fa_fn[0])
    protein_fn = os.path.abspath(parser.parse_args().pro_fn[0])
    
    #print(augustus_gff_fn)
    #print(pfam_fn)
    
    prefix_list = os.path.basename(protein_fn).split('.')

    #if len(prefix_list) > 2:  
        #prefix = '.'.join(prefix_list[:-1])
    #else:
    prefix = os.path.basename(protein_fn).split('.')[0]
    suffix = os.path.basename(augustus_gff_fn).split('.')[-1]
    
    #input sorted move to setting up the tmp directory
    TMP_DIR = os.path.join(os.getcwd(), 'tmp')
    if not os.path.exists(TMP_DIR):
        os.makedirs(TMP_DIR)
    
    ##get initial chrom names from NLR file
    NLR_chrom_list = []
    with open(NLR_loci_fa_fn, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('>'):
                line = line.strip('\n')
                NLR_chrom = line[1:]
                NLR_chrom_list.append(NLR_chrom)
                
    ##split up braker.gtf or augustus.gff3 into tmp dir files one for each NLR_locus
    
    
    tmp_beds = []
    for NLR_chrom in NLR_chrom_list:
        tmp_fn = os.path.join(TMP_DIR, F"{NLR_chrom}.tmp.{suffix}")
        tmp_beds.append(tmp_fn)
        with open (tmp_fn, 'w') as tmp_out_fh:
            with open(augustus_gff_fn, 'r') as  augfh:
                for line in augfh.readlines():
                    line = line.strip('\n')
                    if line.startswith(NLR_chrom) and 'gene' == line.split('\t')[2]:
                        print(line, file = tmp_out_fh)
                        
    #check for bedfiles with overlap
    bed_with_overlap = []
    bedtool_list = []
    for tmp_gff in tmp_beds:
        try:
            bedtool_list.append(BedTool(tmp_gff).merge(c=1, o='count').to_dataframe())
        except:
            print('This file tmp_gff file raises an error', tmp_gff)
           # bedtool_list.append(float('NaN'))
    
    #check if there are any overlapping gene models
    for tmp_df in bedtool_list:
        try:
            if sum(tmp_df.name > 1) > 0:
                bed_with_overlap.append(tmp_df.chrom.unique())
                print(F"Overlapping gene models in {tmp_df.chrom.unique()}")

        except:
            print(F'Having and issue with {tmp_df.chrom.unique()}')
            
    if len(bed_with_overlap):
        print('Script needs updating, as some gene models overlap.')
    
    
    #get all protein names with NB-ARC domain
    filter_pfam_out = pfam_fn.replace('.txt', '.NB_ARC_filtered.txt')
    NLR_gene_names = []
    with open(pfam_fn, 'r') as pfam_fh:

        for line in pfam_fh.readlines():
            line = line.strip('\n')
            try:
                line = [x for x in line.split(' ') if x != '']
                if line[5] == 'PF00931.22':
                    NLR_gene_names.append(line[0])
            except:
                pass
    
    
    
    
    #generate a filtered pfam list only conatining NB-ARC domain containing proteins
    rename_dict = {}
    with open(pfam_fn, 'r') as pfam_fh:
        with open(filter_pfam_out, 'w') as f_pfam_fh:
            for line in pfam_fh.readlines():
                try:
                    line = line.strip('\n')
                    line_split = [x for x in line.split(' ') if x != '']
                    if line_split[0] in NLR_gene_names:

                        line = prefix+'_'+line
                        print(line, file = f_pfam_fh)
                        rename_dict[line_split[0]] = F"{prefix}_{line_split[0]}"
                except:
                    pass
    
    #
    print(F"Filtered the pfam annotation file for proteins that only contain \
    the NB-ARC domain in {filter_pfam_out}.\n")
    
    #number of gene models with NB-ARC domain
    print(F"This is the number of gene models with NB-ARC domains: {len(set(NLR_gene_names))}")
    
    #pull out gene models and save them
    protein_filter_fn = protein_fn.replace('.aa', '.NB_ARC_filtered.aa')
    NLR_protein_seqs = []
    for seq in SeqIO.parse(protein_fn, 'fasta'):
        if seq.id in rename_dict.keys():
            seq.id = rename_dict[seq.id]
            NLR_protein_seqs.append(seq)
    SeqIO.write(NLR_protein_seqs, protein_filter_fn, 'fasta')
    
    ##
    print(F"Pulled out filtered gene models and saved corresponding proteins \
    as {protein_filter_fn}\n")
    
    #now filter the augusts gff file and rename all the gene and transcript ids with file the prefix
    #also adjust the coordinates to fit the original genome again. This was a bit tricky based on the 1 off errors.
    augustus_gff_filter_fn = os.path.join(os.getcwd(), F'{prefix}.NLRs.{suffix}')
    with open(augustus_gff_filter_fn, 'w') as out_fh:
        with open(augustus_gff_fn) as aug_fh:
            for line in aug_fh.readlines():
                if line.startswith('#'):
                    continue
                else:
                    try:
                        line = line.strip('\n')
                        line_list = line.split('\t')

                    except:
                        print(F"Funny line {line}")
                #print(len(line_list), line)
                
                
                if any([ x for x in rename_dict.keys() if F"transcript_id \"{x}\"" in  line]): #catches old transcript id lines
                    try: 
                        for old, new in rename_dict.items():
                            if line_list[8] == F"transcript_id \"{old}\"; gene_id \"{old.split('.')[0]}\";":
                                line_list[8] = F"transcript_id \"{new}\"; gene_id \"{new.split('.')[0]}\";"
                                #line = F'{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{line_list[4]}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                line = F'{line_list[0].split(":")[0]}\t{line_list[1]}\t{line_list[2]}\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[3]) }\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[4])}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                print(line, file = out_fh)

                    except:
                        print('Issue_1')
                    #print(line, file = out_fh)
                elif line_list[2] == 'gene': #catches old gene id lines
                    #print('gene')
                    try: 
                        for old, new in rename_dict.items():
                            if line_list[8] ==  old.split('.')[0]:
                                line_list[8] = new.split('.')[0]
                                #line = F'{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{line_list[4]}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                line = F'{line_list[0].split(":")[0]}\t{line_list[1]}\t{line_list[2]}\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[3]) }\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[4])}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                                print(line, file = out_fh)
                    except:
                        print('Issue')

                elif line_list[2] == 'transcript':
                    for old, new in rename_dict.items():
                        if line_list[8] ==  old:
                            line_list[8] = new
                            #line = F'{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{line_list[4]}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                            line = F'{line_list[0].split(":")[0]}\t{line_list[1]}\t{line_list[2]}\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[3]) }\t{int(line_list[0].split(":")[1].split("-")[0])+int(line_list[4])}\t{line_list[5]}\t{line_list[6]}\t{line_list[7]}\t{line_list[8]}'
                            print(line, file = out_fh)
                            
    command = F"getAnnoFasta.pl {augustus_gff_filter_fn} --seqfile={os.path.join(os.getcwd(), F'{prefix}.fa')}"
    try:
        output = subprocess.check_output(
            command,
            stderr=subprocess.STDOUT,
            shell=True)
    except:
        print(F'{command}.\n This did not work please check.')
    
    
    #translate coding sequence
    NLR_protein_fn = augustus_gff_filter_fn.replace(suffix, 'aa')
    NLR_proteins = []
    for seq in SeqIO.parse(augustus_gff_filter_fn.replace(suffix, 'codingseq'), 'fasta'):
        seq.seq = seq.seq.translate()
        NLR_proteins.append(seq)
    SeqIO.write(NLR_proteins, NLR_protein_fn, 'fasta')
    
    print(F"Generated new translation of proteins in {NLR_protein_fn}")

    shutil.rmtree(TMP_DIR)
    
    print('All done selecting your NLR gene models.')
