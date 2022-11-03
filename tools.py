#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 14:02:56 2020

@author: christospapadopoulos
"""
import sys
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio import SeqIO
from Bio.Blast import NCBIXML

def read_multiFASTA(fasta_file):
    dico = {}
    with open(fasta_file,'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                name = str(line.split()[0])[1:]
                dico[name] = ''
            elif line == '\n':
                continue
            else:
                seq = line.strip()
                dico[name] = dico[name] + seq
    return(dico) 


def read_anc_list(ancestors_list):
    dico = {}
    with open(ancestors_list,'r') as f:
        for line in f:
            dico[line.split()[0]] = line.split()[-1] 
    return dico
            


def cut_into_peaces(sequence):
    dico_aa = {}
    dico_nt = {}
    for i in range(3):
        frame  = i+1
        my_seq = sequence[i:len(sequence) - len(sequence[i:len(sequence)])%3]
        codons = []
        for n in range(0, len(my_seq), 3):
            codons.append(my_seq[n:n + 3])
        #codons = [my_seq[n:n + 3] for n in range(0, len(my_seq), 3)]
        
        
        seq_tmp = ''
        my_seqs_cut = []
        for x,codon in enumerate(codons):
            if codon  not in ["TAA","TAG","TGA"] and x != len(codons)-1:
                seq_tmp = seq_tmp + codon
            elif codon in ["TAA","TAG","TGA"]:
                seq_tmp = seq_tmp + codon
                my_seqs_cut.append(seq_tmp)
                seq_tmp = ''
            else:
                seq_tmp = seq_tmp + codon
                my_seqs_cut.append(seq_tmp)
                seq_tmp = ''
                
        for j,peace in enumerate(my_seqs_cut):
            if len(peace) > 63:
                dico_aa['frame_'+str(frame)+'_part_'+str(j+1)] = str(Seq(peace).translate()).replace("*","")
                dico_nt['frame_'+str(frame)+'_part_'+str(j+1)] = peace
                
    return(dico_aa , dico_nt)
    
def decide_NandC_ter(Nter_tmp,Cter_tmp,Dec_dico,hsp_tmp,name):
    if hsp_tmp.query_start < Nter_tmp:
        Nter_tmp = hsp_tmp.query_start
        Dec_dico['Nter'] = name
    if hsp_tmp.query_end > Cter_tmp:
        Cter_tmp = hsp_tmp.query_end
        Dec_dico['Cter'] = name
    return(Nter_tmp,Cter_tmp,Dec_dico)
    
def read_blast(blast_file,frag_HCA):
    dico = {}
    
    Scer_vs_Mega = open(blast_file,"r")
    Scer_vs_Mega_records= NCBIXML.parse(Scer_vs_Mega)
    
    item = next(Scer_vs_Mega_records)
    gene_name = item.query
    Nter = item.query_length
    Cter = 0
    positions_dec = {'Nter':'','Cter':''}
    
    while item != 'FINISH':
        al_count = 0
        dico[item.query] = {}     
        #print(item.query)
        s = '-'*item.query_length
        for alignment in item.alignments:
            hsp_count = 0
            al_count += 1 
            
            for hsp in alignment.hsps:
                hsp_count += 1
                #print(alignment.hit_def,hsp.expect)
                #print(alignment.hit_def,hsp.expect,hsp.identities,hsp.align_length)
                if hsp.expect <= 1e-02:
                    dico[item.query][alignment.hit_def] = {}
                    dico[item.query][alignment.hit_def]['coverage'] = round(hsp.align_length/item.query_length ,1)
                    dico[item.query][alignment.hit_def]['align'] = '-'*(hsp.query_start-1) + hsp.sbjct + '-'*(item.query_length - hsp.query_end)
                    dico[item.query][alignment.hit_def]['frag_cov'] = round(hsp.align_length / alignment.length,1)
                    dico[item.query][alignment.hit_def]['position'] = 'Center'
                    dico[item.query][alignment.hit_def]['align_HCA'] = '-'*(hsp.query_start-1) + frag_HCA[alignment.hit_def]['sequence'][hsp.sbjct_start-1:hsp.sbjct_end] + '-'*(item.query_length - hsp.query_end)
                    s = s[:hsp.query_start-1] + hsp.query + s[hsp.query_end:] 
                    
                    Nter,Cter,positions_dec = decide_NandC_ter(Nter,Cter,positions_dec,hsp,alignment.hit_def)
                        
                    #hsp.
                else:
                    my_bad_hit = s[hsp.query_start-1:hsp.query_end]
                    if ((hsp.identities/hsp.align_length) + my_bad_hit.count('-') / len(my_bad_hit)) / 2 > 0.5 and len(my_bad_hit)>5:
                        #print(alignment.hit_def,hsp.expect,hsp.identities,hsp.align_length,my_bad_hit.count('-') / len(my_bad_hit),hsp.query)
                        dico[item.query][alignment.hit_def] = {}
                        dico[item.query][alignment.hit_def]['coverage'] = round(hsp.align_length/item.query_length ,1)
                        dico[item.query][alignment.hit_def]['align'] = '-'*(hsp.query_start-1) + hsp.sbjct + '-'*(item.query_length - hsp.query_end)
                        dico[item.query][alignment.hit_def]['frag_cov'] = round(hsp.align_length / alignment.length,1)
                        Nter,Cter,positions_dec = decide_NandC_ter(Nter,Cter,positions_dec,hsp,alignment.hit_def)
                        dico[item.query][alignment.hit_def]['position'] = 'Center'
                        dico[item.query][alignment.hit_def]['align_HCA'] = '-'*(hsp.query_start-1) + frag_HCA[alignment.hit_def]['sequence'][hsp.sbjct_start-1:hsp.sbjct_end] + '-'*(item.query_length - hsp.query_end)

                #dico[item.query][alignment.hit_def]
        
                if hsp_count == 1:
                    break

            
        try:
            item=next(Scer_vs_Mega_records)
        except:
            item='FINISH'
    
    if positions_dec['Nter'] == positions_dec['Cter']:
        dico[gene_name][positions_dec['Nter']]['position'] = 'Long'
        
    else:
        dico[gene_name][positions_dec['Nter']]['position'] = 'Nter'
        dico[gene_name][positions_dec['Cter']]['position'] = 'Cter'
    return(dico)


def read_HCA(hca_file):
    file=open(hca_file,'r').readlines()
    #print(hca_file)
    dico_HCA = {}
    for x,line in enumerate(file):
        found_dom=''
        if line.startswith('>'):
            domain_count = 0
            protein = []
            IGORF = line.split()[0].split('>')[1]
            dico_HCA[IGORF] = {}
            prot_size = line.split()[1]
            prot_score = line.split()[-1]
            # Generate the list of HCA clusters representation
            protein = ['.' for i in range(int(prot_size))]
            #dico_HCA[IGORF]['HCA_score'] = {}
            dico_HCA[IGORF]['HCA_score'] = prot_score
            dico_HCA[IGORF]['sequence'] = ''
            dico_HCA[IGORF]['Domains'] = {}
            
            # When HCA score is -10 and has not domain or cluster:
            if prot_score == '-10.000':
                if file[x+1].startswith('>'):
                    dico_HCA[IGORF]['sequence'] = protein
                    continue
        # We keep the domains  
        if line.startswith('domain'):
            domain_count = domain_count + 1
            domain_name = 'Domain'+str(domain_count)
            domain_start = line.split()[1]
            domain_stop  = line.split()[2]
            domain_pval  = line.split()[3]
            domain_score = line.split()[4]
            ####print ('domain',domain_start,domain_stop,domain_pval)
    
            dico_HCA[IGORF]['Domains'][domain_name] = {}
            dico_HCA[IGORF]['Domains'][domain_name]['start'] = domain_start
            dico_HCA[IGORF]['Domains'][domain_name]['stop'] = domain_stop
            dico_HCA[IGORF]['Domains'][domain_name]['pval'] = domain_pval
            dico_HCA[IGORF]['Domains'][domain_name]['score'] = domain_score
            dico_HCA[IGORF]['Domains'][domain_name]['clusters'] = []
        
        # We assign the clusters in the domains:
        if line.startswith('cluster'):
            cluster_start = line.split()[1]
            cluster_stop  = line.split()[2]
            cluster_motif = line.split()[3]
            # Append the clusters into the protein sequence:
            position = int(cluster_start)-1
            for cl in cluster_motif:
                protein[position] = cl
                position = position + 1
            
            # Assign clusters on the sequence
            for dom in dico_HCA[IGORF]['Domains'].keys():
                if int(cluster_start) >= int(dico_HCA[IGORF]['Domains'][dom]['start']) and \
                   int(cluster_stop)  <= int(dico_HCA[IGORF]['Domains'][dom]['stop']):
                   found_dom = dom
                   break
            if found_dom != '':
                dico_HCA[IGORF]['Domains'][found_dom]['clusters'].append( (cluster_motif,cluster_start,cluster_stop) )
            
            # Finishes the reading of the IGORF and assigns the prot seq 
            if x == len(file)-1:
                dico_HCA[IGORF]['sequence'] = ''.join(protein)
            elif file[x+1].startswith('>'):
                dico_HCA[IGORF]['sequence'] = ''.join(protein)
    
    return(dico_HCA)


def refine_lalign_alignment(st):
    s_init = st
    s      = s_init
    s = s.replace(":.:",":::")
    s = s.replace(": :",":::")
    s = s.replace(" .:"," ::")
    s = s.replace(":. ",":: ")
    s_tmp = ""
    while(s_tmp != s):
        s = s.replace(":.:",":::")
        s = s.replace(": :",":::")
        s = s.replace(" .:"," ::")
        s = s.replace(":. ",":: ")
        s_tmp = s
    return(s)


def parse_Lalign(lalign,name):
    with open(lalign,'r') as f:
        sbjct = 'HELLO'
        label = 'OFF'
        dico = {}
        for n,line in enumerate(f):
            if line.startswith('>>>'+name):
                query = line.split('>>>')[1].split(',')[0]
                query_size = int(line.split()[1].strip())
                continue
            
            if line.startswith('>>') and line.strip() not in ['>>><<<','>>>///']:
                label = 'ON'
                sbjct = line.split()[0].split('>>')[1]
                sbjct_size = int(line.split()[1].split('(')[1])
                dico[sbjct] = {'SEQ':'','Eval':'','Q_start':'','Q_end':'','S_start':'','S_end':'','S_size':sbjct_size,'positions':''}
                continue
            
            if line.startswith(' Waterman-Eggert score:') and label=='ON':
                evalue = float(line.split('<')[1].strip())
                dico[sbjct]['Eval']=evalue
                continue
            
            if '% identity' in line and label == 'ON':
                identity    = float(line.split('% identity')[0])
                overlap     = int(line.split("aa overlap (")[0].split()[-1].strip())
                query_start = int(line.split('overlap (')[1].split(':')[0].split('-')[0])
                query_end   = int(line.split('overlap (')[1].split(':')[0].split('-')[1])
                sbjct_start = int(line.split('overlap (')[1].split(':')[1].split('-')[0])
                sbjct_end   = int(line.strip().split('overlap (')[1].split(':')[1].split('-')[1].split(')')[0])
                dico[sbjct]['identity']= identity
                dico[sbjct]["overlap"]   = overlap 
                dico[sbjct]['Q_start'] = query_start
                dico[sbjct]['Q_end']   = query_end
                dico[sbjct]['S_start'] = sbjct_start
                dico[sbjct]['S_end']   = sbjct_end
    
            if line.startswith(sbjct[0:6]) and label == 'ON':
                dico[sbjct]['SEQ'] = dico[sbjct]['SEQ'] + line.split()[1]
            if line.startswith(name[0:6]) and label == 'ON':
                #dico[query]['SEQ'] = dico[sbjct]['SEQ'] + line.split()[1]
                new_line = str(f.readline().strip())
                positions = refine_lalign_alignment(st=new_line)
                dico[sbjct]['positions'] = dico[sbjct]['positions'] + positions
                
            if line.strip() == '>--':
                label='OFF'
    
    return (dico)



def decide_NandC_ter_Lalign(Nter_tmp,Cter_tmp,Dec_dico,frag_tmp,name):
    if frag_tmp['Q_start'] < Nter_tmp:
        Nter_tmp = frag_tmp['Q_start']
        Dec_dico['Nter'] = name
    if frag_tmp['Q_end'] > Cter_tmp:
        Cter_tmp = frag_tmp['Q_end']
        Dec_dico['Cter'] = name
    return(Nter_tmp,Cter_tmp,Dec_dico)

def ranges(nums):
    '''
    Finds the consecutive indexes of numbers
    in order to find the overlapping region between 2
    series of numbers
    '''
    #nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))

def decrease_aligned_position(dico):
    y = [i+1 for i,pos in enumerate(dico["positions"]) if pos == ":"]
    y_r = ranges(nums=y)
    y_rs = [(len(list(range(i[0],i[1])))) for i in y_r]
    max_y = y_r[y_rs.index(max(y_rs))]
    dico['Q_start'] = dico['Q_start'] + max_y[0]
    dico['Q_end']   = dico['Q_end'] - (len(dico["SEQ"])-max_y[1])
    dico['S_start'] = dico['S_start'] + max_y[0]
    dico['S_end']   = dico['S_end'] - (len(dico["SEQ"])-max_y[1])
    new_SEQ = dico["SEQ"][ max_y[0]: max_y[1]]
    dico["SEQ"] = new_SEQ
    return(dico)


def decide_fragments_aligned(LALIGN,gsize,frag_HCA):
    dico = {}    
    Nter = gsize
    Cter = 0
    positions_dec = {'Nter':'','Cter':''}
    s = '-'*gsize
    keys_to_delete = []
    for frag in LALIGN:
        
        LALIGN[frag] = decrease_aligned_position(dico=LALIGN[frag])
        
        if LALIGN[frag]['Eval'] <= 0.01:
            try:
                #LALIGN[frag]['coverage'] = round(len(LALIGN[frag]['SEQ']) / gsize,2)
                LALIGN[frag]['coverage'] = round(LALIGN[frag]['overlap'] / gsize,2)
                LALIGN[frag]['frag_coverage'] = round(len(LALIGN[frag]['SEQ'].replace('-',''))/LALIGN[frag]['S_size'],2)
                #print(LALIGN[frag]['Q_start'] , LALIGN[frag]['SEQ'] , LALIGN[frag]['Q_end'])
                LALIGN[frag]['align'] =  '-'*(LALIGN[frag]['Q_start']-1) + LALIGN[frag]['SEQ'] + '-'*(gsize - LALIGN[frag]['Q_end'])
                LALIGN[frag]['HCA_align'] = '-'*(LALIGN[frag]['Q_start']-1) + \
                frag_HCA[frag]['sequence'][LALIGN[frag]['S_start']-1:LALIGN[frag]['S_end']] \
                + '-'*(gsize - LALIGN[frag]['Q_end'])
                LALIGN[frag]['position'] = 'Center'
                s = s[:LALIGN[frag]["Q_start"]-1] + LALIGN[frag]["SEQ"] + s[LALIGN[frag]["Q_end"]:]
                
                dico[frag] = LALIGN[frag]
                Nter,Cter,positions_dec = decide_NandC_ter_Lalign(Nter,Cter,positions_dec,LALIGN[frag],frag)
            except:
                keys_to_delete.append(frag)
                continue
        else:
            try:
                my_bad_hit = s[LALIGN[frag]["Q_start"]-1:LALIGN[frag]["Q_end"]]
                if ( LALIGN[frag]["identity"] + (my_bad_hit.count('-') / len(my_bad_hit))*100) / 2 > 50 and len(my_bad_hit)>5:
                    LALIGN[frag]['coverage'] = round(len(LALIGN[frag]['SEQ']) / gsize,2)
                    LALIGN[frag]['frag_coverage'] = round(len(LALIGN[frag]['SEQ'].replace('-',''))/LALIGN[frag]['S_size'],2)
                    LALIGN[frag]['align'] =  '-'*(LALIGN[frag]['Q_start']-1) + LALIGN[frag]['SEQ'] + '-'*(gsize - LALIGN[frag]['Q_end'])
                    LALIGN[frag]['HCA_align'] = '-'*(LALIGN[frag]['Q_start']-1) + \
                    frag_HCA[frag]['sequence'][LALIGN[frag]['S_start']-1:LALIGN[frag]['S_end']] \
                    + '-'*(gsize - LALIGN[frag]['Q_end'])
                    LALIGN[frag]['position'] = 'Center'
                    s = s[:LALIGN[frag]["Q_start"]-1] + LALIGN[frag]["SEQ"] + s[LALIGN[frag]["Q_end"]:]
                    
                    dico[frag] = LALIGN[frag]
                    Nter,Cter,positions_dec = decide_NandC_ter_Lalign(Nter,Cter,positions_dec,LALIGN[frag],frag)
            except:
                keys_to_delete.append(frag)
                continue

    if positions_dec['Nter'] == positions_dec['Cter']:
        dico[positions_dec['Nter']]['position'] = 'Long'
        
    else:
        dico[positions_dec['Nter']]['position'] = 'Nter'
        dico[positions_dec['Cter']]['position'] = 'Cter'
        
    for frag in keys_to_delete:
        del LALIGN[frag]
    return (dico)


