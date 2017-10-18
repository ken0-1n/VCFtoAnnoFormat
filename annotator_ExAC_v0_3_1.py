import sys
import os
import re
import logging
import pysam

vcffile = open(sys.argv[1],'r')

for line in vcffile:

    if line.startswith("#"):continue
    line = line.rstrip()
    record = line.split('\t')

    chr_org = record[0]
    start_org = record[1]
    end_org = record[1]
    rsid_freq = record[2]
    ref_org = record[3]
    alt_org = record[4]
    qual_org = record[5]
    filter_org = record[6]
    info_org = record[7]

    fh ={'AC_Adj':"", 'AN_Adj':"", 'AC_POPMAX':"", 'AN_POPMAX':"", 'POPMAX':""}

    for info_item in info_org.split(';'):
        array = info_item.split('=')
        if len(array) > 1:
            key = array[0]
            value = array[1]
            if key != None and key in fh: 
                fh[key] = value

    idx = 0
    alt_list = alt_org.split(',')
    ac_adj_list = fh['AC_Adj'].split(',')
    an_adj_list = fh['AN_Adj'].split(',')
    ac_popmax_list = fh['AC_POPMAX'].split(',')
    an_popmax_list = fh['AN_POPMAX'].split(',')
    popmax_list = fh['POPMAX'].split(',')
    for alt in alt_list:
        start = int(start_org)
        ref = ref_org
        ref_len = len(ref)
        alt_len = len(alt)

        ac_adj = ac_adj_list[idx]
        an_adj = an_adj_list[idx]
        ac_popmax = ac_popmax_list[idx]
        an_popmax = an_popmax_list[idx]
        popmax = popmax_list[idx]
      
        # SNV
        if alt_len == 1 and ref_len == 1:
            end = int(end_org)

        # deletion
        elif alt_len < ref_len and ref.startswith(alt):
            ref = ref[alt_len:]
            start = start + ref_len - len(ref)
            end = start + len(ref) - 1
            alt = '-'

        # deletion
        elif alt_len > ref_len and alt.startswith(ref):
            alt = alt[ref_len:]
            start = start + alt_len - len(alt)
            end = start
            ref = '-'

        # block substitution 1
        else:
            continue

        freq_adj = '---'
        freq_popmax = '---'
        # try:
        if ac_adj != 'NA' and an_adj != 'NA' and int(an_adj) > 0:
            freq_adj = '{0:.6f}'.format(float(ac_adj) / float(an_adj))
        if ac_popmax != 'NA' and an_popmax != 'NA' and int(an_popmax) > 0:
            freq_popmax = '{0:.6f}'.format(float(ac_popmax) / float(an_popmax))
        # except ValueError:
        #     continue
        print "chr"+ chr_org +"\t"+ str(start-1) +"\t"+ str(end) +"\t"+ ref +"\t"+ alt +"\t"+ filter_org +';;'+ ac_adj +';;'+ an_adj +';;'+ str(freq_adj) +';;'+ ac_popmax +';;'+ an_popmax +';;'+ str(freq_popmax) +';;'+ popmax
                
####
vcffile.close()


