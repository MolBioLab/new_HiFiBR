#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 22:05:58 2021

@author: mobilab5
"""

import os
import pandas as pd
import csv
import re
from Bio import pairwise2
import sys
import subprocess


def cigar_substring_list(cigar_string):
    MDI_list = ["M", "D", "I"]
    cigar_substring_list = []
    cigar_len = len(cigar_string)
    i = 0
    while i < cigar_len:
        sub_cigar = ""
        while cigar_string[i] not in MDI_list:
            sub_cigar = sub_cigar + cigar_string[i]
            i += 1
        sub_cigar = sub_cigar + cigar_string[i]
        cigar_substring_list.append(sub_cigar)
        i += 1
    return cigar_substring_list

def match_mis(read,ref):
    mismatch_L = []
    for i in range(len(read)):
        if read[i] == ref[i]:
            mismatch_L.append('0')
        else:
            mismatch_L.append('1')
    return mismatch_L

def check_unsure(mismatch_L, read, side):
    num_mis = mismatch_L.count("1")
    if num_mis >= 1:
        unsure.write("\t".join(line) + "\n")
        
        with open(sample+"_Unsure_log.txt", "a+") as unsure_log:
            if side == "l": unsure_log.write(line[0] +"\t Side: left\n")
            else: unsure_log.write(line[0] +"\t Side: right\n")
        raise Exception
        
    else:
        num_appear = ref.count(read)
        if num_appear >= 2:
            unsure.write("\t".join(line)  + "\n")
            
            with open(sample+"_Unsure_log.txt", "a+") as unsure_log:
                if side == "l": unsure_log.write(line[0] +"\t Side: left\n")
                else: unsure_log.write(line[0] +"\t Side: right\n")
                
            raise Exception
            
def writing_temp(new_line, mismatch_l, mismatch_r):    
    with open(sample+"_temp_analyse.fq", "a+") as no_auto:
        no_auto.write("\n".join(new_line) + "\n")
        
    with open(sample+"_AutoOrNot_log.txt", "a+") as AutoOrNot:
        if mismatch_l.count("1") > mismatch_r.count("1"):
            AutoOrNot.write(line[0] +"\t Side: left\n")
        else:
            AutoOrNot.write(line[0] +"\t Side: right\n")
            
def check_AutoOrNot(mismatch_l, mismatch_r, read_l, read_r, qual_l, qual_r):     
    ref_start = start_pos + M_l
    ref_end = end_pos - M_r
    pad_seq = ref[ref_start:ref_end]
    new_read = read_l + pad_seq + read_r
    new_qual = qual_l + "G"*len(pad_seq) + qual_r
    new_line = ["@"+line[0], new_read, "+", new_qual]   
    writing_temp(new_line, mismatch_l, mismatch_r)


# RUN 
sample = sys.argv[1]
wd_path = sys.argv[2] + "/" + sample + "/2_mapping"
os.chdir(wd_path) 
sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)
from ExtendedHiFiBR import SeqInfo as si

ref = si.get_refSeq(sys.argv[3])      
reader = csv.reader(open(sample+'_NoHairpin.sam',"r"), dialect="excel-tab")
cmd = "samtools view -H " + sample + '_len.sam'
proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, )
header=proc.communicate()[0].decode("utf-8")
   
with open(sample+"_AutoAnalyse.sam", "a+") as auto:
    auto.write(header)

unsure = open(sample+"_Unsure.sam", "a+") 
unsure.write(header)

#
count_realign = 0  
count_auto = 0

for line in reader:
    if line[0][0] == "@":
        pass    
    else:
        try:
            # A. Get info
            start_pos = int(line[3]) - 1
            end_pos = int(line[-1].split(":")[2])
            cigar_list = cigar_substring_list(line[5]) 
            
            M_l = 11 if len(cigar_list) == 1 else int(cigar_list[0][:-1])
            M_r = 11 if len(cigar_list) == 1 else int(cigar_list[-1][:-1])
            read_seq = line[9]
            qual = line[10]
            
            read_l = read_seq[:M_l]
            read_r = read_seq[-M_r:]  
            qual_l = qual[:M_l]
            qual_r = qual[-M_r:]
            
            ref_l = ref[start_pos:M_l+start_pos]
            ref_r = ref[(end_pos - M_r):end_pos]           

            # B. Mismatch            
            mismatch_l = match_mis(read_l,ref_l)    # 1. List of match (0) / mismatch (1) 
            mismatch_r = match_mis(read_r,ref_r)    
            
            # C. Check unsure
            if len(cigar_list) == 1: # No need to check unsure
                new_line = ["@"+line[0], line[9], "+", line[10]]
                count_realign += 1
                writing_temp(new_line, mismatch_l, mismatch_r)
                raise Exception
            else:    
                if (M_l < 12):                          # 2. Check if M_l < 12 or M_r < 12 (unsure mapped / split)
                    res = check_unsure(mismatch_l, read_l, "l")
                if (M_r < 12):
                    res = check_unsure(mismatch_r, read_r, "r")
                    
            # D. Check cantAuto side  
                if (mismatch_l.count("1") == 0) and (mismatch_r.count("1") == 0):
                    with open(sample+"_AutoAnalyse.sam", "a+") as auto:
                        auto.write("\t".join(line)+"\n")
                    count_auto +=1
                else:
                    check_AutoOrNot(mismatch_l, mismatch_r, read_l, read_r, qual_l, qual_r)
                    count_realign +=1
                
        except Exception:
            continue
print("--- REPORT --- \nCount Auto Analyze Junctions",count_auto)
print("Count Temporary Recontruction Junctions ",count_realign)

