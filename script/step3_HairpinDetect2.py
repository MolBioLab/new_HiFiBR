#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 17:50:05 2021

@author: mobilab5
"""
import os
import csv
import pandas as pd
import subprocess
from Bio.Seq import Seq
import sys


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

def mismatch_num(M_l, M_r):
    read_l = read[:M_l]
    ref_l = ref[start_pos:M_l+start_pos]
    mismatch_L_l = match_mis(read_l, ref_l)
    mismatch_num_l = mismatch_L_l.count("1")
      
    read_r = read[-M_r:]
    ref_r = ref[(end_pos - M_r):end_pos]
    mismatch_L_r = match_mis(read_r, ref_r)
    mismatch_num_r = mismatch_L_r.count("1")    
    return mismatch_num_l, mismatch_num_r
            
def check_hairpin(read_seq, read2_start, side):
    if side == "l":
        i = 11 # read_start_l = nu 12th
        a = read2_start -1
        while True:
            try:
                if Seq(read_seq[i]).reverse_complement() != read_seq[a]:
                    break
                i += 1
                a -= 1
            except IndexError:
                pass
    if side == "r":
        i = -12
        a = read2_start + 11
        while True:
            try:
                if Seq(read_seq[i]).reverse_complement() == read_seq[a]:
                    break
                i -= 1
                a += 1
            except IndexError:
               pass
    return i,a
            
# RUN
sample = sys.argv[1]
wd_path = sys.argv[2] + "/" + sample + "/2_mapping"
os.chdir(wd_path) 
sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)
from ExtendedHiFiBR import SeqInfo as si

ref = si.get_refSeq(sys.argv[3]) 
reader = csv.reader(open(sample+'_len.sam',"r"), dialect="excel-tab")  
cmd = "samtools view -H " + sample + '_len.sam' 
proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, )
header=proc.communicate()[0].decode("utf-8")

hairpin = open(sample+"_Hairpin.sam", "a+")   
no_hairpin = open(sample+"_NoHairpin.sam", "a+")  

hairpin.write(header)
no_hairpin.write(header)

hairpin.close()
no_hairpin.close()

for line in reader:
    if line[0][0] == "@":
        pass    
    else:
        
        # A. Get info   
        start_pos = int(line[3]) - 1           
        end_pos = int(line[-1].split(":")[2]) 
        cigar = line[5]
        cigar_L = cigar_substring_list(cigar)
        read = line[9]
        read_seq = Seq(read)    
        
        read_seed_l = read_seq[:11] # hairpin often happens from the end
        read_seed_r = read_seq[-11:]  
        
        # B. Check if any mismatch both side
        M_l = int(cigar_L[0][:-1])  
        M_r = int(cigar_L[-1][:-1])
        
        M_l = 11 if M_l >= 11 else M_l
        M_r = 11 if M_r >= 11 else M_r
        mismatch_num_l, mismatch_num_r = mismatch_num(M_l, M_r)

        # C. Execute left side          
        read2_start_l = read_seq[11:].find(str(read_seed_l.reverse_complement()))
        if (read2_start_l == -1) or (mismatch_num_l < 3):
            pass
            
        else:        
            read2_start_l = read2_start_l + 11
            read_end_l, read2_end_l = check_hairpin(read_seq, read2_start_l, "l")
            read_hairpin = str(read_seq[:read_end_l])
            read2_hairpin = str(read_seq[read2_end_l + 1: read2_start_l + 11] )
            
            with open(sample+"_Hairpin_log.txt", "a+") as log: 
                log.write(line[0] + "\t Side: left" +
                          "\n Region1: " + read_hairpin + "\n Region2: " + read2_hairpin + "\n")
                # log.write("Reg1_start: " + str(read_end_l) + "\t Reg1_end: " + str(start_pos) + "\t")
                # log.write("Reg2_start: " + str(read2_end_l + 1) + "\t Reg2_end: " + str(read2_start_l + 11) + "\n")
            with open(sample+"_Hairpin.sam", "a+") as hairpin: 
                hairpin.write("\t".join(line)  + "\n")
            continue # Dont need to execute right side
            
        # C. Execute right side
        read2_start_r = read_seq[:-11].find(str(read_seed_r.reverse_complement()))  
        if (read2_start_r == -1) or (mismatch_num_r < 3):
            pass    
        else:      
            read_end_r, read2_end_r = check_hairpin(read_seq, read2_start_r,"r")
            read_hairpin = str(read_seq[read_end_r+1:])
            read2_hairpin = str(read_seq[read2_start_r:read2_end_r])        

            with open(sample+"_Hairpin_log.txt", "a+") as log:
                log.write(line[0] + "\t Side: right" +
                          "\n Region1: " + read_hairpin + "\n Region2: " + read2_hairpin + "\n")
                # log.write("Reg1_start: " + str(read_end_r + 1) + "\t Reg1_end: " + str(end_pos) + "\t")
                # log.write("Reg2_start: " + str(read2_start_r) + "\t Reg2_end: " + str(read2_end_r) + "\n")
            with open(sample+"_Hairpin.sam", "a+") as hairpin: 
                hairpin.write("\t".join(line)  + "\n")
            continue
                
        # D. If neither hairpin left of right side
        
        with open(sample+"_NoHairpin.sam", "a+") as no_hairpin: # fix
            no_hairpin.write("\t".join(line)  + "\n")
                    

