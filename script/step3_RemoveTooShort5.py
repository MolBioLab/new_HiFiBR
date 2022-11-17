#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 00:11:48 2021

@author: mobilab5
"""

import sys
import os
import csv
from diff_match_patch import diff_match_patch
import concurrent.futures as ccf
import collections
from itertools import repeat, filterfalse
import itertools
from tqdm import tqdm

def compare2dict(seq, in_D):
    unmatch_count = 0
    for key in in_D.keys():
        dmp = diff_match_patch()
        diff_L = dmp.diff_main(key, seq)
        stat_L = [k[0] for k in diff_L] # del = -1, ins = 1, equal = 0 of seq(considering) compare to seq(key)
        
    # A. Only change key name & add value   
        # 1. seq(considering) totally inside/exact the same seq(key)
        if stat_L == [-1, 0, -1] or stat_L == [-1, 0] or stat_L == [0, -1] or stat_L == [0]:
            return ["done1", key, None]
        # 2. seq(considering) with longer left-side
        elif stat_L == [1, 0, -1] or stat_L == [1, 0]:
            return ["done2", key, diff_L]                              
        # 3. seq(considering) with longer right-side
        elif stat_L == [-1, 0, 1] or stat_L == [0, 1]:
            return ["done3", key, diff_L]
        # 4. seq(considering) with longer both-side
        elif stat_L == [1, 0, 1]:
            return ["done4", key, None] 
    # B. Add new key to dict    
        else:
            unmatch_count += 1
        
        if unmatch_count == len(in_D.keys()):
            return ["not yet", None, None]   

def alt_dict(seq_L, name_L, res_L, in_D, cate):
# if cate = "single" -> normal. if cate = "multi" -> add extend group
# add cate argu to alt_dict used    
    def work(seq, name, status, key, concat_seq, in_D):
        if status == "done1":
            in_D[key].extend(name)
        elif status == "done2" or status == "done3":
            in_D = {concat_seq if k == key else k:v for k,v in in_D.items()} # change seq(key)
            in_D[concat_seq].extend(name) 
        elif status == "done4":
            in_D = {seq if k == key else k:v for k,v in in_D.items()}
            in_D[seq].extend(name)   
    
    def work2(word, diff, side): # create concat_seq with longest left/right-side
        pos = [j for j,x in enumerate(status) if x == word] 
        if len(pos) > 0:
            temp = 0 if side == "l" else 1
            diffn_temp = max([diff[k][temp][1] for k in pos], key = len)  
            diffn = [diff[k] for k in pos if diff[k][temp][1] == diffn_temp][0]
            concat_seq = "".join( [k[1] for k in diffn] )
        else:
            concat_seq = ""
        return concat_seq
    # Check uniq / dup
    key_L = [i[1] for i in res_L]
    uniq = [key_L.index(i) for i, c in collections.Counter(key_L).items() if c == 1]
    dup = [[j for j,x in enumerate(key_L) if x == item]
           for item, c in collections.Counter(key_L).items() if c > 1]
    
    # Execute uniq
    for i in uniq:
        seq, name, status, key, diff_L = seq_L[i], name_L[i], res_L[i][0], res_L[i][1], res_L[i][2]
        concat_seq = "".join( [k[1] for k in diff_L] ) if diff_L != None else ""
        
        if cate == "single":
            work(seq, [name], status, key, concat_seq, in_D)    
        elif cate == "multi":
            work(seq, name, status, key, concat_seq, in_D)
                        
    # Execute dup
    for i_L in dup:
        seq = [seq_L[k] for k in i_L]  
        name = [name_L[k] for k in i_L] 
        if cate == "multi": name = list(itertools.chain.from_iterable(name))         
        status = [res_L[k][0] for k in i_L]
        key = [res_L[k][1] for k in i_L]
        diff = [res_L[k][2] for k in i_L]
    
        if "done4" in status:
            pos = [j for j,x in enumerate(status) if x == "done4"]
            seq2_L = [seq[k] for k in pos]
            seq_max = max(seq2_L, key=len)
            work(seq_max, name, "done4", key[0], None, in_D)
            
        elif "done2" in status or "done3" in status:
            concat2 = work2("done2", diff, "l")
            concat3 = work2("done3", diff, "r")
            
            dmp = diff_match_patch()    # merge longest left- & right-side seq
            diff2 = dmp.diff_main(concat2, concat3)
            concat_seq = "".join( [k[1] for k in diff2] )
            work(None, name, "done2", key[0], concat_seq, in_D)
        
        elif "done1" in status:
            work(None, name, "done1", key[0], None, in_D)
    return in_D        

def compare_each_other(sub_seq, sub_name, in_D, cate):
    while True:
        if len(sub_name) == 0:
            break
        if len(sub_name) == 1:
            in_D.update({sub_seq[0] : [sub_name[0]] })
            break
                   
        in2_D = {sub_seq[0] : [sub_name[0]] }        
        with ccf.ProcessPoolExecutor() as executor:
            results2 = list(executor.map(compare2dict, sub_seq[1:], repeat(in2_D)))
        done_i = [j for j,x in enumerate(results2) if "done" in x[0]]
            
        if len(done_i) > 0: # if there's any seq similar to [0] item               
            sub_seq_done2 = [sub_seq[1:][k] for k in done_i]
            sub_name_done2 = [sub_name[1:][k] for k in done_i]                
            results_done2 = [results2[k] for k in done_i]
            in2_D = alt_dict(sub_seq_done2, sub_name_done2, results_done2, in2_D, cate)
            
            for k in sub_seq_done2: sub_seq.remove(k) 
            for k in sub_name_done2: sub_name.remove(k)   
            
        in_D.update(in2_D) 
        del sub_seq[0], sub_name[0]  
    return in_D

def grouping(f_name):
    # A. Get name, cigar, seq list & group info of 3 files    
    reader = csv.reader(open(f_name,"r"), dialect="excel-tab")
    name1_L, seq_L = [], []
    for line in reader:
        if line[0][0] == "@":
            pass
        else:
            name1_L.append(line[0])
            seq_L.append(line[9])
            
    name_L = [x for _,x in sorted(zip(seq_L, name1_L), key=lambda pair: pair[0])]        
    seq_L.sort()  
    
    # B. Grouping remain seqs  
    group_D = {seq_L[0]: [name_L[0]]}
    for i in tqdm(range(1,len(name_L), 10)): 
        sub_name, sub_seq = name_L[i:i+10], seq_L[i:i+10] # Get sub_seq (10 seqs)
        
        # B1. Compare sub_seq with key_seq (in group_D)
        with ccf.ProcessPoolExecutor() as executor:
            results = list(executor.map(compare2dict, sub_seq, repeat(group_D)))
            
        sub_name_notYet, sub_seq_notYet, sub_name_done, sub_seq_done, results_done = [], [], [], [], []
        for j in range(len(results)):
            if results[j][0] == "not yet": # Continue with B3
                sub_name_notYet.append(sub_name[j]) 
                sub_seq_notYet.append(sub_seq[j])
            elif "done" in results[j][0]:  # Compare results & write into group_D
                sub_name_done.append(sub_name[j]) 
                sub_seq_done.append(sub_seq[j]) 
                results_done.append(results[j]) 
        
        if len(sub_name_done) > 0:
            group_D = alt_dict(sub_seq_done, sub_name_done, results_done, group_D, "single")
        
        # B2. Compare n seqs (which are not contributed into groups) with each others    
        if len(sub_name_notYet) == 1:
            group_D.update({sub_seq_notYet[0] : [sub_name_notYet[0]] })
        elif len(sub_name_notYet) > 1:             
            group_D = compare_each_other(sub_seq_notYet, sub_name_notYet, group_D, "single")
    
    # Print out info for checking
    print("Groups: ", len(group_D.keys()))
    count = 0
    for k,v in group_D.items():
        count += len(v)
    print("Seq num: ", count)    
    
    # Only get group with more than n seqs (n=2 -> put into argu when finished script)
    filtered_group = {k:v for k,v in group_D.items() if len(v) >= 2}
    
    # Sort dict by length of items in groups   
    sorted_group = {k:v for k, v in sorted(filtered_group.items(), 
                                           key=lambda item: len(item[1]), 
                                           reverse= True )}
    info_group = dict()
    count_line = 0
    for k, v in sorted_group.items():
        count_line +=1
        k2 = (k, "Group " + str(count_line))
        info_group[k2] = v
    return info_group
    
def writing_log(text, f_side, info_group, num = 1):  
    with open(f_side,"r") as f:
        lines = f.read().splitlines()
    info_lines = lines[::num]
    info_side = [k.split("\t") for k in info_lines]
    name_info = [k[0] for k in info_side]
    side_info = [k[1] for k in info_side]
    
    for k, v in info_group.items():
        with open(sample + "_Group.txt","a+") as log: # k[0]: seq, k[1]: group num
            # Get all name (in value) caused-drop side
            index_info = [name_info.index(k) for k in v]
            side_L = [side_info[k] for k in index_info]
            
            if side_L.count("left") > side_L.count("right"):
                log.write(k[1] + "_" + text + "\t Seq: " + k[0] + "\t Count: " + str(len(v)) + 
                          "\tSide: left" + "\n")
            else: 
                log.write(k[1] + "_" + text + "\t Seq: " + k[0] + "\t Count: " + str(len(v)) + 
                          "\tSide: right" + "\n")
            for i in range(len(v)):
                log.write(v[i] + "\t" + side_L[i] + "\n")                  
# RUN: 
sample = sys.argv[1]
wd_path = sys.argv[2] + "/" + sample + "/2_mapping"
os.chdir(wd_path) 
sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)
from ExtendedHiFiBR import SeqInfo as si 
ref = si.get_refSeq(sys.argv[3])
          
# A. Grouping unsure, cantAuto, auto
print("---REPORT---")
print("*** Hairpin ***")
hairpin_group = grouping(sample + "_Hairpin.sam")
writing_log("hairpin", sample + "_Hairpin_log.txt", hairpin_group, 3)

print("*** Unsure ***")
unsure_group = grouping(sample + "_Unsure.sam")
writing_log("unsure", sample + "_Unsure_log.txt", unsure_group)

print("*** CantAuto ***")
cantAuto_group = grouping(sample + "_CantAutoAnalyse.sam")
writing_log("cantAuto", sample + "_AutoOrNot_log.txt", cantAuto_group)



