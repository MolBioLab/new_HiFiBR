sample=$1
ANALYSE_PATH=$2
cutPos1=$3 
cutPos2=$4
ref_path=$5
SOFTWARE_PATH=$6

alias python=python3.6
# A. Fix multiple alignment info from each read
# Input: {sample}_aligned1.sam & {sample}_aligned2.sam
# Output: {sample}_fixed_aligned1.sam & {sample}_fixed_aligned2.sam
echo "** STEP 3.1: FIX SECONDARY ALIGNMENTS **"
python step3_FixAligned.py ${sample} ${ANALYSE_PATH} $SOFTWARE_PATH

# B. GET JUNCTIONS which align over 2 RE positions on reference 
# Input: {sample}_fixed_aligned1.sam & {sample}_fixed_aligned2.sam
# Output: {sample}_FilteredAligned1.sam & {sample}_FilteredAligned2.sam
echo "** STEP 3.2: GET Junctions Aligned over 2 RE positions **"
./step3_FilterMapped.sh ${sample} $cutPos1 $cutPos2


# C. ALTER & ADD alignment info
#	1. Fix NGS Error
echo "** STEP 3.3: ALTER & ADD alignment info **"
samtools merge ${sample}_FilteredAligned.sam ${sample}_FilteredAligned1.sam ${sample}_FilteredAligned2.sam 
python step3_FixedMappedSam2.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH
# 	2. Add End Position to sam file
python step3_Len.py ${sample}_NoNgsError.sam > ${sample}_temp.sam
samtools view -h ${sample}_temp.sam > ${sample}_len.sam
rm ${sample}_temp.sam


# D. HAIRPIN
echo "** STEP 3.4: Detect & remove HAIRPIN **"
python step3_HairpinDetect2.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH


# E. JUNCTION CLASSIFICATION: Auto / CantAuto / unsure
echo "** STEP 3.5: Detect & Remove UNSURE **"
python step3_Extract_CantAutoAnalyse5.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH

	# E.1. RE-MAP temporary reconstruction read with mismatches
echo "*** RE-MAP Temporary Reconstruction Read ***"
bowtie2-build -f ${ref_path} ref
bowtie2 -x ref --local -U ${sample}_temp_analyse.fq -S ${sample}_temp_analyse.sam
cat ${sample}_temp_analyse.fq | grep "^@M" | sort | uniq | sed 's/^.//'> ${sample}_temp_name.txt

	# E.2. If "H"/"S" in cigar -> goes into "cant auto" bin
echo "*** CLASSIFY Junction to CantAuto bin ***"
samtools view ${sample}_temp_analyse.sam | awk -F"\t" '$6~/S/ || $6~/H/ {print $1}' | sort | uniq > ${sample}_cant_auto_name.txt
samtools view -f 4 ${sample}_temp_analyse.sam | awk -F"\t" '{print $1}' | sort | uniq >> ${sample}_cant_auto_name.txt
while read -r line;do echo "${line}\t" >> ${sample}_cant_auto_name2.txt;done < ${sample}_cant_auto_name.txt
samtools view -H ${sample}_NoHairpin.sam > ${sample}_CantAutoAnalyse.sam
cat ${sample}_NoHairpin.sam | rg -j 12 -f ${sample}_cant_auto_name2.txt >> ${sample}_CantAutoAnalyse.sam

	# E.3. If not -> goes in to "auto" bins
echo "*** CLASSIFY Junction to Auto bin ***"
sort ${sample}_cant_auto_name.txt > ${sample}_cant_auto_name3.txt
comm -23 ${sample}_temp_name.txt ${sample}_cant_auto_name3.txt > ${sample}_auto_name.txt
while read -r line;do echo "${line}\t" >> ${sample}_auto_name2.txt;done < ${sample}_auto_name.txt
cat ${sample}_NoHairpin.sam | rg -j 12 -f ${sample}_auto_name2.txt >> ${sample}_AutoAnalyse.sam

${sample}_cant_auto_name.txt ${sample}_auto_name.txt ${sample}_temp_name.txt ${sample}_cant_auto_name3.txt ref.*

# F. Remove tooShort (Hinh 2.4E)
#	python step3_RemoveTooShort5.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH
#	while read -r line;do echo "${line}\t" >> pass_auto2.txt;done < pass_auto.txt
#	while read -r line;do echo "${line}\t" >> dropby_tooShort2.txt;done < dropby_tooShort.txt

#	samtools view -H ${sample}_AutoAnalyse.sam > ${sample}_AutoAnalyse2.sam
#	cat ${sample}_AutoAnalyse.sam | rg -j 12 -f pass_auto2.txt >> ${sample}_auto_analyse2.sam
#	samtools view -H ${sample}_auto_analyse.sam > ${sample}_dropby_tooShort.sam
#	cat ${sample}_auto_analyse.sam | rg -j 12 -f dropby_tooShort2.txt >> ${sample}_dropby_tooShort.sam

#	rm pass_auto2.txt pass_auto.txt dropby_tooShort2.txt dropby_tooShort.txt

	

