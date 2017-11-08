#!/bin/sh

export EXEC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# -----------------------------------------------------------------------------
#                                       Info
# -----------------------------------------------------------------------------

set -e

if [ "$#" -ne 3 ]; then
    echo "usage: CONFIG CHROM SV_START"
    exit 1
fi

CONFIG=$1
CHROM=$2
SV_START=$3

if [ "${CONFIG:0:1}" != "/" ] # no absolute path
    then
    CONFIG="$(pwd)/$CONFIG"
fi

if [ ! -f "$CONFIG" ]
    then
    echo "[ERROR - sweep_variant.sh] Could not find file $CONFIG."
    exit 1
fi

source "$CONFIG"

DATE=$(date +"%d-%B-%Y")
LOG="$(pwd)/sv.$CHROM.$SV_START.log" # logs all the programm output
echo $DATE > $LOG                    # initialize logfile

BAMFILE_ILLUMINA_HEADER="$POLISHING_DIR/BAMFILE_ILLUMINA.header"
samtools view -H "$BAMFILE_ILLUMINA" -o "$BAMFILE_ILLUMINA_HEADER"

BAMFILE_ONT_HEADER="$POLISHING_DIR/BAMFILE_ONT.header"
samtools view -H "$ONT_BAM_FILE" -o "$BAMFILE_ONT_HEADER"


echo "" | tee -a $LOG
echo -e "======================================================================" | tee -a $LOG
echo -e "                       POLISHING VARIANT $CHROM:$SV_START" | tee -a $LOG
echo -e "======================================================================" | tee -a $LOG

SV_END_CHROM=$(grep "^$CHROM\s$SV_START" $VCF_FILE | sed -e "s/.*CHR2=\([a-Z0-9]*\).*/\1/g")
if [ ! "$SV_END_CHROM" == "$CHROM" ]
    then
    echo "[ERROR sweep_variant.sh] polishing variants where CHR != CHR2 ($CHROM != $SV_END_CHROM) is currently not supprted."
    exit 1
fi
SV_END=$(grep "^$CHROM\s$SV_START" $VCF_FILE | sed -e "s/.*END=\([0-9]*\).*/\1/g" )
SV_TYPE=$(grep "^$CHROM\s$SV_START" $VCF_FILE |cut -f 5)
SV_LEN=$(grep "^$CHROM\s$SV_START" $VCF_FILE | sed -e "s/.*SVLEN=\([0-9]*\).*/\1/g" )
if [ -z $SV_LEN ]
    then
    SV_LEN=$(($SV_END-$SV_START))
fi
START=$(($SV_START-$REGION_AROUND_VARIANT > 0 ? $SV_START-$REGION_AROUND_VARIANT : 0))
END=$(($SV_END+$REGION_AROUND_VARIANT))

# -----------------------------------------------------------------------------
#                       Extract ONT reads
# -----------------------------------------------------------------------------

# extract long reads that span either the start or end position of the SV or both
# exclude secondary alignments, low quality or duplicates but inlcude supplementary allignments (-F 1792)
echo  "### Extract ONT reads for left break point. Region $CHROM:$(($SV_START-50))-$(($SV_START+50))" | tee -a $LOG
samtools view -F 1792 "$ONT_BAM_FILE" "$CHROM:$(($SV_START-50))-$(($SV_START+50))" > "chosen.long.reads.sam" 2>> $LOG
echo  "### Extract ONT reads for right break point. Region $CHROM:$(($SV_END-50))-$(($SV_END+50))" | tee -a $LOG
samtools view -F 1792 "$ONT_BAM_FILE" "$CHROM:$(($SV_END-50))-$(($SV_END+50))" >> "chosen.long.reads.sam" 2>> $LOG
sort "chosen.long.reads.sam" | uniq > "chosen.long.reads.sorted_by_name.sam" 2>> $LOG

# -----------------------------------------------------------------------------
#                       Extract Illumina reads
# -----------------------------------------------------------------------------

# Illumina reads should be taken +-$REGION_AROUND_VARIANTS around the breakpoints.
# If start end end point of the SV are further apart then $REGION_AROUND_VARIANT*2 bp's
# then there are two regions to extract reads from.
# Note: To limit computation time, at most 2000 reads are allowed (~150x Coverage).
#       In case of more extracted reads, 'shuf' is used for random down-sampling.
if [ $((SV_END - SV_START)) -ge 401 ]
    then
    echo  "### Extract Illumina short reads for left break point -> region $CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_START + $REGION_AROUND_VARIANT))" | tee -a $LOG
    samtools view "$BAMFILE_ILLUMINA" "$CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_START + $REGION_AROUND_VARIANT))" > "chosen.short.reads.sam" 2>> $LOG
    echo  "### Extract Illumina short reads for right break point -> region $CHROM:$((SV_END - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" | tee -a $LOG
    samtools view "$BAMFILE_ILLUMINA" "$CHROM:$((SV_END - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" >> "chosen.short.reads.sam" 2>> $LOG
    shuf -n 2000 "chosen.short.reads.sam" | sort > "chosen.short.reads.sorted_by_name.sam"
else
    echo  "### Extract Illumina short reads for region $CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" | tee -a $LOG
    samtools view "$BAMFILE_ILLUMINA" "$CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" | shuf -n 2000 | sort > "chosen.short.reads.sorted_by_name.sam" 2>> $LOG
fi

# There might be reads where its mate is not in the bam file now because it mapped to a different chromosome. Those should be included
# Filtering for aberrant read pairs means the flags 2,4 and 8 are NOT SET.
cat "$BAMFILE_ILLUMINA_HEADER" "chosen.short.reads.sorted_by_name.sam" | samtools view -F 14 > "aberrant-short-reads.sam" 2>> $LOG
if [ -s "aberrant-short-reads.sam" ] # if file not empty
    then
    echo "--- ATTENTION: Found $(wc -l aberrant-short-reads.sam) aberrant read pairs. Start hunting their mates.." | tee -a $LOG
    cut -f 1,3,7,8 "aberrant-short-reads.sam" |  awk '{ if ($3 == "=") {print $1,$2,$4;} else {print $1,$3,$4;} }' > matesToHunt.txt 2>> $LOG
    samtools view -H "$BAMFILE_ILLUMINA" -o "dummy.bam" 2>> $LOG # mate hunter needs a bam file to read and copy from which is not neeed here
    # samtools view -H "$BAMFILE_ILLUMINA_HEADER" -o "mates.bam"
    echo "--- -> Mate Hunting.." && "$EXEC_DIR/numtMateHunt" matesToHunt.txt "$BAMFILE_ILLUMINA" "$BAMFILE_ILLUMINA.bai" "dummy.bam" "mates.bam"  | tee -a $LOG
    samtools view mates.bam >> "chosen.short.reads.sorted_by_name.sam" 2>> $LOG
    sort "chosen.short.reads.sorted_by_name.sam" | uniq > "tmp.sam" 2>> $LOG
    rm "chosen.short.reads.sorted_by_name.sam" "aberrant-short-reads.sam" "dummy.bam" 2>> $LOG # "mates.bam" "matesToHunt.txt"
    mv "tmp.sam" "chosen.short.reads.sorted_by_name.sam" 2>> $LOG
fi

cat "$BAMFILE_ILLUMINA_HEADER" "chosen.short.reads.sorted_by_name.sam" | \
    samtools view -b - | \
    samtools fastq -n -1 "illumina.paired1.fastq" -2 "illumina.paired2.fastq" -  2>> $LOG

# -----------------------------------------------------------------------------
#                             Compute conensus
# -----------------------------------------------------------------------------
grep "^$CHROM\s$SV_START" $VCF_FILE > "sv.$CHROM.$SV_START.vcf"
$EXEC_DIR/build_consensus_from_all_reads "chosen.long.reads.sorted_by_name.sam" "sv.$CHROM.$SV_START.vcf"

# Output: consensus.fa

# -----------------------------------------------------------------------------
#                              Polish conensus
# -----------------------------------------------------------------------------
echo "### Perform Pilon polishing rounds until there is no change anymore" >> $LOG
echo "Note: change means that the pilon output contains no 'Corrected X snps/indels/..'" >> $LOG
echo "      with X > 0. And no 'BreakFix' in the output." >> $LOG

REF="consensus.fa" # The current "genome" reference (read to polish)
ILLUMINA1="illumina.paired1.fastq"
ILLUMINA2="illumina.paired2.fastq"

$EXEC_DIR/sweep_sequence.sh $CONFIG $REF $ILLUMINA1 $ILLUMINA2

# Output: polished.fasta

# -----------------------------------------------------------------------------
#                       Align polished read to reference
# -----------------------------------------------------------------------------

# cat variant information into id for easier match up
echo -e ">final_""$CHROM""_""$SV_START""_""$(cut -f 3 sv.$CHROM.$SV_START.vcf)" > "final.fa"
grep -v "^>" "polished.fasta" | tr -d '\n' > "middle.txt"
samtools faidx "$OT/hg38.fa" "$CHROM:$(($START-5000))-$(($START))" | grep -v "^>" | tr -d '\n' > "left-flank.txt"
samtools faidx "$OT/hg38.fa" "$CHROM:$(($END))-$(($END+5000))"     | grep -v "^>" | tr -d '\n' > "right-flank.txt"
cat "left-flank.txt" "middle.txt" "right-flank.txt" | fold -w 70 >> "final.fa"
echo -e "\n" >> "final.fa"

# echo "## align polished sequence with minimap"
# $OT/programs/minimap2-2.3_x64-linux/minimap2 -ax map-ont "$OT/hg38.mmi" "final.fa" > "final.sam"
#
# echo "## Evaluate Final Mapping"
# $EXEC_DIR/evaluate_final_mapping "final.sam"
#
#
# cat "$BAMFILE_ONT_HEADER" > polished.hg38.sam
#
# grep -v "^@" "polished.sam"                 | awk "BEGIN{OFS=\"\t\";} {\$1=\"polished.ngmlr\"; print \$0;}" >> polished.hg38.sam
# grep -v "^@" "polished.minimap.asm.sam"     | awk "BEGIN{OFS=\"\t\";} {\$1=\"polished.minimap.asm\"; print \$0;}" >> polished.hg38.sam
# grep -v "^@" "polished.minimap.default.sam" | awk "BEGIN{OFS=\"\t\";} {\$1=\"polished.minimap.default\"; print \$0;}" >> polished.hg38.sam
# grep -v "^@" "polished.minimap.ont.sam"     | awk "BEGIN{OFS=\"\t\";} {\$1=\"polished.minimap.ont\"; print \$0;}" >> polished.hg38.sam
#
# samtools sort polished.hg38.sam > polished.hg38.bam
# samtools index polished.hg38.bam


echo "--------------------------------------------------------------------------------" >> $LOG
echo "                                     DONE" | tee -a $LOG
echo "--------------------------------------------------------------------------------" | tee -a $LOG
