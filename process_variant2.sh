#!/bin/sh

# -----------------------------------------------------------------------------
#                                       Info
# -----------------------------------------------------------------------------
#
# This script needs to have the vcf file sorted by position!
#

set -e

if [ "$#" -ne 4 ]; then
    echo "usage: CONFIG CHROM SV_START SV_BEFORE"
    exit 1
fi

CONFIG=$1
CHROM=$2
SV_START=$3
SV_BEFORE=$4

source "$(pwd)/$CONFIG"

DATE=$(date +"%d-%B-%Y")
LOG="$(pwd)/sv.$CHROM.$SV_START.log" # logs all the programm output
echo $DATE > $LOG                    # initialize logfile

BAMFILE_ILLUMINA_HEADER="$POLISHING_DIR/BAMFILE_ILLUMINA.header"
samtools view -H "$BAMFILE_ILLUMINA" -o "$BAMFILE_ILLUMINA_HEADER"

BAMFILE_ONT_HEADER="$POLISHING_DIR/BAMFILE_ONT.header"
samtools view -H "$ONT_BAM_FILE" -o "$BAMFILE_ONT_HEADER"


mkdir "polished-reads"
cd "polished-reads"

echo "" | tee -a $LOG
echo -e "======================================================================" | tee -a $LOG
echo -e "                       POLISHING VARIANT $CHROM:$SV_START" | tee -a $LOG
echo -e "======================================================================" | tee -a $LOG

SV_END_CHROM=$(grep "^$CHROM\s$SV_START" $VCF_FILE | sed -e "s/.*CHR2=\([a-Z0-9]*\).*/\1/g")
if [ ! "$SV_END_CHROM" == "$CHROM" ]
    then
    echo "[ERROR process_variant.sh] polishing variants where CHR != CHR2 ($CHROM != $SV_END_CHROM) is currently not supprted."
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
# exclude secondary and supplementary allignments (-F 2304)
echo  "### Extract ONT reads for left break point. Region $CHROM:$SV_START-$SV_START" | tee -a $LOG
samtools view -F 2304 "$ONT_BAM_FILE" "$CHROM:$SV_START-$SV_START" > "chosen.long.reads.sam" 2>> $LOG
echo  "### Extract ONT reads for right break point. Region $CHROM:$SV_END-$SV_END" | tee -a $LOG
samtools view -F 2304 "$ONT_BAM_FILE" "$CHROM:$SV_END-$SV_END" >> "chosen.long.reads.sam" 2>> $LOG
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
#                       Process every single read
# -----------------------------------------------------------------------------

echo "### Get supporting reads by cigar string."
python $OU/polishing/get_supporting_reads.py "chosen.long.reads.sorted_by_name.sam" "$SV_TYPE" "$SV_START" "$SV_END" "$SV_LEN" | tee -a $LOG

echo "### Now process all supporting $(wc -l chosen.long.reads.sorted_by_name.sam.supporting.sam) reads... " | tee -a $LOG
while IFS= read -r line
do
    READ_NAME=$(echo "$line" | cut -f 1)

    echo ""  >> $LOG
    echo -e "\e[1m $READ_NAME" | tee -a $LOG
    echo -e "--------------------------------------------------------------------------------\e[0m" >> $LOG

    if [ ! -d "$READ_NAME" ]; then
        mkdir "$READ_NAME"
    fi

    cd "$READ_NAME"

    echo "### Crop FA region in read that mapps to reference region [$START-$END]" >> $LOG
    echo ">$READ_NAME" > "$READ_NAME.fa"
    READ_MAPPPING_POS=$(echo "$line" | cut -f 4)

    # in case the read was involved in a variant polishing before we want to
    # continue polishing the polished read
    OFFSET=0
    PREV_DIR="../../../sv.$CHROM.$SV_BEFORE/polished-reads/$READ_NAME"
    if [ -f "$PREV_DIR/polished-read.fasta" ]
        then
        echo "ATTENTION: read was already partly polished in sv.$CHROM.$SV_BEFORE." >> $LOG
        OFFSET=$(grep "^>" "$PREV_DIR/polished-read.fasta" | cut -f 2)
        READ_SEQ=$(grep -v "^>" "$PREV_DIR/polished-read.fasta" | tr -d '\n')
        mv "$PREV_DIR/polished-read.fasta" "$PREV_DIR/OLD-polished-read.fasta" 2>> $LOG
        echo "See sv.$CHROM.$POS/polished-reads/$READ_NAME for polished-read.fasta" > "$PREV_DIR/MOVED_POLISHED_READ.txt"
    else
        READ_SEQ=$(echo "$line" | cut -f 10)
    fi

    echo "$EXEC_DIR/compute_region_positions_from_cigar $START $END $(echo "$line" | cut -f 4) $(echo "$line" | cut -f 6)" >> $LOG # for debugging
    READ_REGION=$($EXEC_DIR/compute_region_positions_from_cigar $START $END $(echo "$line" | cut -f 4) $(echo "$line" | cut -f 6) )
    set -- $READ_REGION
    READ_REGION_START=$(($1+$OFFSET))
    READ_REGION_LENGTH="$2"
    READ_REGION_END=$(($READ_REGION_START+$READ_REGION_LENGTH))
    echo "${READ_SEQ:$READ_REGION_START:$READ_REGION_LENGTH}" | fold -w 80 >> "$READ_NAME.fa"
    echo "--- -> Read Span: [$READ_REGION_START-$READ_REGION_END]" >> $LOG

    cd ..

done < "chosen.long.reads.sorted_by_name.sam.supporting.sam"

cat *template/*template.fa > regions.fa

$EXEC_DIR/build_consensus_2 regions.fa > consensus.fa

echo "### Perform Pilon polishing rounds until there is no change anymore" >> $LOG
echo "Note: change means that the pilon output contains no 'Corrected X snps/indels/..'" >> $LOG
echo "      with X > 0. And no 'BreakFix' in the output." >> $LOG

REF="consensus.fa" # The current "genome" reference (read to polish)
ILLUMINA1="illumina.paired1.fastq"
ILLUMINA2="illumina.paired2.fastq"

$EXEC_DIR/process_sequence.sh ../$CONFIG $REF $ILLUMINA1 $ILLUMINA2

samtools faidx $OT/hg38.fa "$CHROM:$(($START-500))-$(($END-500))" > ref-region.fa

# ------------------------------------------------------------------------------
echo "### NGLMR mapping to ref-region.fa " | tee -a $LOG
# ------------------------------------------------------------------------------
/odinn/tmp/svenjam/programs/ngmlr-0.2.6/ngmlr -r ref-region.fa -q "polished.fasta" -o "polished.sam" | tee -a $LOG
grep -v "^@" polished.sam | grep "$CHROM:" | awk -v OFS='\t' '{split($3,chr,":"); split(chr[2],pos,"-"); print $1,$2,chr[1],(pos[1]+$4-1),$5,$6,$7,$8,$9,$10,$11,$12,$13; }' > polished.hg38.sam

cat "$BAMFILE_ONT_HEADER" polished.hg38.sam | samtools sort - > "polished.hg38.bam"
samtools index "polished.hg38.bam"

echo "--------------------------------------------------------------------------------" >> $LOG
echo "                                     DONE" | tee -a $LOG
echo "--------------------------------------------------------------------------------" | tee -a $LOG

cd ..
