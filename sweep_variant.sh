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

# ------------------------------------------------------------------------------
#                            File existence checks
# ------------------------------------------------------------------------------

# Check executable files
# -----------------------
if [ ! -f "$EXEC_DIR/sweep_variant.sh" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find executable $EXEC_DIR/sweep_variant.sh."\
         " Please ensure correct file structure of the polishing toolset."
    exit 1
fi

if [ ! -f "$EXEC_DIR/utilities/compute_region_positions_from_cigar" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find executable $EXEC_DIR/compute_region_positions_from_cigar (CPP)."\
         " Please ensure correct file structure of the polishing toolset."
    exit 1
fi

if [ ! -f "$EXEC_DIR/utilities/numtMateHunt" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find executable $EXEC_DIR/numtMateHunt (CPP)."\
         " Please ensure correct file structure of the polishing toolset."
    exit 1
fi

# Check user provided variables in config file
# --------------------------------------------
if [ ! -f "$REFERENCE" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find file $REFERENCE."\
         " Please correct filename in polishing.config."
    exit 1
fi

if [ ! -f "$REFERENCE.fai" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find index file $ONT_BAM_FILE."\
         " Please ensure that your reference file is indexed."
    exit 1
fi

if [ ! -f "$ONT_BAM_FILE" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find file $ONT_BAM_FILE."\
         " Please correct filename in polishing.config."
    exit 1
fi

if [ ! -f "$ONT_BAM_FILE.bai" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find index file $ONT_BAM_FILE.bai."\
         " Please provide a bam index file."
    exit 1
fi

if [ ! -f "$BAMFILE_ILLUMINA" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find file $BAMFILE_ILLUMINA."\
    " Please correct filename in polishing.config."
    exit 1
fi

if [ ! -f "$BAMFILE_ILLUMINA" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find index file $BAMFILE_ILLUMINA.bai"\
    " Please provide a bam index file."
    exit 1
fi

if [ ! -f "$VCF_FILE" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find file $VCF_FILE."\
    " Please correct filename in polishing.config."
    exit 1
fi

if [ -z "$REGION_AROUND_VARIANT" ]
    then
    echo "[ERROR - sweep_sample.sh] Variable REGION_AROUND_VARIANT cannot be empty."\
    " Please correct variable in polishing.config."
    exit 1
fi

# ------------------------------------------------------------------------------
#                            Start polishing variant
# ------------------------------------------------------------------------------

DATE=$(date +"%d-%B-%Y")
LOG="$(pwd)/sv.$CHROM.$SV_START.log" # logs all the programm output
echo $DATE > $LOG                    # initialize logfile

echo "" | tee -a $LOG
echo -e "======================================================================" | tee -a $LOG
echo -e "                       POLISHING VARIANT $CHROM:$SV_START" | tee -a $LOG
echo -e "======================================================================" | tee -a $LOG

SV_END_CHROM=$(grep "^$CHROM\s$SV_START" $VCF_FILE | sed -e "s/.*CHR2=\([a-zA-Z0-9]*\).*/\1/g")
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
    samtools view -F 1536 -q 2 "$BAMFILE_ILLUMINA" "$CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_START + $REGION_AROUND_VARIANT))" > "chosen.short.reads.sam" 2>> $LOG
    echo  "### Extract Illumina short reads for right break point -> region $CHROM:$((SV_END - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" | tee -a $LOG
    samtools view -F 1536 -q 2 "$BAMFILE_ILLUMINA" "$CHROM:$((SV_END - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" >> "chosen.short.reads.sam" 2>> $LOG
    shuf -n 2000 "chosen.short.reads.sam" | sort > "chosen.short.reads.sorted_by_name.sam"
else
    echo  "### Extract Illumina short reads for region $CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" | tee -a $LOG
    samtools view -F 1536 -q 2 "$BAMFILE_ILLUMINA" "$CHROM:$((SV_START - REGION_AROUND_VARIANT))-$((SV_END + $REGION_AROUND_VARIANT))" | shuf -n 2000 | sort > "chosen.short.reads.sorted_by_name.sam" 2>> $LOG
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
    echo "--- -> Mate Hunting.." && "$EXEC_DIR/utilities/numtMateHunt" matesToHunt.txt "$BAMFILE_ILLUMINA" "$BAMFILE_ILLUMINA.bai" "dummy.bam" "mates.bam"  | tee -a $LOG
    samtools view mates.bam >> "chosen.short.reads.sorted_by_name.sam" 2>> $LOG
    sort "chosen.short.reads.sorted_by_name.sam" | uniq > "tmp.sam" 2>> $LOG
    rm "chosen.short.reads.sorted_by_name.sam" "aberrant-short-reads.sam" "dummy.bam" 2>> $LOG # "mates.bam" "matesToHunt.txt"
    mv "tmp.sam" "chosen.short.reads.sorted_by_name.sam" 2>> $LOG
fi

cat "$BAMFILE_ILLUMINA_HEADER" "chosen.short.reads.sorted_by_name.sam" | \
    samtools view -b - | \
    samtools fastq -n -1 "illumina.paired1.fastq" -2 "illumina.paired2.fastq" -  2>> $LOG

# -----------------------------------------------------------------------------
#                 Compute and polish conensus with cpp programm
# -----------------------------------------------------------------------------
grep "^$CHROM\s$SV_START" $VCF_FILE > "sv.$CHROM.$SV_START.vcf"
$EXEC_DIR/utilities/polishing -v "$REFERENCE" "chosen.long.reads.sorted_by_name.sam" "sv.$CHROM.$SV_START.vcf"

# output: final.fa (which is the consensus flanked by +-5000 bo reference)

echo "--------------------------------------------------------------------------------" >> $LOG
echo "                                     DONE" | tee -a $LOG
echo "--------------------------------------------------------------------------------" | tee -a $LOG
