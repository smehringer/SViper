#!/bin/sh

export EXEC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -ne 1 ]; then
    echo "[ERROR - sweep_sample.sh] Config File must be provided as first argument!"

    if [ ! -f "polishing.config" ] # if default file does not already exist
        then
        echo '# Polishing Pipeline Config File'   > polishing.config
        echo '# =============================='  >> polishing.config
        echo ''                                  >> polishing.config
        echo '# Input parameters for polishing'  >> polishing.config
        echo '# -------------------------------' >> polishing.config
        echo 'export OUTPUT_DIR="/path/to/directory/toStore/the/output"' >> polishing.config
        echo 'export REFERENCE="/path/to/reference/genome/file.fa"' >> polishing.config
        echo 'export ONT_BAM_FILE="/path/to/long/reads/that/shall/be/polished/data.bam"' >> polishing.config
        echo 'export BAMFILE_ILLUMINA="/path/to/accurate/short/read/to/be/used/for/polishiing/data.bam"' >> polishing.config
        echo 'export VCF_FILE="/path/to/variants/or/regions/you/want/to/have/polished/data.vcf"' >> polishing.config
        echo 'export REGION_AROUND_VARIANT=400 #default' >> polishing.config
        echo 'export LOG=polishing.log #default log file name' >> polishing.config
        echo ''                                  >> polishing.config
        echo '# Third party programs'            >> polishing.config
        echo '# -------------------------------' >> polishing.config
        echo '# if you want to use another mapper please adjust the mapper call' >> polishing.config
        echo '# WIHOUT changing the variables MAPPER_IN_READS and MAPPER_OUT but' >> polishing.config
        echo '# by placing them and the appropiate places.' >> polishing.config
        echo 'export MAPPER_CALL="/odinn/tmp/svenjam/programs/minimap2-2.3_x64-linux/minimap2 -ax map-ont /odinn/tmp/svenjam/hg38.mmi $MAPPER_IN_READS > $MAPPER_OUT"' >> polishing.config
        echo 'export PILON_EXECUTABLE="$EXEC_DIR/pilon-1.22.jar"'                           >> polishing.config
        echo ''                                                                             >> polishing.config
        echo '# DO NOT TOUCH'                                                               >> polishing.config
        echo '# -------------------------------'                                            >> polishing.config
        echo '# The following preparations are neccessary once, regardless if executing'    >> polishing.config
        echo '# sweep_sample/variant/sequence so the changes are made here for preparation' >> polishing.config
        echo 'if [ ! -d "$OUTPUT_DIR" ]'                                                    >> polishing.config
        echo '    then'                                                                     >> polishing.config
        echo '    mkdir $OUTPUT_DIR'                                                        >> polishing.config
        echo 'fi'                                                                           >> polishing.config
        echo ''                                                                             >> polishing.config
        echo 'export BAMFILE_ILLUMINA_HEADER="$OUTPUT_DIR/BAMFILE_ILLUMINA.header"'         >> polishing.config
        echo 'samtools view -H "$BAMFILE_ILLUMINA" -o "$BAMFILE_ILLUMINA_HEADER"'           >> polishing.config
        echo ''                                                                             >> polishing.config
        echo 'export BAMFILE_ONT_HEADER="$OUTPUT_DIR/BAMFILE_ONT.header"'                   >> polishing.config
        echo 'samtools view -H "$ONT_BAM_FILE" -o "$BAMFILE_ONT_HEADER"'                    >> polishing.config

        echo "I just dropped a sample config file 'polishing.config' in the "\
             "current directory for you to alter."
    fi

    exit 1
fi

CONFIG=$1
if [ "${CONFIG:0:1}" != "/" ] # no absolute path
    then
    CONFIG="$(pwd)/$CONFIG"
fi

if [ ! -f "$CONFIG" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find file $CONFIG."
    exit 1
fi

# File to source global variables
# prior to sourcing, set variables $MAPPER_IN_REF $MAPPER_IN_READS and $MAPPER_OUT
# such that the mapping call works correctly
export MAPPER_IN_READS="final-all.fa"
export MAPPER_OUT="final-all.sam"
source $CONFIG

# ------------------------------------------------------------------------------
#                            File existence checks
# ------------------------------------------------------------------------------

# Check output dir
# -----------------
if [ ! -d "$OUTPUT_DIR" ]
    then
    mkdir $OUTPUT_DIR
fi

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

if [ ! -f "$PILON_EXECUTABLE" ]
    then
    echo "[ERROR - sweep_sample.sh] Could not find pilon executable $PILON_EXECUTABLE (CPP)."\
         " Please ensure correct executable path or change the config file."
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
#                            Polishing
# ------------------------------------------------------------------------------

cd $OUTPUT_DIR

DATE=$(date +"%d-%B-%Y")
LOG="$(pwd)/polishing.log"    # logs all the programm output
echo $DATE > $LOG      # initialize logfile

echo "----------------------------------------------------------------------" | tee -a $LOG
echo "START polishing variants in $VCF_FILE" | tee -a $LOG
echo "----------------------------------------------------------------------" | tee -a $LOG

grep -v "^#" "$VCF_FILE" > variants.vcf 2>> $LOG # subtract header for line by line processing

COUNT=0
PREV_POS="-1"
while IFS= read -r line
	do
	SECONDS=0
	CHROM=$(echo "$line" | cut -f 1)
	POS=$(echo "$line" | cut -f 2)
	ALT=$(echo "$line" | cut -f 5)

	if [ "$ALT" != "<INS>" ] && [ "$ALT" != "<DEL>" ]
		then
		echo -e "\e[93m[  SKIP  ] \e[39mSkip SV $CHROM:$POS because $ALT is currently not supported." | tee -a $LOG
		continue
	fi

    if [ -d "sv.$CHROM.$POS" ]
        then
        echo -e "\e[93m[  SKIP  ] \e[39mSkip SV $CHROM:$POS because the directory sv.$CHROM.$POS already exists." | tee -a $LOG
        continue
    else
        mkdir "sv.$CHROM.$POS"
        cd "sv.$CHROM.$POS"
    fi

    # store reported supporting reads for downstream analysis like scoring
	RNAMES=$(echo $line | sed -e "s/.*RNAMES=\([^\;]*\);.*/\1/g")
	echo $RNAMES | tr "," "\n" > "supporting_reads_ids.txt"

	# do the polishing
	$EXEC_DIR/sweep_variant.sh "$CONFIG" "$CHROM" "$POS" &>> $LOG

	if [ $? -eq 0 ]
		then
		echo -e "\e[92m[ SUCCESS ] \e[39mSuccessfully polished SV $CHROM:$POS \t \e[93m$ALT\e[39m \t[~ $(($SECONDS / 60))m $(($SECONDS % 60))s]" | tee -a $LOG
	else
		echo -e "\e[91m[  ERROR  ] \e[39mError while polishing SV $CHROM:$POS \t \e[93m$ALT\e[39m \t[~ $(($SECONDS / 60))m $(($SECONDS % 60))s]" | tee -a $LOG
	fi

	cd ..

	let COUNT=COUNT+1

	PREV_POS=$POS #update

done < variants.vcf
echo "----------------------------------------------------------------------"       | tee -a $LOG

cat sv*/final.fa > "final-all.fa" 2>> $LOG

if [ ! -f "final-all.fa" ] || [ -s "final-all.fa" ] # does not exists or is empty
    then
    echo "ATTENTION: No variants were polished."
    echo "----------------------------------------------------------------------"       | tee -a $LOG
    echo "                                 END"                                         | tee -a $LOG
    echo "----------------------------------------------------------------------"       | tee -a $LOG
fi

# ------------------------------------------------------------------------------
echo "### Remapping. "                                                              | tee -a $LOG
# ------------------------------------------------------------------------------
echo "## align polished sequence by calling the following command:"                 | tee -a $LOG
echo "## $MAPPER_CALL"                                                              | tee -a $LOG
SECONDS=0
eval "$MAPPER_CALL 2>> $LOG"
samtools sort "final-all.sam" -o "final-all.bam" 2>> $LOG
samtools index "final-all.bam" 2>> $LOG
echo -e "\t\t\t\t\t\t\t[~ $(($SECONDS / 60))m $(($SECONDS % 60))s]"                 | tee -a $LOG

# ------------------------------------------------------------------------------
echo "### Create NameSorted sam file for evaluation. "                              | tee -a $LOG
# ------------------------------------------------------------------------------
SECONDS=0
grep -v "^@" "final-all.sam" | sort > "final-all.sortedByName.sam" 2>> $LOG
echo -e "\t\t\t\t\t\t\t[~ $(($SECONDS / 60))m $(($SECONDS % 60))s]"                 | tee -a $LOG

# ------------------------------------------------------------------------------
echo "### Evaluation of variants. "                                                 | tee -a $LOG
# ------------------------------------------------------------------------------
SECONDS=0
$EXEC_DIR/utilities/evaluate_final_mapping "final-all.sortedByName.sam" "$VCF_FILE" > "$VCF_FILE.polished.vcf" 2>> $LOG
echo -e "\t\t\t\t\t\t\t[~ $(($SECONDS / 60))m $(($SECONDS % 60))s]"                 | tee -a $LOG

echo "----------------------------------------------------------------------"       | tee -a $LOG
echo "                                 END"                                         | tee -a $LOG
echo "----------------------------------------------------------------------"       | tee -a $LOG
