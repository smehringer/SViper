#!/bin/sh

export EXEC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -ne 4 ]; then
    echo "usage: CONFIG.config REF_SEQUENCE.fa SHORT_READS_1.fastq SHORT_READS_.fasq"
    exit 1
fi

CONFIG=$1
REF=$2
ILLUMINA1=$3
ILLUMINA2=$4

if [ "${CONFIG:0:1}" != "/" ] # no absolute path
    then
    CONFIG="$(pwd)/$CONFIG"
fi

if [ ! -f "$CONFIG" ]
    then
    echo "[ERROR - process_sequence.sh] Could not find file $CONFIG."
    exit 1
fi

source "$CONFIG"

ROUND=1                    # Tracks the current polishing round
PIL="pilon$ROUND"          # The current pilon filename prefix

while true
do
    echo -e "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ROUND $ROUND " >> $LOG
    echo -e "### BWA Indexing. " >> $LOG
    bwa index $REF &>> $LOG

    echo -3 "\ņ### BWA mapping of illumina reads. " >> $LOG
    bwa mem $REF $ILLUMINA1 $ILLUMINA2 2>> $LOG 1> mapped.illumina.reads.sam
    samtools sort mapped.illumina.reads.sam > mapped.illumina.reads.bam 2>> $LOG
    samtools index mapped.illumina.reads.bam 2>> $LOG

    echo -e "ņ### Pilon Piloshing. " >> $LOG # output: pilon.fasta
    java -jar $PILON_EXECUTABLE --genome $REF --frags mapped.illumina.reads.bam > out.pilon
    cat out.pilon >> $LOG

    CORRECTIONS=$(cat out.pilon | grep "Corrected" | grep -E "[123456789]") || true # ignore grep error
    BREAKFIXES=$(cat out.pilon | grep "BreakFix") || true                           # ignore grep error
    echo -e "\nCORRECTIONS: $CORRECTIONS" | fold -w 80 >> $LOG
    echo -e "\nBREAKFIXES: $BREAKFIXES" | fold -w 80 >> $LOG

    # cleaning up a little
    rm $REF.* # remove bwa index files
    #rm mapped.reads* # remove mapping to "the-chosen-one"
    mv pilon.fasta $PIL.fasta

    # check if pilon did some changes
    if [ -z "$CORRECTIONS" ] && [ -z "$BREAKFIXES" ]
    then # grep results are empty => no changes
        echo "~~~~~~~~~~~~~~~~~~~~~~~ No more changes done. Exiting loop. ~~~~~~~~~~~~~~~~~~~" >> $LOG
        break
    fi

    # renew variables for next round
    let ROUND=$ROUND+1
    REF=$PIL.fasta
    PIL="pilon$ROUND"

    echo "Done with round $ROUND." >> $LOG

    if [ "$ROUND" -eq "20" ]
    then # stop after 20 rounds
        echo "Enough now." >> $LOG
        break
    fi
done

mv $PIL.fasta polished.fasta
