NanoSweep
=========

This pipeline polishes deletions and isnertions called on long read data (ONT) by applying [Pilon](https://github.com/broadinstitute/pilon).

Installation
------------

For installation, simply clone the repo and use make to compile any utilities.

~~~~
~$ git clone git@lsource2.decode.is:svenjam/NanoSweep.git
~$ cd NanoSweep
~$ make
~~~~

Note: You need a compiler that supports c++14.

Dependencies
------------

The pipeline depends on the following (unix) programs:

* _grep_, _sed_, _awk_, _cat_, _cut_, _sort_ and _shuf_. Well, and _echo_.
* _samtools_ (1.5, using htslib 1.5)
* _bwa (index/mem)_
* A long read mapper of your choice. Recommended: _minimap_
* _Pilon_ (A released jar file is in this repo ready to use)

- - - -

Using NanoSweep
---------------

This respository contains several executables that can be used to either sweep a whole sample (_sweep_sample.sh_), a single variant (_sweep_variant.sh_) or even a single sequence (_sweep_sequence.sh_).

For every (sweep)script you need to provide a config file.

### The config file

You don't need to copy the text below but just execute _sweep\_variant.sh_ without any arguments and a sample config file is dropped into the current directory.

~~~~
~$ ./sweep_sample
[ERROR - sweep_sample.sh] Config File must be provided as first argument!
I just dropped a sample config file 'polishing.config' in the  current directory for you to alter.
~~~~


```bash

# Polishing Pipeline Config File
# ==============================

# Input parameters for polishing
# -------------------------------
export OUTPUT_DIR="/path/to/directory/toStore/the/output"
export REFERENCE="/path/to/reference/genome/file.fa"
export ONT_BAM_FILE="/path/to/long/reads/that/shall/be/polished/data.bam"
export BAMFILE_ILLUMINA="/path/to/accurate/short/read/to/be/used/for/polishiing/data.bam"
export VCF_FILE="/path/to/variants/or/regions/you/want/to/have/polished/data.vcf"
export REGION_AROUND_VARIANT=400 #default
export LOG=polishing.log #default log file name

# Third party programs
# -------------------------------
# if you want to use another mapper please adjust the mapper call
# WIHOUT changing the variables MAPPER_IN_READS and MAPPER_OUT but
# by placing them and the appropiate places.
export MAPPER_CALL="/odinn/tmp/svenjam/programs/minimap2-2.3_x64-linux/minimap2 -ax map-ont /odinn/tmp/svenjam/hg38.mmi $MAPPER_IN_READS > $MAPPER_OUT"
export PILON_EXECUTABLE="$EXEC_DIR/pilon-1.22.jar"

# DO NOT TOUCH
# -------------------------------
# The following preparations are neccessary once, regardless if executing
# sweep_sample/variant/sequence so the changes are made here for preparation
if [ ! -d "$OUTPUT_DIR" ]
    then
    mkdir $OUTPUT_DIR
fi

export BAMFILE_ILLUMINA_HEADER="$OUTPUT_DIR/BAMFILE_ILLUMINA.header"
samtools view -H "$BAMFILE_ILLUMINA" -o "$BAMFILE_ILLUMINA_HEADER"

export BAMFILE_ONT_HEADER="$OUTPUT_DIR/BAMFILE_ONT.header"
samtools view -H "$ONT_BAM_FILE" -o "$BAMFILE_ONT_HEADER"

```

Simply update the variables to point to your data files.

> IMPORTANT: Use absolute paths!

Some Explanations:

* VCF_FILE
    This file contains all the variants that will be polished in the run.
    If you want to accelarate polishing on a large vcf file, you might want to split the file into several small ones and create a config fie for each.

* REGION_AROUND_VARIANT
    This is a parameter that fixes the size of the flanking area (in bp's, e.g. 50) around each breakpoint (start/end) of a variant.

    ```
                    start x             end y
    ------------------|------------------|----------------
               !------------!      !------------!
              x-50        x+50    y-50         y+50
    ```

* MAPPER_CALL
    This variable is very sensible to changes so touch it with care.
    The polishing pipeline currently requires a remapping of the polished consensus sequence to the reference (A simple pairwise alignment does not work well with standard DP). Therefore a long read mapper is needed after all consensus sequences are collected in a file called _final-all.fa_.
    Change the mapper call string to any mapper of your choice and adjust the commandline arguments **without** changing the variables _$MAPPER_IN_READS_ and _$MAPPER_OUT_. Thos will be replaced in the script such that the mapper gets the correct files for mapping.
    For example if you want to use NGMLR you might want to change the call to:

    ~~~~
    export MAPPER_CALL="/odinn/tmp/svenjam/programs/ngmlr-0.2.6/ngmlr -r $REFERENCE -x ont -q $MAPPER_IN_READS -o $MAPPER_OUT"
    ~~~~

..and don't touch the DO NOT TOUCH part. :-)

### Sweep A Genome

To polish a whole genome you mainly need to adjust the config file and give it to the script _sweep\_sample.sh_ as a first argument.

~~~~
~$ ./sweep_sample polishing.config
~~~~

This will automatically check the programs needed and launch _sweep\_variant.sh_ for each line in the vcf file (VCF_FILE in config).

You will see the progress of the sweeping on your terminal. There is no need to pipe the output because it will also appear (in more detail) in the logfile (LOG_FILE in config) in your output directory (OUTPUT_DIR in config).

### Sweep A Variant

You can also sweep a single variant by giving the script _sweep\_variant.sh_ the config file, as well as the chromosome/reference name and the start position of your variant of choice. The script will find the corresponding variant in the vcf file (VCF_FILE in config) with the help of grep.

~~~~
~$ ./sweep_variant polishing.config chr1 993459
~~~~

Note: Of course you can also create a vcf file of only one variant and run _sweep\_sample.sh_.

### Sweep A Sequence

In case you want to only polish a specific sequence for another application, you can use the _sweep\_sequence.sh_ script.

This a needs a little more prior work though: The script expects a sequence to polish in fasta format and two files with illumina paired reads in fasta/fastq format that you wish to polish with. For ecample you could extract your sequence and reads of choice using _samtools view_:

~~~~
~$ samtools view myAssemblySequence.fa "seq:100-900" > sequence.fa
~$ samtools view -hb myBamFile.bm "ref:5500-6300" | \
samtools fastq - -1 illumina1.fq -2 illumina2.fq
~$ ./sweep_sequence polishing.config sequence.fa illumina1.fq illumina2.fq
~~~~

### Contact
If you have any questions write me a mail: svenja.mehringer[AT]gmail.com
