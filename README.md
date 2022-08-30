SViper
=======

This pipeline polishes deletions and insertions called on long read data (ONT) using short exact reads for refinement.

Installation
------------

For installation, simply clone the repo and use cmake to compile sviper or any utilities.

~~~~
~$ git clone --recursive git@github.com:smehringer/SViper.git
~$ cd SViper
~$ mkdir build
~$ cd build
~$ cmake ..
~$ make sviper
~~~~

Note: You need a compiler that supports c++14.

Dependencies
------------

* Linux operating system
* Code so far only tested for gcc 5.4.0
* C++14 support
* The SeqAn C++ library (included as a submodule, no installation required)
* CMake (minimum version 3.6)

- - - -

Using SViper
---------------

You can look at all the input requirements by calling the sviper help page:

~~~~
~$ ./sviper -h
~~~~

Examples:

Call sviper
~~~~
~$ ./sviper -s short-reads.bam -l long-reads.bam -r ref.fa -c variants.vcf -o polished_variants
~~~~
This will output a `polished_variants.vcf` file, that contains all the refined variants.

Sometimes it is helpful to look at the polished sequence, e.g. with the IGV browser.
In that case you want SViper to output the polished and aligned sequences in a bam file via the option `--output-polished-bam`:
~~~~
~$ ./sviper -s short-reads.bam -l long-reads.bam -r ref.fa -c variants.vcf -o polished_variants --output-polished-bam
~~~~

### IMPORTANT

There are several requirements for using the polishing:

1. The vcf file must be a structural variant format (tags instead of sequences, e.g. `<DEL>`). ALso the INFO field must include the END tag, giving the end position of the variant, as well as the SVLEN tag in case of insertions.
2. The bam files must be indexed.
3. The reference sequence (FASTA) must be indexed.
4. (Obvisoulsy, the bam files should correspond to the same individual/sample mapped to the given reference.)

### Input arguments:

* `-c, --candidate-vcf` The candidate vcf file
    This file contains all the variants that will be polished in the run.
    If you want to accelarate polishing on a large vcf file, you might want to split the file into several small ones and call sviper each one seperately.

* `-s, --short-read-bam` (BAM_FILE)
          The indexed bam file containing short used for polishing at variant sites.

* `-l, --long-read-bam` (BAM_FILE)
          The indexed bam file containing long reads to be polished at variant sites.

* `-r, --reference` (FA_FILE)
          The indexed (fai) reference file.

* `-o, --output-prefix` (PREFIX)
          A name for the output files. The current output is a log file and vcf file, that contains the polished
          sequences for each variant.

* `-g, --log-file` (TXT_FILE)
          A filename for the log file. Default: polishing.log.

* `-k, --flanking-region` (INT)
          The flanking region in bp's around a breakpoint to be considered for polishing In range [50..1000]. Default: 400.

            ~~~~
                            start x             end y
            ------------------|------------------|----------------
                       !------------!      !------------!
                     x-400       x+400   y-400        y+400
            ~~~~

* `-x, --coverage-short-reads` (INT)
          The original short read mean coverage. This value is used to restrict short read coverage on extraction to
          avoid mapping bias

* `--output-polished-bam` (FLAG)
          For debugging or manual inspection the polished reads can be written to a file.

### VCF Output - the FAIL codes explained

* **FAIL1**: There are no long reads at all in the reference region indicated by the variant. Neither at the start location (nor at the end of the deletion).
* **FAIL2**: There are no supporting long reads for a specific variant. This means that if you are for example searching for a deletion of size 300 hundred, there are no long reads that have a deletion around the location of interest or the deletion is too large or too small.
* **FAIL3**: After selection supporting long reads, those reads are cropped around the potential variant. For a good consensus it is necessary that the variant does not lie at the end of a long read. Cropped reads (regions) that are too small because the flanking region is not large enough are discarded. If all reads have been discarded, `FAIL3` appears.
* **FAIL4**: There are not enough short reads in the area of the variant. SViper requires at least 20 short reads to be found in the area of the variant which corresponds to a short read median coverage of about 10.
* **FAIL5**: After polishing the long-read-consensus with the short reads, a variant similar to the one of interest could not be reconstructed. 

Note that a good start to check for issues in your data, interesting variants, or possible bugs in SViper is to use the [IGV Browser](https://software.broadinstitute.org/software/igv/). Upload your short and long read BAM files and check the variant locations of interest. With the option `--output-polished-bam` you can even output the polished reads and load them in IGV too. This could look like the following:

[sviper_IGV.pdf](https://github.com/smehringer/SViper/files/9451618/sviper_IGV.pdf)

### Utilities

There are some utilities that come in handy when you want to look at specific variants and understand the substeps of polishing (or the preparation). You can make the utilities with:

~~~~
~$ make utilities
~~~~

The executables can then be found in the `utilities` directory.

Some short explanation for each:

* `utilities_merge_split_alignments`
    This little program will take a bam/sam file as first argument. It will then output a file called 'your-input-filename.merged.sam' which includes all reads but all supplementary alignments are now wither  merged to it's primary or discarded. **IMPORTANT: The input BAM file needs to be sorted by name!** This is especially useful if you want to look at the sviper-output-bam file in the IGV browser (the merging is done automatically when evaluating).

* `utilities_get_supporting_reads`
    This little program will take a bam/sam file as first argument and a vcf file as second (Note that it only reads the first variant because it is made for one variant only). It will then output each reads name and wether the read is supporting or not.

* `utilities_compare_vcf`
    This little program will take a vcf file as first argument and a TRUTH vcf file as second. This script is not debugged enough so you might run into errors. If so please report them to me and I will try to fix them as soon as possible.


### Contact
If you have any questions don't hesitate to contact me: svenja.mehringer[AT]gmail.com
