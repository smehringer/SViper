import re
import sys

# call by
# cat chosen.long.reads.sam | python ../get-names.py "$(grep 933969 ../IWXJYQX_sniffles_s2.vcf)"

# first argument is the sv in vcf format
if (len(sys.argv) != 6):
    print("USAGE ERROR: SAM_FILE SV_TYPE START END LENGTH")

sv_type = sys.argv[2] # <DEl> or <INS>
if (sv_type != '<DEL>' and sv_type != '<INS>'):
    sys.exit("Currently, we do not support SV type '" + sv_type +"'. Only <DEL> and <INS>")
elif (sv_type == '<DEL>'):
    sv_type = 'D'
else:
    sv_type = 'I'

sam_file_in = open(sys.argv[1], 'r')
sam_file_out = open((sys.argv[1] + ".supporting.sam"), 'w')
sv_start = int(sys.argv[3])
sv_end   = int(sys.argv[4])
sv_len   = int(sys.argv[5])

#dev_pos = bp deviation in position of the SV
if (sv_type == 'D'): # deletion
    dev_pos = sv_len - 1 # deletion must overlap
else: # Insertion
    dev_pos = 500 + sv_len # insertions are often very imprecise
dev_size = 0.8*sv_len  #bp deviation in length of the SV

# print("Searching for " + str(sv_len - dev) + "-" + str(sv_len + dev) + sv_type + " in range [" + str(sv_start - dev) + "," + str(sv_end + dev) + "]")

# pipe in reads in sam format
for line in sam_file_in:

    read = line.split()

    sv_name  = read[0]
    cigar    = read[5]

    numbers = re.split("[A-Z]",  cigar)[:-1] # without last empty field
    identis = re.split("[0-9]+", cigar)[1:]  # wihtout first empty field

    ref_pos   = int(read[3]) # mapping position
    index = 0

    # only check the cigar string within the sv region [sv_start-dev, sv_end+dev]

    # go to startof range (TODO this could be done by a binary search)
    # go to startof range (TODO this could be done by a binary search)
    while (ref_pos <= sv_start - dev_pos) & (index < len(numbers)):
        if (identis[index] != 'I') & (identis[index] != 'S'):
            ref_pos += int(numbers[index])

        index += 1

    # check within range
    while (ref_pos < sv_end + dev_pos) & (index < len(numbers)):
        if (identis[index] == sv_type): # correct sv
            if ((int(numbers[index]) <= sv_len + 2*dev_size) & (int(numbers[index]) >= dev_size)):
                print(sv_name + " ref_pos:" + str(ref_pos) + " len:" + numbers[index])
                sam_file_out.write(line)

        if (identis[index] != 'I') & (identis[index] != 'S'):
            ref_pos += int(numbers[index])

        index += 1
