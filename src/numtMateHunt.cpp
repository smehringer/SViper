// slightly altered from Snædis
#include <iostream> 
#include <seqan/bam_io.h>
#include <map>
#include <seqan/sequence.h>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>

using namespace std;
using namespace seqan;

struct FindAndPrintMateRet {
    String<BamAlignmentRecord> readsToPrint;
    Pair<unsigned> foundAndSearched;
} ;
FindAndPrintMateRet findAndPrintMate(Pair<string, unsigned> chrAndPos, 
                                     vector<string> readNames,
                                     CharString bamPathIn,
                                     BamIndex<Bai> const & baiIndex)
{
    //cout << "processing " << chrAndPos.i1 << "\t" << chrAndPos.i2 << endl;
    FindAndPrintMateRet returnValue;
    returnValue.foundAndSearched = Pair<int>(0,readNames.size());
    CharString rName = chrAndPos.i1;
    BamFileIn bigBam;
    if (!open(bigBam, toCString(bamPathIn)))
    {
        std::cerr << "ERROR: Could not open " << bamPathIn << std::endl;
        return returnValue;
    }
    BamHeader header;
    readHeader(header, bigBam);

    // Translate from reference name to rID.
    int rID = 0;
    if (!getIdByName(rID, contigNamesCache(context(bigBam)), rName))
    {
        std::cerr << "ERROR: Reference sequence named " << rName << " not known.\n";
        return returnValue;
    }
    int beginPos = chrAndPos.i2-1;
    bool hasAlignments = false;
    if (!jumpToRegion(bigBam, hasAlignments, rID, beginPos, beginPos+1, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << beginPos << "\n";
        return returnValue;
    }
    if (!hasAlignments)
    {
        cout << "No alignments for: " << chrAndPos.i1 << "\t" << chrAndPos.i2 << endl;
        return returnValue;  // No alignments here.
    }
    BamAlignmentRecord record;
    vector<string>::iterator low;
    while (!atEnd(bigBam))
    {
        readRecord(record, bigBam);
        low = lower_bound(readNames.begin(), readNames.end(), toCString(record.qName));
        if (low != readNames.end())
        {
            if ((*low).compare(toCString(record.qName))==0)
            {
                appendValue(returnValue.readsToPrint,record);
                ++returnValue.foundAndSearched.i1;
            }
        }
        if (record.beginPos > beginPos+1)
        {
            break;
        }
    }
    return returnValue;
}

int main(int argc, char const ** argv)
{
    if (argc != 6)
    {
        cerr << "USAGE: " << argv[0] << " readsMissingMate bigBam.bam bigBamIndex.bai mt.bam newMt.bam\n";
        return 1;
    }
    //Open bamfile to extend and bamFile to write extension to
    BamFileIn oldBam;
    if (!open(oldBam, argv[4]))
    {
        std::cerr << "ERROR: Could not open " << argv[4] << std::endl;
        return 1;
    }
    BamHeader header;
    readHeader(header, oldBam);
    BamFileOut extendedBam(context(oldBam), argv[5]);
    writeHeader(extendedBam, header);
    BamAlignmentRecord record;
    while (!atEnd(oldBam))
    {
        readRecord(record, oldBam);
        writeRecord(extendedBam, record);
    }
    //List of read names and locations to find mates for.
    ifstream missingMateFile(argv[1]);
    map<Pair<string, unsigned> , vector<string> > chrAndPosToNames;
    string readName, chr;
    unsigned pos;
    while (!missingMateFile.eof())
    {
        missingMateFile >> readName;
        missingMateFile >> chr;
        missingMateFile >> pos;
        chrAndPosToNames[Pair<string, unsigned>(chr,pos)].push_back(readName);
    }
    CharString bamPathIn = argv[2];
    CharString baiIndexPath = argv[3];
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(baiIndexPath)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiIndexPath << "\n";
        return 1;
    }
    FindAndPrintMateRet findAndPrintMateRet;
    findAndPrintMateRet.foundAndSearched = Pair<unsigned>(0,0);
    unsigned searched = 0, found = 0;
    map<Pair<string, unsigned> , vector<string> >::iterator itEnd = chrAndPosToNames.end();
    for (map<Pair<string, unsigned> , vector<string> >::iterator it=chrAndPosToNames.begin(); it != itEnd; ++it)
    {
        Pair<string, unsigned> chrAndPos = it->first;
        vector<string> readNames = it->second;
        std::sort(readNames.begin(), readNames.end());
        findAndPrintMateRet = findAndPrintMate(chrAndPos, readNames, bamPathIn, baiIndex);
        for (unsigned i=0; i<length(findAndPrintMateRet.readsToPrint); ++i)
            writeRecord(extendedBam, findAndPrintMateRet.readsToPrint[i]);
        clear(findAndPrintMateRet.readsToPrint);
        found += findAndPrintMateRet.foundAndSearched.i1;
        searched += findAndPrintMateRet.foundAndSearched.i2;
    }
    //cout << "Found mates for " << found << " reads out of " << searched << " reads." << endl; 
    return 0;
}
