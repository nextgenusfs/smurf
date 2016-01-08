#!/usr/bin/env python

#script to mimic JCVI htab.pl, taking
import sys, os

#can construct dictionary of cutoffs - since this is apparently "secret" as it is hidden, will have to do trial and error to find right numbers
NoiseCutoff = {}
TrustCutoff = {}

def parseHmmer(input, output):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        from Bio import SearchIO
    with open(output, 'w') as output:
        with open(input, 'rU') as results:
            for qresult in SearchIO.parse(results, 'hmmer3-text'):
                hits = qresult.hits
                num_hits = len(hits)
                if len(hits) > 0:
                    try:
                        col1 = qresult.accession
                    except AttributeError:
                        col1 = qresult.id[0]
                    col2 = "01/01/2000"
                    col3 = qresult.seq_len
                    col4 = qresult.program
                    col5 = qresult.target
                    for i in range(0,num_hits):
                        col6 = hits[i].id
                        num_hsps = len(hits[i].hsps)
                        for x in range(0,num_hsps):
                            col7 = hits[i].hsps[x].query_start
                            col8 = hits[i].hsps[x].query_end
                            col9 = hits[i].hsps[x].hit_start
                            col10 = hits[i].hsps[x].hit_end
                            col11 = ""
                            col12 = hits[i].hsps[x].bitscore
                            col13 = hits[i].bitscore
                            col14 = hits[i].hsps[x].domain_index
                            col15 = hits[i].domain_obs_num
                            if qresult.description == '<unknown description>':
                                col16 = qresult.id
                            else:
                                col16 = qresult.description
                            col17 = hits[i].description
                            col18 = 1e-10
                            col19 = 1e-5
                            col20 = hits[i].evalue
                            col21 = hits[i].hsps[x].evalue
                            output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21))


#use sys/os to find list of files ending with out from directory, then parse each and finally write to file?
input = []
for file in os.listdir(sys.argv[1]):
    if file.endswith(sys.argv[2]):
        file = os.path.join(sys.argv[1], file)
        input.append(file)

for f in input:
    Output = f + '.htab'
    print "working on %s" % f
    parseHmmer(f, Output)

'''
 Description of the output format (tab-delimited, one line per domain hit)
    col  perl-col   description
    1      [0]      HMM accession
    2      [1]      Date of the htab job
    3      [2]      Length of the HMM (from database)
    4      [3]      Search method (typically hmmsearch or hmmpfam)
    5      [4]      Name of the database file
    6      [5]      Protein accession
    7      [6]      Start position of HMM match - hmm-f
    8      [7]      End position of HMM match - hmm-t
    9      [8]      Start position of the protein - seq-f
    10     [9]      End position of the protein  - seq-t
    11     [10]     (unassigned)
    12     [11]     Domain score
    13     [12]     Total score
    14     [13]     Domain number
    15     [14]     Total number of domains
    16     [15]     Biological description of the HMM (maybe truncated by hmmsearch or hmmpfam, but
                               can be recovered from DB using hmm_com_name if -s is not used)
    17     [16]     Biological description of the protein (maybe truncated by hmmsearch or hmmpfam)
    18     [17]     Trusted cutoff of the HMM (from database)
    19     [18]     Noise cutoff of the HMM (from database)
    20     [19]     Expect value for the whole match
    21     [20]     Expect value for the domain match
EXAMPLE
   ls  xxx/xxxx/xxxx/file_name* | htab.pl -f
   cat xxx/xxxx/xxxx/file_name* | htab.pl
'''
