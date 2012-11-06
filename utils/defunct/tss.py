import pymysql
import os
import pybedtools
import matplotlib.pylab as pylab
import numpy

outputbed = 'tss.bed'
fairedir = 'FAIRESEQ_seq1'
tsswindowsize = 2500
bins = 50
# make a bed file consisting of chromosome and tss start and tss start + 1
if not os.path.exists(outputbed):
    obed = file(outputbed,'w')
    conn = pymysql.connect(host='genome-mysql.cse.ucsc.edu', port=3306, user='genome', db='hg19')
    cur = conn.cursor()
    cur.execute('select distinct chrom, txStart, txStart+1 from knownGene order by chrom,txStart')
    for line in cur.fetchall():
        chrom = line[0].strip("chr")
        start = str(line[1])
        end = str(line[2])
        obed.write("\t".join([chrom, start, end]) + "\n")

# overlap the tss bed plus 2500 bases upstream and downstream with the faireseq peaks
# only record distance between tss and fseq peak
a = pybedtools.BedTool(outputbed)
totaltss = a.count()
totalpeaks = 0
distances = []
for f in os.listdir(fairedir):
    b = pybedtools.BedTool(fairedir + "/" + f)
    totalpeaks = totalpeaks + b.count()
    for overlap in a.window(b, w=tsswindowsize):
        distances.append(float(overlap[1]) - float(overlap[4]))

print "1. Peaks within %ibp of TSS: %i" % (tsswindowsize,len(distances))
print "2. Total Peaks: %i" % totalpeaks
print "3. Total TSS: %i" % totaltss

# plot the distances from tss
binsIn = numpy.arange(-2500,2500,bins)
pylab.subplot(111)
pylab.hist(distances, binsIn)
pylab.show()