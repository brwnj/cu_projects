from toolshed import nopen

def fastqparser(fastq):
    """yields name, seq, qual from fastq file"""
    fastq = nopen(fastq)
    line_num = -1
    record = []
    for line in fastq:
        line_num += 1
        if line_num == 4:
            yield record[0][1:], record[1], record[3]
            line_num = 0
            record = []
        record.append(line.strip())
    
    if record:
        if record[0]:
            yield record[0][1:], record[1], record[3]


class FastqReader(object):
    """Untested fastq class. Returns name, seq, qual."""
    def __init__(self, fastq):
        super(FastqReader, self).__init__()

    def __iter__(self):
        fq = nopen(self.fastq)
        id1  = fq.next().strip()
        seq  = fq.next().strip()
        id2  = fq.next().strip()
        qual = fq.next().strip()
        if qual == "":
            if id1 != "":
                sys.stderr.write(">> Incomplete fastq... skipping.\n")
            break

        yield id1[1:], seq, qual


def fastaparser(fasta):
    """yields name, seq or name, qual"""
    fasta = nopen(fasta)
    line_num = -1
    record = []
    for line in fasta:
        line_num += 1
        if line_num == 2:
            yield record[0][1:], record[1]
            line_num = 0
            record = []
        record.append(line.strip())

    if record:
        if record[0]:
            yield record[0][1:], record[1]