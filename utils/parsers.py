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


class ReadFastq(object):
    """Untested fastq class. Yields name, seq, qual."""
    def __init__(self, fq):
        super(FastqReader, self).__init__()

    def __iter__(self):
        id1  = fq.next().strip()
        seq  = fq.next().strip()
        id2  = fq.next().strip()
        qual = fq.next().strip()
        if qual == "":
            if id1 != "":
                sys.stderr.write(">> Incomplete fastq... skipping.\n")
            break
        yield id1[1:], seq, qual


def read_fasta(fa):
    r"""yields name and seq from fasta.
    
    with open('f.fasta') as fp:
        for name, seq in read_fasta(fp):
            print(name, seq)
    """
    name, seq = None, []
    for line in fa:
        line = line.rstrip()
        if line.startswith('>'):
            if name: yield (name, ''.join(seq))
            name, seq = line.lstrip('>'), []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))