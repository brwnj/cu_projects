from toolshed import nopen

class FastqReader(object):
    """Yields name, seq, qual."""
    def __init__(self, fastq):
        self.fastq = nopen(fastq)

    def __iter__(self):
        fq = self.fastq
        while True:
            id1 = fq.next().strip()
            seq = fq.next().strip()
            id2 = fq.next().strip()
            qual = fq.next().strip()
            if qual == "":
                if id1 != "":
                    sys.stderr.write(">> Incomplete record... skipping.\n")
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