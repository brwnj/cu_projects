#! /usr/bin/env python

'''seqtools.formats: tools for reading file formats '''

import csv
import pysam
from collections import OrderedDict, namedtuple
from ._utils import line_to_fields

__version__ = '$Revision: 357 $'

class Reader(object):
    ''' Wrapper for the fields of a file specification '''
    def __init__(self, file, trackline=None, native=True, *args, **kwargs):
        self.file = file
        self.trackline = trackline
        self.native = native

    def __iter__(self):
        reader = csv.reader(self.file,
                            delimiter=self.delimiter)
        for line in reader:
            if len(line) == 0: continue

            if line[0].startswith('track') and not self.trackline: continue

            if self.native:
                yield DatumNative(line, self.field_names, self.field_types)
            else:
                yield Datum(line, self.field_names)

class DatumNative(object):

    def __init__(self, line, names, format, validate=None, *args, **kwds):

        self._data = []

        # map the formats to each line
        for idx, field in enumerate(line):
            name = names[idx]

            # if the value is empty (default is '.'), it's automatically a
            # string
            if field == '.':
                func = str
            else:
                func = format[name]

            try:
                self.__setattr__(name, func(field))
                self._data.append(func(field))
            except ValueError:
                raise TypeError, "Error converting %s to %s"%(field, name)

        return None 

    def __str__(self):
        return '\t'.join([str(i) for i in self._data])

class read_bed(Reader):

    ''' Parses the BED file format.  The first three fields (chrom,
    chromStart, chromEnd) are required.  The remaining nine are optional.
    Set the delimiter (default is tab) of the BED file manually by passing a delimiter
    value.

    Names and Format specs are:

        1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 

        The 9 additional optional BED fields are:

        4. name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        5. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
        6. strand - Defines the strand - either '+' or '-'.
        7. thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
        8. thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
        10. blockCount - The number of blocks (exons) in the BED line.
        11. blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        12. blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
    '''

    delimiter = '\t'

    field_types = OrderedDict([('chrom',str),('chromStart',int),
                               ('chromEnd',int),('name',str),
                               ('score',float), ('strand',str),
                               ('thickStart',str), ('thickEnd',str),
                               ('itemRgb',str), ('blockCount',int),
                               ('blockSizes',str),('blockStarts',str)])

    field_names = field_types.keys()

    def __init__(self, file, *args, **kwds):
        Reader.__init__(self, file=file)
    
        # option to set the delimiter manually
        if 'delimiter' in kwds:
            self.delimiter = kwds['delimiter']

class read_bedgraph(Reader):

    ''' Parses BEDGRAPH file format. '''

    delimiter = '\t'

    field_types = OrderedDict([('chrom',str),('chromStart',int),
                               ('chromEnd',int),('score',int)])
   
    field_names = field_types.keys()

    def __init__(self, filename, *args, **kwds):
        Reader.__init__(self, filename)
    
        # option to set the delimiter manually
        if 'delimiter' in kwds:
            self._delimiter = kwds['delimiter']

class read_mapview(Reader):
    ''' Reader for Maqview output. From the Maq documentation
        (http://maq.sourceforge.net/maq-manpage.shtml):
    
        Display the read alignment in plain text. For reads aligned before
        the Smith-Waterman alignment, each line consists of read name,
        chromosome, position, strand, insert size from the outer
        coorniates of a pair, paired flag, mapping quality, single-end
        mapping quality, alternative mapping quality, number of mismatches
        of the best hit, sum of qualities of mismatched bases of the best
        hit, number of 0-mismatch hits of the first 24bp, number of
        1-mismatch hits of the first 24bp on the reference, length of the
        read, read sequence and its quality. '''
        
    delimiter = "\t"

    field_types = OrderedDict([('readname',str),('chrom',str),
                               ('pos',int),('strand',str),('insert_size',int),
                               ('paired_flag',int),('map_qual',int),('se_map_qual',int),
                               ('alt_map_qual',int),
                               ('num_mismatch',int),('sum_qual',int),('num_zero_match',int),
                               ('num_one_match',int),('read_len',int),('read_seq',str),('read_qual',str)])

    field_names = field_types.keys()

    def __init__(self, filename, *args, **kwds):
        Reader.__init__(self, filename)
        # option to set the delimiter manually
        if 'delimiter' in kwds:
            self._delimiter = kwds['delimiter']

class read_gff(Reader):
    '''
    Names and formats are:
       1. seqname - The name of the sequence. Must be a chromosome or scaffold.
       2. source - The program that generated this feature.
       3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
       4. start - The starting position of the feature in the sequence. The first base is numbered 1.
       5. end - The ending position of the feature (inclusive).
       6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
       7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
       8. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
       9. group - All lines with the same group are linked together into a single item. 
    '''

    delimiter = '\t'

    field_types = OrderedDict([('seqname',str),('source',str),
                               ('feature',str),('start',int),
                               ('end',int),('score',float),
                               ('strand',str),('frame',int),
                               ('group',str)])

    field_names = field_types.keys()

    def __init__(self, filename, *args, **kwds):
        Reader.__init__(self, filename)
    
        # option to set the delimiter manually
        if 'delimiter' in kwds:
            self._delimiter = kwds['delimiter']

class read_psl(Reader):
    '''
    Names and formats are:
       1. matches - Number of bases that match that aren't repeats
       2. misMatches - Number of bases that don't match
       3. repMatches - Number of bases that match but are part of repeats
       4. nCount - Number of 'N' bases
       5. qNumInsert - Number of inserts in query
       6. qBaseInsert - Number of bases inserted in query
       7. tNumInsert - Number of inserts in target
       8. tBaseInsert - Number of bases inserted in target
       9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
      10. qName - Query sequence name
      11. qSize - Query sequence size
      12. qStart - Alignment start position in query
      13. qEnd - Alignment end position in query
      14. tName - Target sequence name
      15. tSize - Target sequence size
      16. tStart - Alignment start position in target
      17. tEnd - Alignment end position in target
      18. blockCount - Number of blocks in the alignment (a block contains no gaps)
      19. blockSizes - Comma-separated list of sizes of each block
      20. qStarts - Comma-separated list of starting positions of each block in query
      21. tStarts - Comma-separated list of starting positions of each block in target
    '''

    delimiter = '\t'

    field_types = OrderedDict([('matches',int),('misMatches',int),
                               ('repMatches',int),('nCount',int),
                               ('qNumInsert',int),('qBaseInsert',int),
                               ('tNumInsert',int),('tBaseInsert',int),
                               ('strand',str),('qName',str),('qSize',int),
                               ('qStart',int),('qEnd',int),('tName',str),
                               ('tSize',int),('tStart',int),('tEnd',int),
                               ('blockCount',int),('blockSizes',str),
                               ('qStarts',str),('tStarts',str)])

    field_names = field_types.keys()

    def __init__(self, filename, *args, **kwds):
        Reader.__init__(self, filename)
    
        # option to set the delimiter manually
        if 'delimiter' in kwds:
            self._delimiter = kwds['delimiter']

class read_bowtie(Reader):
    ''' 
    From http://bowtie-bio.sourceforge.net/manual.shtml#algn_reporting:

     The bowtie aligner outputs each alignment on a separate line. Each line is a collection of 8 fields separated by tabs; from left to right, the fields are:

       1. Name of read that aligned (readname)
       2. Orientation of read in the alignment, - for reverse complement,
       + otherwise (strand)
       3. Name of reference sequence where alignment occurs, or ordinal ID
       if no name was provided (chrom)
       4. 0-based offset into the forward reference strand where leftmost
       character of the alignment occurs (pos)
       5. Read sequence (reverse-complemented if orientation is -)
       (read_seq)
       6. ASCII-encoded read qualities (reversed if orientation is -). The
       encoded quality values are on the Phred scale and the encoding is
       ASCII-offset by 33 (ASCII char !). (read_qual)
       7. Number of other instances where the same read aligns against the
       same reference characters as were aligned against in this
       alignment. This is not the number of other places the read aligns
       with the same number of mismatches. The number in this column is
       generally not a good proxy for that number (e.g., the number in
       this column may be '0' while the number of other alignments with
       the same number of mismatches might be large). This column was
       previously described as "Reserved". (num_other_match)
       8. Comma-separated list of mismatch descriptors. If there are no
       mismatches in the alignment, this field is empty. A single
       descriptor has the format offset:reference-base>read-base. The
       offset is expressed as a 0-based offset from the high-quality (5')
       end of the read. (match_descrip)

    '''    
    delimiter = "\t"

    field_types = OrderedDict([('readname',str),('strand',str),
                               ('chrom',str),('pos',int),
                               ('read_seq',str),('read_qual',str),
                               ('num_other_match',int),
                               ('match_descrip',str)])

    field_names = field_types.keys()

    def __init__(self, filename, *args, **kwds):
        Reader.__init__(self, filename)

class read_scarf(Reader):

    ''' Reader for SCARF output from Illumina GA '''

    delimiter = ":"

    field_types = OrderedDict([('cell_id',str),('lane',int),
                               ('block',int),('region',int),
                               ('cluster',str),('sequence',str),
                               ('quals',str)])

    field_names = field_types.keys()

class read_crossmatch(object):
    ''' Parser for Crossmatch output.

    N.B. Must be used with the option "-tags" to cross_match to report
    parsable alignment and discrepancy information

    Add -discrep_lists to get discrepancies in the alignment.

    Available fields are:

    Alignment: 'sw_score','perc_sub','perc_del','perc_ins',
               'query_id','query_start','query_end','num_base_extend',
               'target_id','num_base_prior','target_start','target_end'

    Discrepancy: 'type','query_pos','nuc','target_pos', 'seq'
    '''
    def __init__(self, filehandle):
        self.filehandle = filehandle

    def __iter__(self):

        cur_alignment = None
        self.total_seqs = 0

        for line in self.filehandle:

            if line.startswith('Query file'):
                # get the query file
                fields = line_to_fields(line)
                self.query_file = fields[2]
                continue

            if line.startswith('Subject file'):
                # get the subject file
                fields = line_to_fields(line)
                self.subject_file = fields[2]
                continue

            if line.startswith('Sequence file'):
                # get the number of seqs
                fields = line_to_fields(line)
                self.total_seqs += int(fields[-2])
                continue

            # Alignment lines start with 2 spaces
            if line.startswith('ALIGNMENT'):

                # alignment line
                if cur_alignment:
                    yield cur_alignment
                    cur_alignment = None

                # removes the 'ALIGNMENT' field
                fields = line_to_fields(line)[1:]

                # find whether seq is complemented
                comp = False
                if 'C' in fields:
                    fields.pop(8)
                    comp = True
                    # flip the last 3 fields
                    rev_fields = fields[-3:]
                    rev_fields.reverse()
                    fields[-3:] = rev_fields

                cur_alignment = CrossmatchAlignmentInfo(fields)
                cur_alignment.complement = comp

            # Discrepancy types are D (deletion), S (transversion) and I
            # (transition)
            elif line.startswith('DISCREPANCY'):

                # removes the 'DISCREPANCY' field
                fields = line_to_fields(line)[1:]
              
                # update formatting
                type = fields[0]
                if '-' in type:
                    type, num_types = type.split('-')
                    num_types = int(num_types)
                else:
                    num_types = 1
              
                # XXX: check for ob1 error here
                query_pos = int(fields[1]) 
                nuc = fields[2]
                target_pos = int(fields[3])
                seq = fields[-1]

                # add multiples in case there are multiple deletions
                for i in range(num_types):

                    query_pos += i
                    target_pos += i

                    discrep_fields = [type, query_pos, nuc, target_pos,
                                      seq]

                    cur_alignment.add_discrep(discrep_fields)

discrepancy_fields = ('type','query_pos','nuc','target_pos','seq')
CrossmatchDiscrepancy = namedtuple('CrossmatchDiscrepancy',discrepancy_fields)

crossmatch_alignment_fields = ('sw_score','perc_sub','perc_del','perc_ins',
                               'query_id','query_start','query_end','num_base_extend',
                               'target_id','target_start','target_end','num_base_prior')
CrossmatchAlignment = namedtuple('CrossmatchAlignment',crossmatch_alignment_fields)

class CrossmatchAlignmentInfo(object):

    def __init__(self, fields):
        self.alignment = CrossmatchAlignment(*fields) 
        self.discreps = []
        self.complement = False

    def add_discrep(self,fields):
        self.discreps.append(CrossmatchDiscrepancy(*fields))

    def __str__(self):
        fields = (self.alignment.query_id, len(self.discreps),
                  self.complement)
        return 'CrossmatchAlingmentInfo: id=%s num_discreps=%d comp=%s' % fields

class read_miranda(object):
    ''' Parser for miranda output '''
    def __init__(self, filehandle):
        self.filehandle = filehandle

    def __iter__(self):

        cur_alignment = None

        for line in self.filehandle:

            if line.startswith('Query Filename'):
                # get the query file
                fields = line_to_fields(line,sep='\t')
                self.query_file = fields[1]
                continue

            if line.startswith('Reference Filename'):
                # get the subject file
                fields = line_to_fields(line,'\t')
                self.reference_file = fields[1]
                continue

            # Alignment lines start with 2 spaces
            if line.startswith('>') and not line.startswith('>>'):

                # alignment line
                if cur_alignment:
                    yield cur_alignment
                    cur_alignment = None

                # remove the leading '>' char
                line = line[1:]

                fields = line_to_fields(line,sep='\t')

                query_start, query_end = fields[6].split(' ')
                ref_start, ref_end = fields[7].split(' ')

                info = []
                info.extend(fields[:6])
                info.extend([query_start,query_end,
                             ref_start,ref_end])
                info.extend(fields[8:])

                cur_alignment = MirandaAlignmentInfo(info)

miranda_alignment_fields = ('query_id','ref_id','score','free_energy','z_score',
                            'p_value','query_start','query_end','ref_start','ref_end',
                            'aln_length','perc_id','perc_sim')
MirandaAlignment = namedtuple('MirandaAlignment',miranda_alignment_fields)

class MirandaAlignmentInfo(object):

    def __init__(self, fields):
        self.alignment = MirandaAlignment(*fields) 

    def __str__(self):
        fields = (self.alignment.query_id, self.alignment.ref_id)
        return 'MirandaAlingmentInfo: query_id=%s ref_id=%s' % fields