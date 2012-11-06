from itertools import izip
import argparse
import os
import re
import sys


def nopen(f, mode="rb"):
    if not isinstance(f, basestring):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r": return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
         else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
         else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
         else urllib.urlopen(f) if f.startswith(("http://", "https://",
             "ftp://")) \
        else open(f, mode)


def reader(fname, header=True, sep="\t"):
    line_gen = (l.rstrip("\r\n").split(sep) for l in nopen(fname))
    if header == True:
        header = line_gen.next()
        header[0] = header[0].lstrip("#")
    if header:
        for toks in line_gen:
            yield dict(izip(header, toks))
    else:
        for toks in line_gen:
            yield toks


def make_dict(fname, key_name, sep="\t"):
    r"""Returns dictionary of file `fname` and assumes a header is present"""
    file_d = {}
    for d in reader(fname, sep=sep):
        file_d[d[key_name]] = d
    return file_d


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def rna_reference(fname):
    """From Ensembl ncRNA Fasta, returns dictionary of transcript
    biotype by Ensembl transcript name
    """
    transcript_bio = {}
    for f in reader(fname, header=False, sep=" "):
        if(f[0].startswith(">")):
            # dict[transcript_name] = {biotype, [chr, start, end, strand]}
            transcript_bio[f[0][1:]] = [f[-1].split(":")[-1],f[2].split(":")[2:]]
    return transcript_bio


def tracking_processing(paths):
    """Takes a dictionary of cuffdiff file paths and writes an annotated
    isoform tracking sheet with expression and significance testing results
    as well as a transcript's biotype.
    """
    # sort isoforms.fpkm_tracking
    if not os.path.exists("%s.sorted" % paths['isotracking']):
        os.system("head -n1 %s > %s.sorted" 
                    % (paths['isotracking'], paths['isotracking']))
        os.system("""awk 'BEGIN{OFS="\t"}\
                            (/^TCONS/){print $0}' %s \
                     | sort -f -k5,5 -k6,6 \
                     >> %s.sorted""" \
                    % (paths['isotracking'], paths['isotracking']))
    f = open('%s.anno' % paths['isotracking'], 'w')

    # Create dictionaries for each file to compare to isoforms.fpkm_tracking
    disot = make_dict(paths['isotracking'], 'tracking_id')
    dgene = make_dict(paths['geneexp'], 'test_id')
    diso = make_dict(paths['isoformexp'], 'test_id')
    dtss = make_dict(paths['tssexp'], 'test_id')
    dspl = make_dict(paths['splicing'], 'test_id')
    dpro = make_dict(paths['promoter'], 'test_id')
    # drna = rna_reference(paths['reference'])

    # Use reader to pull first line and create header
    mainheader = reader(paths['isotracking'], header=False).next()
    mainheader.insert(mainheader.index('tss_id'), 'tss_chart_id')
    expheader = reader(paths['isoformexp'], header=False).next()
    # remove trailing "significant"
    expheader.pop()
    # rename identical fields of 'expression' header
    for i, v in enumerate(expheader):
        expheader[i] = "%s_de" % v
    mainheader.extend(expheader)
    # mainheader.extend(['transcript_bt', 'bt_chromosome', 'bt_strand'])
    if "ln(fold_change)_de" in expheader:
        mainheader.extend(['gene_ln(fold_change)'])
    else:
        mainheader.extend(['gene_log2(fold_change)'])
    sigtesting = ['gene_sig', 'isoform_sig', 'tss_sig', 'splicing_sig', 'promoter_sig']
    mainheader.extend(sigtesting)
    f.write("\t".join(mainheader))
    f.write("\n")

    iter1 = reader(paths['isotracking'] + '.sorted')
    iter2 = reader(paths['isotracking'] + '.sorted')
    iter2.next()
    i = 1
    for line, nextline in izip(iter1, iter2):
        # retrieving dictionary entries for current transcript
        gene = dgene.get(line['gene_id'])
        isoform = diso.get(line['tracking_id'])
        tss = dtss.get(line['tss_id'])
        splicing = dspl.get(line['tss_id'])
        promoter = dpro.get(line['gene_id'])
        # try:
        #     reference = drna.get(line['nearest_ref_id'])
        # except:
        #     reference = False
        # append isoform expression table
        for v in expheader:
            line[v] = isoform[v[:-3]]
        # labeling tss with something shorter and easier identifiable
        line['tss_chart_id'] = str(i)
        # gene fold change
        try:
            line['gene_log2(fold_change)'] = gene['log2(fold_change)']
        except:
            line['gene_ln(fold_change)'] = gene['ln(fold_change)']
        # significance testing
        for item in sigtesting:
            d = eval(item.split("_")[0])
            line[item] = str((d and d['status'] == 'OK' and \
                            d['significant'] == 'yes') or 'false').lower()
        # handling tss bar id increments
        try:
            next_id = nextline['gene_short_name']
        except:
            next_id = False
        if(next_id and line['gene_short_name'] == next_id and \
                line['gene_short_name'] != '-'):
            i += 1
        else:
            i = 1
        # if reference: # lincRNA annotation library
        #     line['transcript_bt'] = reference[0]
        #     line['bt_chromosome'] = '%s:%s-%s' % (reference[1][0],\
        #                                             reference[1][1],\
        #                                             reference[1][2])
        #     line['bt_strand'] = reference[1][3]
        # else:
        #     line['transcript_bt'] = "-"
        #     line['bt_chromosome'] = "-"
        #     line['bt_strand'] = "-"
        f.write("\t".join(line[h] for h in mainheader))
        f.write("\n")
    f.close()
    return True


def cuffdiff(cuffdiffdir):
    cdiff = {}
    cdiff['isotracking'] = "%s/isoforms.fpkm_tracking" % cuffdiffdir
    cdiff['isoformexp'] = "%s/isoform_exp.diff" % cuffdiffdir
    cdiff['geneexp'] = "%s/gene_exp.diff" % cuffdiffdir
    cdiff['tssexp'] = "%s/tss_group_exp.diff" % cuffdiffdir
    cdiff['splicing'] = "%s/splicing.diff" % cuffdiffdir
    cdiff['promoter'] = "%s/promoters.diff" % cuffdiffdir
    if(tracking_processing(cdiff)):
        sys.stderr.write('isoforms.fpkm_tracking has been successfully annotated.\n')


def __main__():
    p = argparse.ArgumentParser(description=__doc__)
    required = p.add_argument_group("required arguments")
    required.add_argument("-cd", "--cuffdiff-dir", dest='cuffdiff_dir', help="Directory containing cuffdiff files")
    args = p.parse_args()
    if args.cuffdiff_dir and os.path.isdir(args.cuffdiff_dir):
        cuffdiff(os.path.abspath(args.cuffdiff_dir))
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        __main__()
