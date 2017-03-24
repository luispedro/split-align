import subprocess
import os
class BWA(object):
    '''BWA Caller'''
    def __init__(self):
        pass


    def index(self, fname):
        subprocess.call(['bwa', 'index', fname])

    def align(self, fname, fqnames):
        return subprocess.Popen(['bwa', 'mem', fname] + fqnames, stdout=subprocess.PIPE).stdout

class FASequence(object):
    '''FASTA Sequence'''
    def __init__(self, h, content):
        self.header = h
        self.content = content

    def nbases(self):
        return len(self.content)

    def write_to(self, output):
        output.write(self.header)
        output.write('\n')
        output.write(self.content)
        output.write('\n')

def read_fasta(fname):
    '''Iterator of FASequence over `fname`'''
    if fname.endswith('.gz'):
        import gzip
        ifile = gzip.open(fname, 'rt')
    else:
        ifile = open(fname, 'rt')
    cur_h = None
    content = []
    for line in ifile:
        if line[0] == '>':
            if cur_h is not None:
                yield FASequence(cur_h, ''.join(content))
            cur_h = line.strip()
            content = []
        else:
            content.append(line.strip())


def split(fname, max_size):
    '''Split a FASTA file into blocks of `max_size`
    '''
    padding = 150
    splits = []
    ix = 0
    output = None
    current = max_size + 1
    for seq in read_fasta(fname):
        if seq.nbases() + current + padding > max_size:
            if output is not None:
                output.close()
            output = '{}.split.{}'.format(fname, ix)
            splits.append(output)
            output = open(output, 'wb')
            ix += 1
            current = 0
        seq.write_to(output)
        current += seq.nbases() + padding
    print("Finished splits ({} splits in total)".format(ix))
    return splits

def index(base, fname, max_size):
    '''Split and index `fname`
    '''
    splits = split(fname, max_size)
    for sp in splits:
        base.index(sp)

def find_splits(fname):
    from glob import glob
    import re
    names = [f for f in glob('{}.split.*'.format(fname))
            if re.match(r'{}\.split\.\d+$'.format(fname), f)]
    try:
        import natsort
        names.sort(key=natsort.natsort_keygen())
    except ImportError:
        names.sort()
    return names

def strip_seq_quals(line):
    tokens = line.split(b'\t')
    tokens[9] = b'*'
    tokens[10] = b'*'
    return b'\t'.join(tokens)


def align_to(base, fname, fqfiles, ofname):
    import tempfile
    temporaries = []
    splits = find_splits(fname)
    for i, sp in enumerate(splits):
        curtmp = tempfile.NamedTemporaryFile(delete=False)
        for line in base.align(sp, fqfiles):
            if i > 0 and line[0] != b'@'[0]:
                line = strip_seq_quals(line)
            curtmp.write(line)
        print("Aligned {}/{}".format(i + 1, len(splits)))
        temporaries.append(curtmp.name)
        curtmp.close()
    print("Merging outputs...")
    merge_temporaries(temporaries, ofname)
    for t in temporaries:
        os.unlink(t)

def merge_temporaries(temporaries, ofname):
    print(temporaries, ofname)
    inputs = [read_samgroup(open(t, 'rt')) for t in temporaries]
    with open(ofname, 'wt') as output:
        for i,gs in enumerate(zip(*inputs)):
            if i == 0:
                # headers
                for g in gs:
                    for line in g:
                        output.write(line)
            else:
                gs = merge1(gs)
                for g in gs:
                    output.write(g)

def read_samgroup(ifile):
    curg = []
    did_header = False
    for line in ifile:
        if not did_header:
            if line[0] == '@':
                curg.append(line)
            else:
                yield curg
                did_header = True
                curg = [line]
                actid = line.split('\t', 2)[0]
            continue
        rid = line.split('\t', 2)[0]
        if rid != actid:
            yield curg
            curg = [line]
            actid = rid
        else:
            curg.append(line)
    if curg:
        yield curg

def is_mate1(g):
    tokens = g.split('\t')
    flags = int(tokens[1])
    return (flags & 0x40)

def merge1(gs):
    gs = [g for gss in gs for g in gss]
    ids = ([g.split('\t')[0] for g in gs])
    if len(set(ids)) > 1:
        raise ValueError('Mismatch: {}'.format(ids))
    mate1 = allbest([g for g in gs if is_mate1(g)])
    mate2 = allbest([g for g in gs if not is_mate1(g)])
    return mate1 + mate2

def allbest(gs):
    if not gs: return gs
    def score(g):
        tokens = g.split('\t')
        if tokens[4] == '*':
            return -1
        return int(tokens[4])
    scores = [score(g) for g in gs]
    max_score = max(scores)
    return [g for g,s in zip(gs, scores) if s == max_score]

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Split short read align')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    parser_index = subparsers.add_parser('index')
    parser_index.add_argument('database', metavar='DB', type=str,
                        help='Database to index')
    parser_index.add_argument('--max-size', metavar='MEGABASEPAIRS', type=int,
                        dest='size',default=2000,
                        help='Maximum number of mega-basepairs')

    parser_align = subparsers.add_parser('align')
    parser_align.add_argument('database', metavar='DB', type=str,
                        help='Database to align against')
    parser_align.add_argument('input1', metavar='FQ1', type=str, default=None,
                        help='Input file (first mate)')
    parser_align.add_argument('input2', metavar='FQ2', type=str, default=None,
                        help='Input file (second mate, optional)')
    parser_align.add_argument('extra', metavar='EXTRA-ARGS', type=str, default='',
                        help='Extra arguments (passed verbatim to bwa)')
    args = parser.parse_args()
    base = BWA()
    if args.subcommand == 'index':
        index(base, args.database, args.size * 1000 * 1000)
    elif args.subcommand == 'align':
        align_to(base, args.database, [args.input1, args.input2] + args.extra.split(), 'output.sam')
    else:
        raise ValueError("Unknown subcommand '{}'".format(args.subcommand))

if __name__ == '__main__':
    main()
