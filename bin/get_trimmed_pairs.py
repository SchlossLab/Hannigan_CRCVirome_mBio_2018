#!/usr/bin/env python

import sys
import optparse


def get_intersection(fwd, rev):
    s1, s1_dict = set_dict(fwd)
    s2, s2_dict = set_dict(rev)
    s1_s2 = s1.intersection(s2)
    return s1_s2, s1_dict, s2_dict

def set_dict(inf):
    s = set([])
    identifier = ''
    s_dict = {}
    seq = ''
    parser =  ParseFastQ(inf)
    for record in parser:
        #print record[0]
        all = record[0].split(' ')
        identifier = all[0]
        s.add(identifier)
        s_dict[identifier] = (all[1]+'\n'+record[1]+'\n'+record[2]+'\n'+record[3])
        
    return s,s_dict
    
def make_new_fastq(inters, dict1, dict2, o1, o2):
    fd1 = open(o1, 'w')
    fd2 = open(o2, 'w')
    for header in inters:
        seq_for = dict1[header]
        seq_rev = dict2[header]
        fd1.write('%s %s\n' % (header, seq_for))
        fd2.write('%s %s\n' % (header, seq_rev))
    
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self

    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else:
                elemList.append(None)

        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line #%s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file and try again **" % (self._hdSyms[0])
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file and try again **" % (self._hdSyms[1])

        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)

    def getNextReadSeq(self):
        """Convenience method: calls self.getNext and returns only the readSeq."""
        try:
            record = self.next()
            return record[1]
        except StopIteration:
            return None


if __name__ == '__main__':
    '''
    This program will get an intersection of two fastq files. It will take two fasta files
    as input and output the intersection of both.
    '''

    usage = '%prog {options}'
    p = optparse.OptionParser(usage=usage)
    p.add_option('-f', '--input_one', default=None,
        help='First input fastq file.')
    p.add_option('-s', '--input_two', default=None,
        help='Second input fastq file.')
    p.add_option('-o', '--output_one', default=None,
        help='First output fastq file.')
    p.add_option('-t', '--output_two', default=None,
        help='Secound output fastq file.')
    opts, args = p.parse_args()
    
    if opts.input_one is None:
        p.error('Must provide input filepath')
        
    if opts.input_two is None:
        p.error('Must provide input filepath')

    s1_s2, s1_dict, s2_dict = get_intersection(opts.input_one, opts.input_two)
    make_new_fastq(s1_s2, s1_dict, s2_dict, opts.output_one, opts.output_two)