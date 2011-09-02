import sys

#
# poor man's fastalib.
# 

class SequenceSource:
    """
        This class is to read ids and sequences from FASTA formatted files
        in a structured manner. For now it expects sequnces to be only one
        line (in this sense it doesn't meet the FASTA standards). Every id
        sequence pair should look like this in the file:

        (...)
        >sequenceidforasequence
        ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
        >sequenceidforanothersequence
        ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
        (...)
    """
    def __init__(self, f_name):
        self.pos = 0
        self.id  = None
        self.seq = None
        self.file_content = [x.strip() for x in open(f_name).readlines()]
        self.total_seq = len(self.file_content) / 2
    def next(self):
        if (self.pos + 1) > self.total_seq:
            return False
        if self.file_content[self.pos * 2][0] != '>':
            print '''
                     this primitive fasta library assumes that sequencs are not split into multiple lines..
                     sorry for this inconvenience. obviously someone was very very lazy :/
                     
                     and now this primitive library will go beyond this and call a sys.exit on your client.
                     
                     please remove extra new line characters from your sequences. to do that, you can use this
                     pesky script:

                     ---->8---->8---->8---->8---->8---->8---->8---->8---->8---->8---->8---->8-----
                     #!/usr/bin/python
                     # -*- coding: utf-8 -*-
                     
                     import os
                     import sys
                     
                     fasta_file = sys.argv[1]
                     
                     if not os.path.exists(fasta_file):
                         print "No such file.."
                         sys.exit(0)
                     
                     fasta = open(fasta_file)
                     
                     print 'har har'
                     
                     new_fasta = open(fasta_file + ".new", "w")
                     sequence = []
                     while 1:
                         line = fasta.readline()
                     
                         if not line:
                             if len(sequence):
                                 new_fasta.write(''.join(sequence) + '\n')
                             break
                     
                         if line.startswith('>'):
                             if len(sequence):
                                 new_fasta.write(''.join(sequence) + '\n')
                             new_fasta.write(line)
                             sequence = []
                     
                         else:
                             sequence.append(line.strip())
                         
                     new_fasta.close()
                     ----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<-----

                     thanks.

            '''
            sys.exit(-1)
        self.id = self.file_content[self.pos * 2][1:].strip()
        self.seq = self.file_content[self.pos * 2 + 1].strip()
        self.pos += 1
        return True
    def load_id(self, id):
        while self.next():
            if self.id == id:
                break
    def reset(self):
        self.pos = 0
        self.id  = None
        self.seq = None
