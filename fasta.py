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
                     sorry for this inconvenience. obviously someone hired a very lazy programmer :/
                     
                     and now this primitive library will go beyond this and call a sys.exit on your client.
                     
                     please remove extra new line characters from your sequences. to do that, you can use this
                     pesky script:

                     ---->8---->8---->8---->8---->8---->8---->8---->8---->8---->8---->8---->8-----
                     #!/usr/bin/python
                     # -*- coding: utf-8 -*-
                     
                     import sys
                     
                     fasta_file = sys.argv[1]
                     
                     fasta_lines = open(fasta_file).readlines()
                     fasta_lines.append("#")
                     
                     new_fasta = open(fasta_file + "-NEW-LINES-REMOVED", "w")
                     
                     for i in range(0, len(fasta_lines) - 1):
                         if fasta_lines[i].startswith('#') or fasta_lines[i].startswith('>'):
                             new_fasta.write(fasta_lines[i])
                         elif fasta_lines[i + 1].startswith('#') or fasta_lines[i + 1].startswith('>'):
                             new_fasta.write(fasta_lines[i])
                         elif fasta_lines[i][-1] == '\\n':
                             new_fasta.write(fasta_lines[i].strip())
                     
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
