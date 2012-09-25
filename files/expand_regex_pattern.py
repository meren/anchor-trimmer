# sorry about the terrible code. example use:
#
#     python expand_regex_pattern.py A.T[C,T]AAA.[A,G]AAT[A,T]GACGG
#

import sys

rp = sys.argv[1]
expant = ['']
valid_bases = ['A', 'T', 'C', 'G']

pos = 0
while 1:
    if rp[pos] in valid_bases:
        for i in range(0, len(expant)):
            expant[i] += rp[pos]
        pos += 1
    
    elif rp[pos] == '.':
        new_stuff = []
        for i in range(0, len(expant)):
            for base in valid_bases:
                new_stuff.append(expant[i] + base)
        expant = new_stuff
        pos += 1

    elif rp[pos] == '[':
        start = pos
        end = pos
        while rp[end] != ']':
            end += 1
        bases = rp[start+1:end].split(',')
        new_stuff = []
        for i in range(0, len(expant)):
            for base in bases:
                new_stuff.append(expant[i] + base)
        expant = new_stuff

        pos = end + 1

    else:
        for i in range(0, len(expant)):
            expant[i] += rp[pos]
        pos += 1

    if pos == len(rp):
        break
    else:
        continue

for e in expant:
    print e
