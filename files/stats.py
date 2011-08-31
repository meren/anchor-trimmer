# -*- coding: utf-8 -*-

# a very simple script to get basic stats of a fasta file.
# helps to see if everything is all right with the sequence
# lenght and length variaiton. an example output looks like this:
#
# meren SSH://MBL $ python stats.py test-hmp-v3v5.fasta
# number of sequences    : 25000
# average sequence length: 448.03
# standard deviation     : 10.99
# minimum sequence length: 431
# maximum sequence length: 475
# meren SSH://MBL $
# 
# it also creates a very basic length distribution image and tries to save
# the figure in "sys.argv[1] + '.png'".
# 

import sys
import numpy

sys.path.append('..')
import fasta as u

fas = u.SequenceSource(sys.argv[1])

length_distro = [0] * 1000
seq_lengths = []
counter = 0

while fas.next():
    counter += 1
    seq_lengths.append(len(fas.seq))
    try:
        length_distro[len(fas.seq)] += 1
    except:
        pass

print 'number of sequences    : %d' % counter
print 'average sequence length: %.2f' % numpy.mean(seq_lengths or [0]) 
print 'standard deviation     : %.2f' % numpy.std(seq_lengths or [0]) 
print 'minimum sequence length: %d' % numpy.min(seq_lengths or [0]) 
print 'maximum sequence length: %d' % numpy.max(seq_lengths or [0])

try:
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize = (25, 10))
    
    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True)
    
    plt.plot(length_distro)
    plt.title(sys.argv[1])
    plt.xticks(range(0, 1000, 50))
    plt.xlabel('Sequence Length')
    plt.ylabel('Abundance')
    plt.savefig(sys.argv[1] + '.png')
except:
    pass
    # 'please go away' solution. this is how I usually deal
    # with my problems :/
