# -*- coding: utf-8 -*-
#
# this is a script that was used to generate valid anchor sequences file for v3v5
# region. can be modified and used for different regions.
#
# expected input, is the output generated from SILVA alignments. RADME file may give 
# more ideas about how to get to this point. Neverthless, the input actually should
# look like this:
#
# A14565_8_1539_223180    GTAGCGGTGAAAT   223180
# AABF01000111_2354_3854_33662    GTAGAGGTGAAAT   33662
# AAAL02000003_178628_180172_28190    GTAGCAGTGAAAT   28190
# AACY020108912_397_1908_26802    GTAGCGGTGGAAT   26802
# AACY020323954_4416_5934_5376    GTAGCGGTAAAAT   5376
# AANZ01000021_102750_104251_3717 GGAGCGGTGAAAT   3717
# AAFJ01000001_39328_40836_3714   GTAGGGGTAAAAT   3714
# AB022035_1_1452_2994    GTAGTGGTGAAAT   2994
# AB027017_1_1482_2953    GTAGGGGTGAAAT   2953
# AB015558_1_1608_2493    GTAGGAGTGAAAT   2493
# AB062809_1_1375_1433    GGAGCGGTGGAAT   1433
# AAQL01003796_110_1491_977   GTAACGGTGGAAT   977
# AB010425_1_1521_636 GTAGTGGTAAAAT   636
# AACY024057770_1_1461_494    GTAGCGGTGACAT   494
# AACY020179327_219_1708_453  GTAGTGGTGGAAT   453
# AB016974_1_1522_430 GTAAGGGTGGAAT   430
# AB034115_1_1506_381 GTAGTAGTGAAAT   381
# AB025965_1_1426_377 GTAGAGGTAAAAT   377
# AB006998_1_1445_352 GTAACGGTGAAAT   352
# AB179511_1_1496_343 GAAGCGGTGAAAT   343
# AB176200_1_1407_333 GTAGCGGTGAAAAT  333
# AB089089_1_1318_328 GTAGAGGTGGAAT   328
# AAYO02000002_77347_78873_252    GGAGCGGTAAAAT   252
# AB195927_1_1478_242 GTAGGGGTGGAAT   242
# AB088949_1_1330_242 GTAGCGGCGAAAT   242
# AB039006_1_1441_241 GTAGCGGTGAGAT   241
# (...) 
#
# 
# first column is the label, which is not important and could be anything. second column
# is the anchor sequence, a window that was cut from SILVA alignments and then gaps were
# removed. third column is the frequency of the given anchor among those sequences in
# SILVA.
#
# this app will create 'v3v5-valid_anchor_sequences.txt' for anchortrimmer.py to use.
#

import numpy as np
import sys

# minimum expected frequency of a nucleotide for a given position in the anchor (percent):
minimum_frequency = 1

silva_results = [line.strip() for line in open(sys.argv[1]).readlines()]

with_unexpected_base = []
with_unexpected_length = []
valid_but_singleton = []
valid_sequences = []

for line in silva_results:
    try:
        _, sequence, frequency = line.split()
        frequency = int(frequency)
    except:
        pass

    if len([c for c in sequence if c not in "ATCG"]):
        with_unexpected_base.append((sequence, frequency),)
        continue

    if len(sequence) != 13:
        with_unexpected_length.append((sequence, frequency),)
        continue

    if frequency == 1:
        valid_but_singleton.append((sequence, frequency),)
        continue

    valid_sequences.append((sequence, frequency),)

print 'sequences eliminated due to having unexpected bases: %d (mean frequency: %.2f)' % (len(with_unexpected_base), np.mean([x[1] for x in with_unexpected_base]))
print 'sequences eliminated for not having 13 nucleotides : %d (mean frequency: %.2f, min length: %d, max length: %d)' % (len(with_unexpected_length), np.mean([x[1] for x in with_unexpected_length]), np.min([len(x[0]) for x in with_unexpected_length]), np.max([len(x[0]) for x in with_unexpected_length]))
print 'singletons that passed first two (and eliminated)  : %d (mean frequency: %.2f)' % (len(valid_but_singleton), np.mean([x[1] for x in valid_but_singleton]))
print 'sequences survived for further analysis            : %d (mean frequency: %.2f)' % (len(valid_sequences), np.mean([x[1] for x in valid_sequences]))


# create consensus sequence from valid sequences
consensus_sequence = ''
for pos in range(0, len(valid_sequences[0][0])):
    bases = {}
    for sequence, frequency in valid_sequences:
        base = sequence[pos]
        if bases.has_key(base):
            bases[base] += frequency
        else:
            bases[base] = frequency
    
    frequent_bases = []
    for base in "ATCG":
        if bases.has_key(base):
            if bases[base] * 100.0 / sum(bases.values()) >= minimum_frequency:
                frequent_bases.append(base)
    
    if len(frequent_bases) == 1:
        consensus_sequence += frequent_bases[0]
    elif len(frequent_bases) == 4:
        consensus_sequence += '.'
    else:
        consensus_sequence += '[%s]' % ','.join(frequent_bases)

print '\n\tconsensus anchor sequence: ', consensus_sequence

# store sequences
print '\n(storing valid anchor sequences in the order of frequency)'
f = open('v3v5-valid_anchor_sequences.txt', 'w')
for sequence in [s[0] for s in sorted(valid_sequences, key=lambda x: x[1], reverse = True)]:
    f.write(sequence + '\n')
f.close()

try:
    ###############################
    # visualize base frequencies..
    ###############################
    
    import pylab
    import matplotlib.pyplot as plt
    
    print
    print '(generating base frequencies image)'
    base_pos = {'A': 4, 'T': 3, 'C': 2, 'G': 1}
    base_colors = {'A': 'red', 'T': 'green', 'C': 'blue', 'G': 'yellow'}
    
    total_reads = sum([x[1] for x in valid_sequences])
    
    N = len(valid_sequences[0][0])
    ind = np.arange(N) + 1
    
    fig = plt.figure(figsize = (N, 5))
    
    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True)
    
    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.025, top = 0.92, bottom = 0.05)
    
    ax = fig.add_subplot(111)
    ax.plot([0], [0], visible = False)
    ax.plot([N], [0], visible = False)
    ax.plot([0], [5], visible = False)
    
    for pos in range(0, len(valid_sequences[0][0])):
        bases = {}
        for sequence, frequency in valid_sequences:
            base = sequence[pos]
            if bases.has_key(base):
                bases[base] += frequency
            else:
                bases[base] = frequency
        for base in bases:
            ratio = bases[base] * 1.0 / total_reads
            ax.plot([pos + 1], [base_pos[base]], 'o', c = 'white', lw = 1, ls="-", alpha = 0.75, ms = ratio * 100)
            ax.plot([pos + 1], [base_pos[base]], 'o', c = base_colors[base], lw = 1, alpha = ratio / 5, ms = ratio * 100)
    
    
    for sequence, frequency in valid_sequences:
        ratio = frequency * 1.0 / total_reads
        ray = [base_pos[p] for p in sequence]
        ax.plot(ind, ray, c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 8, alpha = ratio / 2)
        ax.plot(ind, ray, c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 6, alpha = ratio / 2)
        ax.plot(ind, ray, c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 4, alpha = ratio / 2)
        ax.plot(ind, ray, c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 1, alpha = ratio / 2)
        ax.plot(ind, ray, c = 'white', solid_capstyle = "round", solid_joinstyle = "round", lw = 1, alpha = ratio * 2)
    
    pylab.yticks(np.arange(6), ('', 'G', 'C', 'T', 'A'), size = 'x-large')
    pylab.title("Base Frequencies")
    
    locs = range(0, N + 2)
    pylab.xticks(locs, [''] + ["P" + str(x) for x in range(0, len(locs))[1:-1]] + [''])
    
    try:
        plt.show()
    except:
        pass
    
    plt.savefig('v3v5-frequency-distribution.png')
    print '("v3v5-frequency-distribution.png" has been saved)'
except:
    pass
