Anchor trimming for amplicon sequences.
-----------------------------------------------------

This is a solution developed in MBL to originally trim long 454 sequences from somewhere between the forward and distal primers, when the amplicon length of PCR process is too much for 454 machine to sequence the entire molecule. By using a semi-conserved anchor somewhere in the reads, it was possible to trim sequences not based on the length (which could vary drastically from taxon to taxon) but based on a biologically viable place which could be prefferable for the consistensy of downstream analyses.

It turns out it also works quite beautifully to trim primers as well (as I recently started using it to trim EMP primers from MiSeq runs).

This application requires a couple of things. When it is run from the command line, besides the obvious things such as input file or output destination, it expects 2 things that require some explanation:

* (1) valid anchors for given region (denoted with parameter '-a'): This is a file that must be generated before hand. There may be different ways, but this is how we did it: look at the secondary structure of 16S rRNA molecule (http://www.rna.icmb.utexas.edu could be a nice place to find it) to decide which region between your primers seem to be conserved enough to use as an anchor (anchor needs to be as conserved as possible so it could be found in every sequence in your bacterial samples for the trimming purpose). Once the place of conserved oligonucleotide is known, universal SILVA alignments can be used to extract EVERY possible anchors from every bacterial species in SILVA database, gaps would be removed from those alignments, sequneces would be sorted by their frequency, singletons / ones that have ambigious bases / shorter or longer than expected anchor sequence length would be eliminated. What is left would be a list of sequences that are valid anchors, and the file would look like this:

             meren SSH://MBL $ head v3v5-valid_anchor_sequences.txt 
             GGATTAGATACCC
             GGATTAGAGACCC
             GGCTTAGATACCC
             AGATTAGATACCC
             GGATTAGAAACCC
             GAATTAGATACCC
             GGATTAGATACCT
             GAATTAGATACTC
             GAATTAGAGACTC
             GGATTAAATACCC
             meren SSH://MBL $ wc -l v3v5-valid_anchor_sequences.txt 
             226 v3v5-valid_anchor_sequences.txt

* (2) maximum divergence (denoted with parameter '-d'): Due to natural variation and/or sequencing error, some sequences may have slightly different sequences than the valid anchor sequences observed in the database. Levenshtein edit distance is being used here to compensate that and minimize the number of raw sequneces that cannot be trimmed because of that. 0.9 is the default value of maximum divergence. which means, if our anchor sequence template is 13 nucleotides long, sequence found in the anchor spot could be 1 nucleotide different than a valid anchor sequence (that difference could be an insertion, deletion or substitution). if 0.8 would be used as maximum divergence, than 2 nucleotide difference would be allowed for a sequence to be kept for investigation. This is not a greedy algorithm, search would stop ONLY if perfect match for a valid anchor sequence is found. Otherwise, every valid anchor sequences for every oligonucleotide in a given window is being tested for their divergences and the best one among all would be picked (then the list would be reorganized and program would optimize itself during the runtime). If you define an output file, a file with a '-FAILED' prefix to your output file will be automatically generated. You are encouraged to examine that file to see how many of your sequences are failed due to a failed anchor search, and tweak your max divergence parameter (and even maybe your valid anchor sequences file) according to that output, if necessary.

 
Class 'Settings' contains some predefined regions and values. New ones could be added, and when a correctly formatted valid anchor sequences file is provided, this application should work on any region in 16S rRNA gene.

If you don't provide an 'output' file name, program basically will print sequences with anchor matches higlighted, which could help to debug or test for different anchor sites or sequences.

If you have questions, please don't hesitate to send an e-mail.


Questions / remarks can be sent to A. Murat Eren, <meren / uchicago.edu>

