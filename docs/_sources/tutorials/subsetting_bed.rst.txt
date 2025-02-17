
.. _subsetting_bed:

Getting read pileup inside peaks (intervals)
---------------------------------------------

This is a brief example of how you can you can analyse alignments (pileup) inside a given set of regions. This example assumes that:

* You have a set of alignments as a bam or bed file
* You have a set of regions in another bed-file

In these scenarios it is common to have a bed-file that easily fits into memory (e.g. binding sites) and a large BAM file that we don't want to keep in memory. As long as these are sorted, BioNumPy will handle streaming of the files and make sure that only data for a single chromosome is processed in memory at once. We can specify the behaviour by using `stream=True/False` when reading the data.

We start by reading our genome as the rest of the data depends on the genome (chromosome sizes).

>>> import numpy as np
>>> import bionumpy as bnp
>>> genome = bnp.Genome.from_file("example_data/hg38.chrom.sizes")
>>> print(genome)
           Chromosome                     Size
                 chr1                248956422
                 chr2                242193529
                 chr3                198295559
                 chr4                190214555
                 chr5                181538259
                 chr6                170805979
                 chr7                159345973
                 chr8                145138636
                 chr9                138394717
                chr10                133797422

Using this genome object, we can now read our peaks and alignemnts. We add `stream=True` to tell BioNumPy to *stream* all data, i.e. only keep the data for a single chromosome in memory at the same time.

>>> peaks = genome.read_intervals("example_data/ctcf_chr21-22.bed.gz", stream=True)
>>> reads = genome.read_intervals("example_data/ctcf_chr21-22.bam", stream=True)

The `peaks` and `reads` are now intervals. Using the `.get_pileup()` method, we can create a `GenomeArray` pileup, which lets us easily query the pileup value at any position along the genome. A `GenomeArray` can be thought of as an array along the genome.

>>> read_pileup = reads.get_pileup()

We can then subset this pileup where we have peaks

>>> peaks_pileup = read_pileup[peaks]

... get e.g. the max pileup value inside each peak:

>>> max_pileup_value_per_peak = np.max(peaks_pileup, axis=1)

If you try to print `max_ipleup_value_per_peak`, you will see that no data is shown. Instead you get a ComputationNode object. This is because since we used `stream=True`, no computations have been performed yet. BioNumPy waits until you ask for results by calling `.compute()` on a ComputationNode object, and then finds out how to correctly stream the data and do the necessary jobs for creating the data you ask for.

Having the max peak values, we can easily filter the original peaks, e.g. get the peaks with max pileup value above some threshold:

>>> high_peaks = peaks[max_pileup_value_per_peak > 4]

If we now call `.compute()` on the final peaks, BioNumPy performs all the necessary steps and gives us a final set of peaks:

>>> print(high_peaks.compute())
    Genomic Intervals on ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', '...']:
    Interval with 362 entries
                   chromosome                    start                     stop
                        chr21                 13980409                 13980679
                        chr21                 37460251                 37460521
                        chr21                 31933717                 31933987
                        chr21                 29164384                 29164654
                        chr21                 34351806                 34352076
                        chr21                 25708519                 25708789
                        chr21                 32403198                 32403468
                        chr21                 30434098                 30434368
                        chr21                 30044653                 30044724
                        chr21                 31874156                 31874275
