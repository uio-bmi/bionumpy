Genomic Data
============

Analysis of genomic data is best done with the knowledge of how the genome looks like. At the very least, we should know how the names of each chromosome in the genome, and how long the sequence of each chromosome is. This serves several purposes:

* We (and the computer) become aware of which reference genome we are using, which can prevent errors from mismatched reference genomes
* Chromosomes without any data gets handled correctly
* We get a speed increase in the algorithms

For these reasons, we always recomend using the genomic data framework for handling genomic data in BioNumPy, not using raw bed, vcf or bdg files directly from the `bnp.open` interface. The main access point of the genomic data interface is the genome object, which can be created directly from a `dict` of chromosome names and sizes; or, more commonly, from a `.chrom.sizes` or `.fa.fai` file, which are tab separated lines of chromosome names and sizes: 

.. testcode::

    import bionumpy as bnp
    genome = bnp.Genome.from_file("example_data/small.chrom.sizes")
    print(genome)

.. testoutput::

                   Chromosome                     Size
                         chr1                    10000
                         chr2                    20000
                         chr3                    30000

From this genome object we can create genomic data from reading files, or from bnpdataclass objects directly. For instance, we can now read a set of intervals from a `.bed` file and get a `GenomicIntervals` object. Because  `GenomicIntervals` objects knows which reference genome they belong to, they are easier and faster to work with than raw 'bnp.Interval' objects. 

.. testcode::

    intervals = genome.read_intervals('example_data/small_peaks.narrowPeak')
    print(intervals)
    mask = intervals.get_mask()
    print(mask)

.. testoutput::

    Genomic Intervals on ['chr1', 'chr2', 'chr3']:
   Interval with 13 entries
                   chromosome                    start                     stop
                         chr1                      482                      790
                         chr1                     5873                     6183
                         chr1                     6969                     7289
                         chr2                      691                     1007
                         chr2                     6155                     6459
                         chr2                     8337                     8641
                         chr2                    11170                    11481
                         chr2                    12375                    12677
                         chr2                    12918                    13232
                         chr2                    18804                    19134
    chr1: [ False False False ... False False False]
    chr2: [ False False False ... False False False]
    chr3: [ False False False ... False False False]

The mask here is a boolean array which is `True` wherever any of the intervals overlap. We can also get a `GenomicArray` from a bedgraph file by using `genome.read_track()`

.. testcode::

    pileup = genome.read_track('example_data/small_treat_pileup.bdg')
    print(pileup)

.. testoutput::

    chr1: [ 0.0 0.0 1.0 ... 0.0 0.0 0.0]
    chr2: [ 0.0 0.0 0.0 ... 0.0 0.0 0.0]
    chr3: [ 0.0 0.0 0.0 ... 0.0 0.0 0.0]
   
A benefit of working with genomic data in this way is that they cooperate in a consistent way. We can use GenomicIntervals as indexes to GenomicData, and can again use the extracted data to filter the intervals. Let's see how the `treat_pileup` looks in peak areas (we get the max and mean pileup value for each peak):

.. testcode::

    peak_pileups = pileup[intervals]
    print(peak_pileups.max(axis=-1))
    print(peak_pileups.mean(axis=-1))

.. testoutput::

    [ 227.  231.  412.  296.  165.  163.  271.  148.  268.  568. 1901.   90.
      236.]
    [145.68181818 147.86129032 256.671875   187.71202532 106.79934211
     106.23355263 173.04501608  96.27152318 169.59235669 346.37272727
     593.44084507  60.25614035 150.72580645] 

we can now again use these values to filter the intervals based on the treat_pileup values in each interval. For instance only keep peaks with a max treatment_value above 200:


.. testcode::

    high_pileup_peaks = intervals[peak_pileups.max(axis=-1)>200]
    print(high_pileup_peaks)

.. testoutput::

   Genomic Intervals on ['chr1', 'chr2', 'chr3']:
    Interval with 9 entries
                   chromosome                    start                     stop
                         chr1                      482                      790
                         chr1                     5873                     6183
                         chr1                     6969                     7289
                         chr2                      691                     1007
                         chr2                    11170                    11481
                         chr2                    12918                    13232
                         chr2                    18804                    19134
                         chr3                    10677                    11387
                         chr3                    27057                    27367

A further way to analyze these peaks is to check the sequence in the peaks for motifs. We can load the reference sequence using the `genome.read_sequence()` method:


.. testcode::

    genome_sequence = genome.read_sequence('example_data/small_sequence.fa')
    print(genome_sequence)
   

.. testoutput::

   GenomicSequence over chromosomes: ['chr1', 'chr2', 'chr3']

Now we can use our intervals as indexed to the reference sequence in much the same way as with genomic arrays. This we can use to get the sequences of the peaks and check them for motifs:

.. testcode::
   :skipif: True

    peak_sequences = genome_sequence[high_pileup_peaks]
    print(peak_sequences)
    from pyjaspar import jaspardb
    from bionumpy.sequence import PWM
    jaspar_object = jaspardb(release="JASPAR2020")
    ctcf_motif = jaspar_object.fetch_motifs_by_name('CTCF')[0]
    motif = PWM.from_dict(ctcf_motif.pwm)
    print(bnp.get_motif_scores(peak_sequences, motif).max(axis=-1))

.. testoutput::
   :skipif: True
    
    AGAGCCGGACCGAATGACGT...
    TCAGGTAGAACTCGCATTTC...
    AAGACTTATTTGATGGCCGG...
    ATAGAGAGCGTTCGGCGCTA...
    CTCCTTAGCATACAAACGGG...
    TGCCCCATCTCTACACAATT...
    GCGCTGCCGTCACGGCGGGG...
    TTACATCCTGACGAAATACA...
    GACGGGTAGGCGATTTTTAT...
    [-0.82852725 -4.78231564 -0.90188366 -4.97871502 -3.06703885 -1.00359223
      7.75920829  0.19382436  2.5188205 ]
   
Further reading
-----------------
    * :ref:`Tutorail: Genomic data from multiple sources <multi_genomic_data>`
    * :ref:`Tutorial: Analysing the read pileup within peaks (intervals) <subsetting_bed>`
    * :ref:`Tutorial: Computing the similarity between to bed files <similarity_measures_tutorial>`
    * :ref:`More about intervals <intervals>`
    * :ref:`API documentation on the arithmetics module <arithmetics_api>`

