Simulating sequence datasets
------------------------------

Simulating sequences is very straightforward in BioNumPy. Since a sequence arrays have underlying numeric representations that are easy to relate to, one can directly simulate such an underlying representation using any numeric simulation procedure. In addition, bionumpy provides convenience functions that allows to directly simulate sequence data without having to think about any underlying numeric representation. This is the focus of the current tutorial. The function simulate_sequences is an easy way to simulate a set of sequences, by simply specifying a desired alphabet and a dictionary with desired sequence ids as keys and the desired length of each such sequence as value (here simulating 20 sequences with length 10..30):

    >>> import numpy as np
    >>> rng = np.random.default_rng(seed=1)
    >>> from bionumpy.simulate import simulate_sequences
    >>> from bionumpy.sequence import match_string
    >>> named_seqs = simulate_sequences('ACGT', {f's{i}':10+i for i in range(20)}, rng)
    >>> named_seqs
    SequenceEntry with 20 entries
                         name                 sequence
                           s0               CGTTAATTAC
                           s1              TCCTCCGGAAT
                           s2             TTGTCCTACACT
                           s3            ACCTAGCATACCC
                           s4           ATGTAGCGTCGACT
                           s5          CGCACGCTCGTTCAG
                           s6         GTCCACGTTAGTCCTG
                           s7        GGGTTAAGTAGTTTAGT
                           s8       CACAATGTTTCCGCTATG
                           s9      CGCTTCCAGGTTTTTAACC



One can now easily do a variety of analyses on these simulated sequences, e.g. compute the GC content per simulated sequence:

    >>> seqs = named_seqs.sequence
    >>> gc_content_per_seq = np.mean((seqs=='C')|(seqs=='G'), axis=1)
    >>> gc_content_per_seq
    array([0.3       , 0.54545455, 0.41666667, 0.53846154, 0.5       ,
           0.66666667, 0.5625    , 0.35294118, 0.44444444, 0.47368421,
           0.5       , 0.38095238, 0.68181818, 0.30434783, 0.66666667,
           0.52      , 0.46153846, 0.44444444, 0.53571429, 0.51724138])



If desired, such computed values per sequence can easily be added back as an additional column of the bionumpy data structure:

    >>> named_seqs = named_seqs.add_fields({'gc':gc_content_per_seq}, {'gc':float})
    >>> named_seqs
    DynamicSequenceEntry with 20 entries
                         name                 sequence                       gc
                           s0               CGTTAATTAC                      0.3
                           s1              TCCTCCGGAAT       0.5454545454545454
                           s2             TTGTCCTACACT       0.4166666666666667
                           s3            ACCTAGCATACCC       0.5384615384615384
                           s4           ATGTAGCGTCGACT                      0.5
                           s5          CGCACGCTCGTTCAG       0.6666666666666666
                           s6         GTCCACGTTAGTCCTG                   0.5625
                           s7        GGGTTAAGTAGTTTAGT      0.35294117647058826
                           s8       CACAATGTTTCCGCTATG       0.4444444444444444
                           s9      CGCTTCCAGGTTTTTAACC      0.47368421052631576



We can also easily apply a variety of built-in bionumpy functionality on our simulated sequences:
    >>> ac_hits = match_string(seqs, "AC")
    >>> ac_hit_sums = np.sum(ac_hits,axis=1)
    >>> ac_hit_sums
    array([1, 0, 2, 2, 1, 1, 1, 0, 1, 1, 1, 1, 2, 1, 2, 2, 2, 0, 1, 1])
    >>> named_seqs = named_seqs.add_fields({'ac_hits':ac_hit_sums}, {'ac_hits':int})
    >>> named_seqs
    DynamicSequenceEntry with 20 entries
                         name                 sequence                       gc                  ac_hits
                           s0               CGTTAATTAC                      0.3                        1
                           s1              TCCTCCGGAAT       0.5454545454545454                        0
                           s2             TTGTCCTACACT       0.4166666666666667                        2
                           s3            ACCTAGCATACCC       0.5384615384615384                        2
                           s4           ATGTAGCGTCGACT                      0.5                        1
                           s5          CGCACGCTCGTTCAG       0.6666666666666666                        1
                           s6         GTCCACGTTAGTCCTG                   0.5625                        1
                           s7        GGGTTAAGTAGTTTAGT      0.35294117647058826                        0
                           s8       CACAATGTTTCCGCTATG       0.4444444444444444                        1
                           s9      CGCTTCCAGGTTTTTAACC      0.47368421052631576                        1

