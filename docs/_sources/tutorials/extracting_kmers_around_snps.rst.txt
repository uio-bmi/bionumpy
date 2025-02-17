Extracting kmers around snps
-----------------------------

This is a small example of how you can use BioNumPy to extract kmers around SNPs.

We have a set of SNPs in a vcf-file and a reference genome, and for every SNP we want to extract kmers
that "overlap" with that SNP both for the reference allele and the variant allele.

For instance, assume we have a reference genome `AAAATAAA` and a SNP `T/G` at position 4.

Then BioNumPy lets us crate windows around this SNP and extract sequences from those windows. If `k=3`, we would expect to get the reference kmers
`AAA`, `ATA` and `TAA` and the "variant allele" kmers kmers `AAG`, `AGA`, `GAA`.

The following is code for doing this. First we read the reference genome and variants, and filter so that we only have SNPs and not indels:

.. testcode::

    import bionumpy as bnp
    k = 5

    # Read genome and variants
    genome = bnp.Genome.from_file("example_data/sacCer3.fa", filter_function=None)
    variants = genome.read_locations("example_data/sacCer3_sample_variants.vcf.gz", has_numeric_chromosomes=False)

    # get only snps from vcf
    snp_mask = (variants.get_data_field('ref_seq').shape[-1] == 1) & (variants.get_data_field('alt_seq').shape[-1] == 1)
    variants = variants[snp_mask]

We  then extract windows around the variants and reads sequences from the reference genome in these windows:

.. testcode::

    # Get windows around these snps
    windows = variants.get_windows(flank=k-1)

    # Use the windows to extract sequences (kmers)
    sequences = genome.read_sequence()[windows]
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    reference_kmers = bnp.get_kmers(sequences, k)

Now `reference_kmers` is a RaggedArray where each row contains kmers for a variant. For instance, the reference kmers for the first three variants are:

.. testcode::

 print(reference_kmers[0:3])

.. testoutput::

    [TTACC, TACCC, ACCCA, CCCAT, CCATA]
    [ATAAC, TAACG, AACGC, ACGCC, CGCCC]
    [TTTTG, TTTGA, TTGAT, TGATA, GATAT]


For getting kmers from the variant alleles, we replace the allele bases with the variant allele bases:

.. testcode::

    # Extract kmers with alternative allele on SNPs
    sequences[:, k-1] = variants.get_data_field('alt_seq').ravel()
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    alt_kmers = bnp.get_kmers(sequences, k)
    print(alt_kmers[0:3])

The `alt_kmers` now contains kmers for the alt alleles for each variant.

.. testoutput::

    [TTACG, TACGC, ACGCA, CGCAT, GCATA]
    [ATAAA, TAAAG, AAAGC, AAGCC, AGCCC]
    [TTTTA, TTTAA, TTAAT, TAATA, AATAT]

Note that these kmers are encoded using a `KmerEncoding`. We can get the raw encoded kmers (int64) values:

.. testcode::

    raw_kmers = alt_kmers.raw()
    print(raw_kmers[0:5])

.. testoutput::

    [591 403 100 793 198]
    [ 12 515 384 352 344]
    [255  63 783 195 816]
    [ 66 528 132 801 712]
    [912 228 825 718 435]

We can also call `.ravel()` on the RaggedArray to get one flat array with all the kmers, and write these to file.
We use `to_string()` on each kmer to get the string kmer. Note that this works well on small datasets, but for large
amounts of kmers, iterating each kmer and calling `to_string()` will be slow:

.. testcode::

    with open("my_kmers.txt", "w") as f:
        # write the first 100 kmers to file
        kmers_as_strings = (kmer.to_string() for kmer in reference_kmers.ravel()[0:100])
        f.writelines((k + "\n" for k in kmers_as_strings))

