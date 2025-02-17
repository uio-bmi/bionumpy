Using BioNumPy with Biopython
==============================

Biopython is a powerful library for bioinformatics, and can in many typical workflows easily be used together with BioNumPy.

For instance, Biopython sequences can be converted to BioNumPy array-based sequence structures by converting them to strings:

.. code-block:: python

	import Bio
	import bionumpy as bnp

	seq = Bio.Seq.Seq("ATGCGTACGTA")
	seq = bnp.as_encoded_array(str(seq.seq))

This becomes more useful when using Biopython to e.g. download sequence data. The following example is based on section 5.3.1 in the Biopython cookbook:

.. code-block:: python

	from Bio import Entrez
	from Bio import SeqIO

	Entrez.email = "A.N.Other@example.com"
	with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="6273291") as handle:
		seq_record = SeqIO.read(handle, "fasta")

	bnp_sequence = bnp.as_encoded_array(str(seq_record.seq))


Note that BioNumPy becomes useful when working with large data. In such cases, converting sequences to BioNumPy arrays from Biopython Seq objects will be slow, and we instead
recommend reading data directly from files (e.g. fasta files).
