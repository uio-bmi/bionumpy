import bionumpy as bnp
from bionumpy.alignments import alignment_to_interval
import typer

def bam2bed(input_file: str, output_file: str):
    alignments_iter = bnp.open(input_file).read_chunks()
    with bnp.open(output_file, "w") as f:
        for alignments in alignments_iter:
            f.write(alignment_to_interval(alignments))

if __name__ == "__main__":
    typer.run(bam2bed)