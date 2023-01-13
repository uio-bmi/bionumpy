rule callpeak:
    input:
        treatment="results/intervals/{a}_reads.bed",   # required: treatment sample(s)
        # control="samples/b.bam"      # optional: control sample(s)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs2/{a}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed",
                 '_control_lambda.bdg',
                 '_treat_pileup.bdg'
                 )
    log:
        "logs/macs2/callpeak_{a}.log"
    params:
        "-f BED -g 6000 --nomodel --extsize 200 -p 0.001 --bdg"
    wrapper:
        "v1.21.1-3-g1bd26948/bio/macs2/callpeak"

rule bionumpy_callpeak:
    input:
        treatment="results/intervals/{a}_reads.bed",
        chrom_sizes='results/small.chrom.sizes'        
    output:
        'results/bionumpy/callpeak/{a}.narrowPeak'
    shell:
        'python ../scripts/macs2.py {input.treatment} {input.chrom_sizes} --fragment-length 200 --p-value-cutoff 0.001 --outfilename {output}'

rule simulate_chipseq:
    output:
        "results/small.chrom.sizes",
        'results/intervals/small_reads.bed'
    script:
        '../scripts/simulate_chipseq.py'
