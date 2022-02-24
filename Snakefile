SEEDS=["123", "124", "125", "126", "127", "128", "129"]
UNPAIRED_TOOLS=["Minimap2_nanopore", "Kraken2Minimap2_nanopore"]
PAIRED_TOOLS=["Bowtie2_end-to-end", "Bowtie2_local", "HISAT2", "Kraken2", "BBMap", "Kraken2Bowtie2", "Kraken2HISAT2"]
MICROBIOME_TYPES=["gut", "oral"]

"""
# draw distribution plot with matplotlib (runtime)
rule draw_distribution_runtime:
    input:
        results="{tech}_reads/{seed}/results/{tool}/{type}_microbiome.txt"
    output:
        ""
    run:
        import json
"""

rule all:
    threads:
        workflow.cores
    input:
        expand("hiseq_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_paired.png", seed=SEEDS, tool=PAIRED_TOOLS, type=MICROBIOME_TYPES),
        expand("miseq_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_paired.png", seed=SEEDS, tool=PAIRED_TOOLS, type=MICROBIOME_TYPES),
        expand("nanopore_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_unpaired.png", seed=SEEDS, tool=UNPAIRED_TOOLS, type=MICROBIOME_TYPES)

# draw confusion matrix with matplotlib
rule draw_confusion_matrix:
    threads:
        workflow.cores
    input:
        results="{tech}_reads/{seed}/results/{tool}/{type}_microbiome_{paired_unpaired}.txt"
    output:
        results="{tech}_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_{paired_unpaired}.png"
    run:
        import json
        from mlxtend.plotting import plot_confusion_matrix
        import matplotlib.pyplot as plt
        import numpy as np
        with open(input.results) as f:
            confusion_matrix = json.loads(f.readlines()[0])
            print(confusion_matrix)

            cm = np.array([[confusion_matrix['TP'], confusion_matrix['FN']],
                           [confusion_matrix['FP'], confusion_matrix['TN']]])

            classes = ['human', 'other']

            figure, ax = plot_confusion_matrix(conf_mat=cm,
                                               class_names=classes,
                                               show_absolute=True,
                                               show_normed=True,
                                               colorbar=True)
            plt.tight_layout()
            plt.savefig(output.results)


# compare to original files with fastq-compare and write results (confusion_matrix) to file
rule fastq_compare:
    threads:
        workflow.cores
    input:
        f1="{tech}_reads/{seed}/{type}_microbiome_R1_{paired_unpaired}.fastq.gz",
        f2="{tech}_reads/{seed}/results/{tool}/{type}_microbiome_clean_1_{paired_unpaired}.fastq.gz"
    output:
        results="{tech}_reads/{seed}/results/{tool}/{type}_microbiome_{paired_unpaired}.txt"
    run:
        import fastq_compare.compare as compare
        import json
        entries1, total1 = compare.count(input.f1, '_')
        entries2, total2 = compare.count(input.f2, '_')
        print(entries1, total1)
        print(entries2, total2)
        confusion_matrix = {'TP': 0, # removed human
                            'FP': 0, # removed other
                            'TN': 0, # not removed other
                            'FN': 0, # not removed human
                            'T1': total1, # total of file 1
                            'T2': total2} # total of file 2
        # if all reads of one class were removed, set to 0 (no remaining)
        if 'human' not in entries2: entries2['human'] = 0
        if 'other' not in entries2: entries2['other'] = 0
        # 
        confusion_matrix['TP'] = entries1['human'] - entries2['human']
        confusion_matrix['FP'] = entries1['other'] - entries2['other']
        confusion_matrix['TN'] = entries2['other']
        confusion_matrix['FN'] = entries2['human']
        with open(output.results, 'w') as f:
            f.write(json.dumps(confusion_matrix))

rule hocort_map_Kraken2Minimap2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/Kraken2/",
        "indexes/Minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Minimap2_{wildcards.tech_brand} && "
        "hocort map Kraken2Minimap2 -t {workflow.cores} -k indexes/Kraken2/human -m indexes/Minimap2_{wildcards.tech_brand}/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} -p {wildcards.tech_brand} && "
        "rm {params.R2}"

rule hocort_map_Kraken2Minimap2_unpaired:
    input:
        "indexes/Kraken2/",
        "indexes/Minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Minimap2_{wildcards.tech_brand} && "
        "hocort map Kraken2Minimap2 -t {workflow.cores} -k indexes/Kraken2/human -m indexes/Minimap2_{wildcards.tech_brand}/human -i {input.R1} -o {output.R1} -p {wildcards.tech_brand}"

rule hocort_map_Kraken2Bowtie2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/Bowtie2/",
        "indexes/Kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Bowtie2 && "
        "hocort map Kraken2Bowtie2 -t {workflow.cores} -b indexes/Bowtie2/human -k indexes/Kraken2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Kraken2Bowtie2_unpaired:
    input:
        "indexes/Bowtie2/",
        "indexes/Kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Bowtie2 && "
        "hocort map Kraken2Bowtie2 -t {workflow.cores} -b indexes/Bowtie2/human -k indexes/Kraken2/human -i {input.R1} -o {output.R1}"

rule hocort_map_Kraken2HISAT2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/HISAT2/",
        "indexes/Kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2HISAT2 && "
        "hocort map Kraken2HISAT2 -t {workflow.cores} -s indexes/HISAT2/human -k indexes/Kraken2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Kraken2HISAT2_unpaired:
    input:
        "indexes/HISAT2/",
        "indexes/Kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2HISAT2 && "
        "hocort map Kraken2HISAT2 -t {workflow.cores} -s indexes/HISAT2/human -k indexes/Kraken2/human -i {input.R1} -o {output.R1}"

# Kraken2 does not compress output
rule hocort_map_Kraken2_paired:
    params:
        out="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean#_paired.fastq.gz",
        R2="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/Kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2 && "
        "hocort map Kraken2 -t {workflow.cores} -x indexes/Kraken2/human -i {input.R1} {input.R2} -o {params.out} && "
        "rm {params.R2}"

# Kraken2 does not compress output
rule hocort_map_Kraken2_unpaired:
    input:
        "indexes/Kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2 && "
        "hocort map Kraken2 -t {workflow.cores} -x indexes/Kraken2/human -i {input.R1} -o {output.R1}"

rule hocort_map_Minimap2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/Minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Minimap2_{wildcards.tech_brand} && "
        "hocort map Minimap2 -t {workflow.cores} -x indexes/Minimap2_{wildcards.tech_brand}/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} -p {wildcards.tech_brand} && "
        "rm {params.R2}"

rule hocort_map_Minimap2_unpaired:
    input:
        "indexes/Minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Minimap2_{wildcards.tech_brand} && "
        "hocort map Minimap2 -t {workflow.cores} -x indexes/Minimap2_{wildcards.tech_brand}/human -i {input.R1} -o {output.R1} -p {wildcards.tech_brand}"

rule hocort_map_HISAT2_paired:
    params:
        R2="{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/HISAT2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/HISAT2 && "
        "hocort map HISAT2 -t {workflow.cores} -x indexes/HISAT2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_HISAT2_unpaired:
    input:
        "indexes/HISAT2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/HISAT2 && "
        "hocort map HISAT2 -t {workflow.cores} -x indexes/HISAT2/human -i {input.R1} -o {output.R1}"

rule hocort_map_BWA_MEM2_paired:
    params:
        R2="{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/BWA_MEM2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BWA_MEM2 && "
        "hocort map BWA_MEM2 -t {workflow.cores} -x indexes/BWA_MEM2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_BWA_MEM2_unpaired:
    input:
        "indexes/BWA_MEM2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BWA_MEM2 && "
        "hocort map BWA_MEM2 -t {workflow.cores} -x indexes/BWA_MEM2/human -i {input.R1} -o {output.R1}"

rule hocort_map_BBMap_paired:
    params:
        R2="{tech}_reads/{seed}/results/BBMap/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/BBMap/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BBMap/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BBMap/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BBMap && "
        "hocort map BBMap -t {workflow.cores} -x indexes/BBMap/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_BBMap_unpaired:
    input:
        "indexes/BBMap/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BBMap/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BBMap/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BBMap && "
        "hocort map BBMap -t {workflow.cores} -x indexes/BBMap/human -i {input.R1} -o {output.R1}"

rule hocort_map_Bowtie2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/Bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2_{wildcards.mode} && "
        "hocort map Bowtie2 -m {wildcards.mode} -t {workflow.cores} -x indexes/Bowtie2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Bowtie2_unpaired:
    input:
        "indexes/Bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_benchmark.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2_{wildcards.mode} && "
        "hocort map Bowtie2 -m {wildcards.mode} -t {workflow.cores} -x indexes/Bowtie2/human -i {input.R1} -o {output.R1}"

rule hocort_index_gen_minimap2:
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("indexes/Minimap2_{tech}/")
    threads:
        workflow.cores
    shell:
        "mkdir -p indexes/Minimap2_{wildcards.tech}/ && "
        "hocort index Minimap2 -i {input.genome} -o indexes/Minimap2_{wildcards.tech}/human -p {wildcards.tech}"

rule hocort_index_gen:
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("indexes/{mapper}/")
    wildcard_constraints:
        mapper="^(Bowtie2|HISAT2|Kraken2|BBMap)"
    threads:
        workflow.cores
    shell:
        "mkdir -p indexes/{wildcards.mapper}/ && "
        "hocort index {wildcards.mapper} -i {input.genome} -o indexes/{wildcards.mapper}/human"

# mix reads to simulate oral microbiome (50% human/host, 50% other)
rule simulate_oral_microbiome_paired:
    params:
        temp_R1="{tech}_reads/{seed}/oral_microbiome_R1_paired.fastq",
        temp_R2="{tech}_reads/{seed}/oral_microbiome_R2_paired.fastq"
    input:
        R1_human="{tech}_reads/{seed}/human/{tech}_human_reads_R1.fastq.gz",
        R1_other="{tech}_reads/{seed}/other/{tech}_other_reads_R1.fastq.gz",
        R2_human="{tech}_reads/{seed}/human/{tech}_human_reads_R2.fastq.gz",
        R2_other="{tech}_reads/{seed}/other/{tech}_other_reads_R2.fastq.gz"
    threads:
        workflow.cores
    output:
        R1="{tech}_reads/{seed}/oral_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/oral_microbiome_R2_paired.fastq.gz"
    shell:
        "fastq-mixshuffle -1 {input.R1_human} {input.R1_other} -2 {input.R2_human} {input.R2_other} -o {params.temp_R1} {params.temp_R2} -a 50 50 -k 5000000 && "
        "gzip {params.temp_R1} && "
        "gzip {params.temp_R2}"

# mix reads to simulate oral microbiome (50% human/host, 50% other)
rule simulate_oral_microbiome_unpaired:
    params:
        temp_R1="{tech}_reads/{seed}/oral_microbiome_R1_unpaired.fastq"
    input:
        R1_human="{tech}_reads/{seed}/human/{tech}_human_reads_R1.fastq.gz",
        R1_other="{tech}_reads/{seed}/other/{tech}_other_reads_R1.fastq.gz"
    threads:
        workflow.cores
    output:
        R1="{tech}_reads/{seed}/oral_microbiome_R1_unpaired.fastq.gz"
    shell:
        "fastq-mixshuffle -1 {input.R1_human} {input.R1_other} -o {params.temp_R1} -a 50 50 -k 5000000 && "
        "gzip {params.temp_R1}"

# mix reads to simulate gut microbiome (20% human/host, 80% other)
rule simulate_gut_microbiome_paired:
    params:
        temp_R1="{tech}_reads/{seed}/gut_microbiome_R1_paired.fastq",
        temp_R2="{tech}_reads/{seed}/gut_microbiome_R2_paired.fastq"
    input:
        R1_human="{tech}_reads/{seed}/human/{tech}_human_reads_R1.fastq.gz",
        R1_other="{tech}_reads/{seed}/other/{tech}_other_reads_R1.fastq.gz",
        R2_human="{tech}_reads/{seed}/human/{tech}_human_reads_R2.fastq.gz",
        R2_other="{tech}_reads/{seed}/other/{tech}_other_reads_R2.fastq.gz"
    threads:
        workflow.cores
    output:
        R1="{tech}_reads/{seed}/gut_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/gut_microbiome_R2_paired.fastq.gz"
    shell:
        "fastq-mixshuffle -1 {input.R1_human} {input.R1_other} -2 {input.R2_human} {input.R2_other} -o {params.temp_R1} {params.temp_R2} -a 1 99 -k 5000000 && "
        "gzip {params.temp_R1} && "
        "gzip {params.temp_R2}"

# mix reads to simulate gut microbiome (20% human/host, 80% other)
rule simulate_gut_microbiome_unpaired:
    params:
        temp_R1="{tech}_reads/{seed}/gut_microbiome_R1_unpaired.fastq"
    input:
        R1_human="{tech}_reads/{seed}/human/{tech}_human_reads_R1.fastq.gz",
        R1_other="{tech}_reads/{seed}/other/{tech}_other_reads_R1.fastq.gz"
    threads:
        workflow.cores
    output:
        R1="{tech}_reads/{seed}/gut_microbiome_R1_unpaired.fastq.gz"
    shell:
        "fastq-mixshuffle -1 {input.R1_human} {input.R1_other} -o {params.temp_R1} -a 1 99 -k 5000000 && "
        "gzip {params.temp_R1}"

rule generate_miseq_human_reads:
    params:
        temp_R1="miseq_reads/{seed}/human/miseq_human_reads_R1.fastq",
        temp_R2="miseq_reads/{seed}/human/miseq_human_reads_R2.fastq"
    threads:
        workflow.cores
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("miseq_reads/{seed}/human/"),
        "miseq_reads/{seed}/human/miseq_human_reads_abundance.txt",
        R1="miseq_reads/{seed}/human/miseq_human_reads_R1.fastq.gz",
        R2="miseq_reads/{seed}/human/miseq_human_reads_R2.fastq.gz"
    shell:
        "mkdir -p miseq_reads/{wildcards.seed}/human/ && "
        "iss generate --cpus {workflow.cores} --genomes {input.genome} --seed {wildcards.seed} --n_reads 6000000 --model miseq --output miseq_reads/{wildcards.seed}/human/miseq_human_reads && "
        "cutadapt --cores {workflow.cores} --max-n 0 -o {output.R1} -p {output.R2} {params.temp_R1} {params.temp_R2} && "
        "rm {params.temp_R1} {params.temp_R2} && "
        "fastq-rename -i {output.R1} -o {params.temp_R1} -b human_ -e /1 && "
        "fastq-rename -i {output.R2} -o {params.temp_R2} -b human_ -e /2 && "
        "rm {output.R1} {output.R2} && "
        "gzip {params.temp_R1} && "
        "gzip {params.temp_R2}"

rule generate_miseq_other_reads:
    params:
        temp_R1="miseq_reads/{seed}/other/miseq_other_reads_R1.fastq",
        temp_R2="miseq_reads/{seed}/other/miseq_other_reads_R2.fastq"
    threads:
        workflow.cores
    input:
        genome="bacteria_viral_genomes.fna"
    output:
        directory("miseq_reads/{seed}/other/"),
        "miseq_reads/{seed}/other/miseq_other_reads_abundance.txt",
        R1="miseq_reads/{seed}/other/miseq_other_reads_R1.fastq.gz",
        R2="miseq_reads/{seed}/other/miseq_other_reads_R2.fastq.gz"
    shell:
        "mkdir -p miseq_reads/{wildcards.seed}/other/ && "
        "iss generate --cpus {workflow.cores} --genomes {input.genome} --seed {wildcards.seed} --n_reads 12000000 --model miseq --output miseq_reads/{wildcards.seed}/other/miseq_other_reads && "
        "cutadapt --cores {workflow.cores} --max-n 0 -o {output.R1} -p {output.R2} {params.temp_R1} {params.temp_R2} && "
        "rm {params.temp_R1} {params.temp_R2} && "
        "fastq-rename -i {output.R1} -o {params.temp_R1} -b other_ -e /1 && "
        "fastq-rename -i {output.R2} -o {params.temp_R2} -b other_ -e /2 && "
        "rm {output.R1} {output.R2} && "
        "gzip {params.temp_R1} && "
        "gzip {params.temp_R2}"

rule generate_hiseq_human_reads:
    params:
        temp_R1="hiseq_reads/{seed}/human/hiseq_human_reads_R1.fastq",
        temp_R2="hiseq_reads/{seed}/human/hiseq_human_reads_R2.fastq"
    threads:
        workflow.cores
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("hiseq_reads/{seed}/human/"),
        "hiseq_reads/{seed}/human/hiseq_human_reads_abundance.txt",
        R1="hiseq_reads/{seed}/human/hiseq_human_reads_R1.fastq.gz",
        R2="hiseq_reads/{seed}/human/hiseq_human_reads_R2.fastq.gz"
    shell:
        "mkdir -p hiseq_reads/{wildcards.seed}/human/ && "
        "iss generate --cpus {workflow.cores} --genomes {input.genome} --seed {wildcards.seed} --n_reads 6000000 --model hiseq --output hiseq_reads/{wildcards.seed}/human/hiseq_human_reads && "
        "cutadapt --cores {workflow.cores} --max-n 0 -o {output.R1} -p {output.R2} {params.temp_R1} {params.temp_R2} && "
        "rm {params.temp_R1} {params.temp_R2} && "
        "fastq-rename -i {output.R1} -o {params.temp_R1} -b human_ -e /1 && "
        "fastq-rename -i {output.R2} -o {params.temp_R2} -b human_ -e /2 && "
        "rm {output.R1} {output.R2} && "
        "gzip {params.temp_R1} && "
        "gzip {params.temp_R2}"

rule generate_hiseq_other_reads:
    params:
        temp_R1="hiseq_reads/{seed}/other/hiseq_other_reads_R1.fastq",
        temp_R2="hiseq_reads/{seed}/other/hiseq_other_reads_R2.fastq"
    threads:
        workflow.cores
    input:
        genome="bacteria_viral_genomes.fna"
    output:
        directory("hiseq_reads/{seed}/other/"),
        "hiseq_reads/{seed}/other/hiseq_other_reads_abundance.txt",
        R1="hiseq_reads/{seed}/other/hiseq_other_reads_R1.fastq.gz",
        R2="hiseq_reads/{seed}/other/hiseq_other_reads_R2.fastq.gz"
    shell:
        "mkdir -p hiseq_reads/{wildcards.seed}/other/ && "
        "iss generate --cpus {workflow.cores} --genomes {input.genome} --seed {wildcards.seed} --n_reads 12000000 --model hiseq --output hiseq_reads/{wildcards.seed}/other/hiseq_other_reads && "
        "cutadapt --cores {workflow.cores} --max-n 0 -o {output.R1} -p {output.R2} {params.temp_R1} {params.temp_R2} && "
        "rm {params.temp_R1} {params.temp_R2} && "
        "fastq-rename -i {output.R1} -o {params.temp_R1} -b other_ -e /1 && "
        "fastq-rename -i {output.R2} -o {params.temp_R2} -b other_ -e /2 && "
        "rm {output.R1} {output.R2} && "
        "gzip {params.temp_R1} && "
        "gzip {params.temp_R2}"

rule generate_nanopore_human_reads:
    params:
        temp_R1="nanopore_reads/{seed}/human/temp",
        R1="nanopore_reads/{seed}/human/nanopore_human_reads_R1.fastq"
    threads:
        workflow.cores
    input:
        "nanopore_reads/human/nanopore_training/",
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("nanopore_reads/{seed}/human/"),
        R1="nanopore_reads/{seed}/human/nanopore_human_reads_R1.fastq.gz"
    shell:
        "mkdir -p nanopore_reads/{wildcards.seed}/human/ && "
        "simulator.py genome -t {workflow.cores} -rg {input.genome} --seed {wildcards.seed} -n 3000000 -o {params.temp_R1} --fastq -b guppy --perfect -c nanopore_reads/human/nanopore_training/training && "
        "fastq-rename -i {params.temp_R1}_aligned_reads.fastq -o {params.R1} -b human_ && "
        "rm {params.temp_R1}_aligned_reads.fastq && "
        "gzip {params.R1}"

rule generate_nanopore_other_reads:
    params:
        temp_R1="nanopore_reads/{seed}/other/temp",
        R1="nanopore_reads/{seed}/other/nanopore_other_reads_R1.fastq"
    threads:
        workflow.cores
    input:
        "nanopore_reads/human/nanopore_training/",
        genome="bacteria_viral_genomes.fna"
    output:
        directory("nanopore_reads/{seed}/other/"),
        R1="nanopore_reads/{seed}/other/nanopore_other_reads_R1.fastq.gz"
    shell:
        "mkdir -p nanopore_reads/{wildcards.seed}/other/ && "
        "simulator.py genome -t {workflow.cores} -rg {input.genome} --seed {wildcards.seed} -n 6000000 -o {params.temp_R1} --fastq -b guppy --perfect -c nanopore_reads/human/nanopore_training/training && "
        "fastq-rename -i {params.temp_R1}_aligned_reads.fastq -o {params.R1} -b other_ && "
        "rm {params.temp_R1}_aligned_reads.fastq && "
        "gzip {params.R1}"

rule generate_nanopore_read_profile:
    params:
        accession="ERR3279199"
    threads:
        workflow.cores
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("nanopore_reads/human/nanopore_training/"),
    shell:
        "mkdir -p nanopore_reads/human/ && "
        # download dataset for read analysis
        "prefetch {params.accession} --output-directory nanopore_reads/ && "
        "fastq-dump nanopore_reads/{params.accession}/{params.accession}.sra -O nanopore_reads/{params.accession}/ && "
        # characterization stage
        "read_analysis.py genome -t {workflow.cores} -rg {input.genome} -i nanopore_reads/{params.accession}/{params.accession}.fastq -o nanopore_reads/human/nanopore_training/training"

rule merge_bacteria_viral_genomes:
    params:
        in_dir="bacteria_viral_genomes/"
    threads:
        workflow.cores
    input:
        "bacteria_viral_genomes/"
    output:
        "bacteria_viral_genomes.fna"
    shell:
        "find {params.in_dir} -name '*.fna.gz' -print0 | xargs -0 cat > bacteria_viral_genomes.fna.gz && gunzip -d bacteria_viral_genomes.fna.gz"

rule download_bacteria_viral_genomes:
    params:
        outdir="bacteria_viral_genomes/"
    threads:
        workflow.cores
    input:
        "accession_ids.txt"
    output:
        directory("bacteria_viral_genomes/")
    shell:
        "ncbi-genome-download -N -p 4 --flat-output -v -P -F fasta -A {input} -o {params.outdir} bacteria,viral"

"""
# already done, commented out for reproducibility sake as the accession ids are uploaded
rule select_bacteria_viral_accession_ids:
    threads:
        workflow.cores
    input:
        "metadata.txt"
    output:
        "accession_ids.txt"
    shell:
        "shuf -n 100 -o {output} --random-source={input} {input}"

rule download_bacteria_viral_metadata:
    threads:
        workflow.cores
    output:
        "metadata.txt"
    shell:
        "ncbi-genome-download -N -p 4 -F fasta -R reference,representative -n bacteria,viral | awk '{{print $1}}' | tail -n +2 > {output}"
"""

rule download_human_genome:
    threads:
        workflow.cores
    params:
        outdir="human_genome/"
    output:
        "human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    shell:
        "ncbi-genome-download -N -p 4 --flat-output -v -P -F fasta -o {params.outdir} -A GCF_000001405.39 all && gunzip -d human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
