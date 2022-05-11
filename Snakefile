SEEDS=["123", "124", "125", "126", "127", "128", "129"]
UNPAIRED_TOOLS=["Minimap2_nanopore", "Kraken2Minimap2_nanopore", "Kraken2"]
PAIRED_TOOLS=["Bowtie2_end-to-end", "Bowtie2_local", "Bowtie2_ad_hoc_un_conc", "HISAT2", "Kraken2", "BBMap_default", "BBMap_fast", "BWA_MEM2", "Kraken2Bowtie2", "Kraken2HISAT2", "Bowtie2Bowtie2", "Bowtie2HISAT2", "Minimap2_illumina", "Kraken2Minimap2_illumina"]
MICROBIOME_TYPES=["gut", "oral"]


rule all:
    threads:
        workflow.cores
    input:
        expand("hiseq_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_paired.png", seed=SEEDS, tool=PAIRED_TOOLS, type=MICROBIOME_TYPES),
        expand("miseq_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_paired.png", seed=SEEDS, tool=PAIRED_TOOLS, type=MICROBIOME_TYPES),
        expand("nanopore_reads/{seed}/results/{tool}/{type}_microbiome_confusion_matrix_unpaired.png", seed=SEEDS, tool=UNPAIRED_TOOLS, type=MICROBIOME_TYPES),
        expand("hiseq_reads/{type}_microbiome_paired_table.txt", type=MICROBIOME_TYPES),
        expand("miseq_reads/{type}_microbiome_paired_table.txt", type=MICROBIOME_TYPES),
        expand("nanopore_reads/{type}_microbiome_unpaired_table.txt", type=MICROBIOME_TYPES),
        "hiseq_reads/bowtie2_endtoend_deconseq_comparison.txt",
        "miseq_reads/bowtie2_endtoend_deconseq_comparison.txt"

rule bowtie2_endtoend_deconseq_comparison:
    threads:
        workflow.cores
    params:
        deconseq_benchmark="{tech}_reads/123/results/deconseq_gut_unpaired/gut_microbiome_benchmark_unpaired.txt",
        bowtie2_benchmark="{tech}_reads/123/results/Bowtie2_end-to-end/gut_microbiome_benchmark_unpaired.txt"
    input:
        results_deconseq="{tech}_reads/123/results/deconseq_gut_unpaired/gut_microbiome_unpaired.txt",
        results_bowtie2="{tech}_reads/123/results/Bowtie2_end-to-end/gut_microbiome_unpaired.txt"
    output:
        results="{tech}_reads/bowtie2_endtoend_deconseq_comparison.txt"
    run:
        import json
        def parse(tech, tool):
            path_results = f'{tech}_reads/123/results/{tool}/gut_microbiome_unpaired.txt'
            path_benchmark = f'{tech}_reads/123/results/{tool}/gut_microbiome_benchmark_unpaired.txt'
            with open(path_results, 'r') as f:
                data = json.loads(f.readlines()[0])
            with open(path_benchmark, 'r') as f:
                bits = f.readlines()[1].split()
            data['runtime'] = float(bits[0])
            return data
        def get_results(tech, tools):
            results = {}
            for tool in tools:
                data = parse(tech, tool)
                results[tool] = data
            return results
        results = get_results(wildcards.tech,
                              ['Bowtie2_end-to-end',
                               'deconseq_gut_unpaired'])
        with open(output.results, "w") as f:
            f.write(json.dumps(results))

rule draw_latex_table:
    threads:
        workflow.cores
    input:
        hiseq_gut="hiseq_reads/gut_microbiome_paired_computed.txt",
        hiseq_oral="hiseq_reads/oral_microbiome_paired_computed.txt",
        miseq_gut="miseq_reads/gut_microbiome_paired_computed.txt",
        miseq_oral="miseq_reads/oral_microbiome_paired_computed.txt",
        nanopore_gut="nanopore_reads/gut_microbiome_unpaired_computed.txt",
        nanopore_oral="nanopore_reads/oral_microbiome_unpaired_computed.txt"
    output:
        hiseq_gut="hiseq_reads/gut_microbiome_paired_table.txt",
        hiseq_oral="hiseq_reads/oral_microbiome_paired_table.txt",
        miseq_gut="miseq_reads/gut_microbiome_paired_table.txt",
        miseq_oral="miseq_reads/oral_microbiome_paired_table.txt",
        nanopore_gut="nanopore_reads/gut_microbiome_unpaired_table.txt",
        nanopore_oral="nanopore_reads/oral_microbiome_unpaired_table.txt"
    run:
        import json
        def build_tables(path):
            limits = {}
            with open(path, 'r') as f:
                data = json.loads(f.readlines()[0])
                for tool in data:
                    for header in data[tool]:
                        if header not in limits: limits[header] = {}
                        for elem in data[tool][header]:
                            val = data[tool][header][elem]
                            if elem not in limits[header]:
                                limits[header][elem] = {}
                                limits[header][elem]['max'] = val
                                limits[header][elem]['min'] = val
                            if limits[header][elem]['max'] < val: limits[header][elem]['max'] = val
                            if limits[header][elem]['min'] > val: limits[header][elem]['min'] = val
            tables = {'table_headers': {},
                      'tables': {}}
            with open(path, 'r') as f:
                data = json.loads(f.readlines()[0])
                # extract header
                first = list(data)[0]
                for header in data[first]:
                    tables['table_headers'][header] = 'Pipeline'
                    tables['tables'][header] = ''
                    for elem in data[first][header]:
                        tables['table_headers'][header] += f' & {header}\_{elem}'
                    tables['table_headers'][header] += ' \\\ \n'
                for tool in data:
                    for header in data[tool]:
                        tables['tables'][header] += tool.replace('_', '\_')
                        for elem in data[tool][header]:
                            # if element is largest or smallest add color
                            val = data[tool][header][elem]
                            decimals = 4
                            # round
                            rounded_val = round(val, decimals)
                            rounded_max = round(limits[header][elem]['max'], decimals)
                            rounded_min = round(limits[header][elem]['min'], decimals)
                            if rounded_val == rounded_max:
                                string = ' & ' + '\cellcolor{cyan!35}' + str(rounded_val)
                            elif rounded_val == rounded_min:
                                string = ' & ' + '\cellcolor{red!35}' + str(rounded_val)
                            else:
                                string = f' & {rounded_val}'
                            tables['tables'][header] += string
                        tables['tables'][header] += ' \\\ \n'
            return tables
        def write_to_file(in_path, out_path):
            tables = build_tables(in_path)
            with open(out_path, 'w') as f:
                for header in tables['table_headers']:
                    f.write(f'{header}:\n')
                    f.write(tables['table_headers'][header])
                    f.write(tables['tables'][header])
                    f.write('\n')
        write_to_file(input.hiseq_gut, output.hiseq_gut)
        write_to_file(input.hiseq_oral, output.hiseq_oral)
        write_to_file(input.miseq_gut, output.miseq_gut)
        write_to_file(input.miseq_oral, output.miseq_oral)
        write_to_file(input.nanopore_gut, output.nanopore_gut)
        write_to_file(input.nanopore_oral, output.nanopore_oral)

rule compute_all:
    threads:
        workflow.cores
    input:
        hiseq=expand("hiseq_reads/{type}_microbiome_paired_parsed.txt", type=MICROBIOME_TYPES),
        miseq=expand("miseq_reads/{type}_microbiome_paired_parsed.txt", type=MICROBIOME_TYPES),
        nanopore=expand("nanopore_reads/{type}_microbiome_unpaired_parsed.txt", type=MICROBIOME_TYPES)
    output:
        hiseq_gut="hiseq_reads/gut_microbiome_paired_computed.txt",
        hiseq_oral="hiseq_reads/oral_microbiome_paired_computed.txt",
        miseq_gut="miseq_reads/gut_microbiome_paired_computed.txt",
        miseq_oral="miseq_reads/oral_microbiome_paired_computed.txt",
        nanopore_gut="nanopore_reads/gut_microbiome_unpaired_computed.txt",
        nanopore_oral="nanopore_reads/oral_microbiome_unpaired_computed.txt"
    run:
        import json
        import numpy as np
        # read data from files
        def parse(tech, type, paired_unpaired):
            path = f'{tech}_reads/{type}_microbiome_{paired_unpaired}_parsed.txt'
            with open(path, 'r') as f:
                data = json.loads(f.readlines()[0])
            return data
        # calculate means, averages
        def compute_vals(data):
            out = {}
            for d in data:
                arr = np.array(data[d])
                mean = np.mean(arr).item()
                std = np.std(arr).item()
                median = np.median(arr).item()
                min = np.amin(arr).item()
                max = np.amax(arr).item()
                out[d] = {'mean': mean,
                          'std': std,
                          'median': median,
                          'min': min,
                          'max': max}
            return out
        def compute(tech, tools, paired_unpaired):
            out = {}
            for microbiome_type in MICROBIOME_TYPES:
                tech_data = parse(tech, microbiome_type, paired_unpaired)
                out[microbiome_type] = {}
                for tool in tools:
                    data = tech_data[tool]
                    final = compute_vals(data)
                    out[microbiome_type][tool] = final
            return out
        def compute_all():
            out = {}
            hiseq = compute('hiseq', PAIRED_TOOLS, 'paired')
            out['hiseq'] = hiseq
            miseq = compute('miseq', PAIRED_TOOLS, 'paired')
            out['miseq'] = miseq
            nanopore = compute('nanopore', UNPAIRED_TOOLS, 'unpaired')
            out['nanopore'] = nanopore
            return out
        final = compute_all()
        # write to output file
        with open(output.hiseq_gut, 'w') as f:
            f.write(json.dumps(final['hiseq']['gut']))
        with open(output.hiseq_oral, 'w') as f:
            f.write(json.dumps(final['hiseq']['oral']))
        with open(output.miseq_gut, 'w') as f:
            f.write(json.dumps(final['miseq']['gut']))
        with open(output.miseq_oral, 'w') as f:
            f.write(json.dumps(final['miseq']['oral']))
        with open(output.nanopore_gut, 'w') as f:
            f.write(json.dumps(final['nanopore']['gut']))
        with open(output.nanopore_oral, 'w') as f:
            f.write(json.dumps(final['nanopore']['oral']))

# retrieve all datapoints
rule parse_all:
    threads:
        workflow.cores
    input:
        hiseq=expand("hiseq_reads/{seed}/results/{tool}/{type}_microbiome_paired.txt", seed=SEEDS, tool=PAIRED_TOOLS, type=MICROBIOME_TYPES),
        miseq=expand("miseq_reads/{seed}/results/{tool}/{type}_microbiome_paired.txt", seed=SEEDS, tool=PAIRED_TOOLS, type=MICROBIOME_TYPES),
        nanopore=expand("nanopore_reads/{seed}/results/{tool}/{type}_microbiome_unpaired.txt", seed=SEEDS, tool=UNPAIRED_TOOLS, type=MICROBIOME_TYPES)
    output:
        hiseq_gut="hiseq_reads/gut_microbiome_paired_parsed.txt",
        hiseq_oral="hiseq_reads/oral_microbiome_paired_parsed.txt",
        miseq_gut="miseq_reads/gut_microbiome_paired_parsed.txt",
        miseq_oral="miseq_reads/oral_microbiome_paired_parsed.txt",
        nanopore_gut="nanopore_reads/gut_microbiome_unpaired_parsed.txt",
        nanopore_oral="nanopore_reads/oral_microbiome_unpaired_parsed.txt"
    run:
        import json
        import numpy as np
        # read data from files
        def parse(tech, seed, tool, type, paired_unpaired):
            path_results = f'{tech}_reads/{seed}/results/{tool}/{type}_microbiome_{paired_unpaired}.txt'
            path_benchmark = f'{tech}_reads/{seed}/results/{tool}/{type}_microbiome_benchmark_{paired_unpaired}.txt'
            with open(path_results, 'r') as f:
                data = json.loads(f.readlines()[0])
            with open(path_benchmark, 'r') as f:
                bits = f.readlines()[1].split()
            data['runtime'] = float(bits[0])
            return data
        # read data for each seed
        def parse_tool(tech, seeds, tool, type, paired_unpaired):
            header = ['TP', 'FP', 'TN', 'FN', 'P', 'N', 'TPR', 'TNR', 'PPV', 'ACC', 'BA', 'runtime']
            out = {}
            for h in header:
                out[h] = []
            for seed in seeds:
                data = parse(tech, seed, tool, type, paired_unpaired)
                for h in header:
                    out[h].append(data[h])
            return out
        def compute(tech, tools, paired_unpaired):
            out = {}
            for microbiome_type in MICROBIOME_TYPES:
                out[microbiome_type] = {}
                for tool in tools:
                    data = parse_tool(tech, SEEDS, tool, microbiome_type, paired_unpaired)
                    out[microbiome_type][tool] = data
            return out
        def parse_all():
            out = {}
            hiseq = compute('hiseq', PAIRED_TOOLS, 'paired')
            out['hiseq'] = hiseq
            miseq = compute('miseq', PAIRED_TOOLS, 'paired')
            out['miseq'] = miseq
            nanopore = compute('nanopore', UNPAIRED_TOOLS, 'unpaired')
            out['nanopore'] = nanopore
            return out
        final = parse_all()
        with open(output.hiseq_gut, 'w') as f:
            f.write(json.dumps(final['hiseq']['gut']))
        with open(output.hiseq_oral, 'w') as f:
            f.write(json.dumps(final['hiseq']['oral']))
        with open(output.miseq_gut, 'w') as f:
            f.write(json.dumps(final['miseq']['gut']))
        with open(output.miseq_oral, 'w') as f:
            f.write(json.dumps(final['miseq']['oral']))
        with open(output.nanopore_gut, 'w') as f:
            f.write(json.dumps(final['nanopore']['gut']))
        with open(output.nanopore_oral, 'w') as f:
            f.write(json.dumps(final['nanopore']['oral']))

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
        import fastq_tools.fastq_compare.compare as compare
        import json
        entries1, total1 = compare.count(input.f1, '_')
        entries2, total2 = compare.count(input.f2, '_')
        print(entries1, total1)
        print(entries2, total2)
        # if all reads of one class were removed, set to 0 (no remaining)
        if 'human' not in entries2: entries2['human'] = 0
        if 'other' not in entries2: entries2['other'] = 0
        total_human = entries1['human']
        total_other = entries1['other']
        tp = entries1['human'] - entries2['human']
        fp = entries1['other'] - entries2['other']
        tn = entries2['other']
        fn = entries2['human']
        p = total_human
        n = total_other
        tpr = tp / p
        tnr = tn / n
        confusion_matrix = {'TP': tp, # removed human
                            'FP': fp, # removed other
                            'TN': tn, # not removed other
                            'FN': fn, # not removed human
                            'T1': total1, # total of file 1
                            'T2': total2, # total of file 2
                            'P': p, # total number of real positive cases
                            'N': n, # total number of real negative cases
                            'TPR': tpr, # sensitivity/recall/true positive rate
                            'TNR': tnr, # specificity/selectivity/true negative rate
                            'PPV': tp / (tp + fp), # precision/positive predictive value
                            'ACC': (tp + tn) / (p + n), # accuracy
                            'BA': (tpr + tnr) / 2} # balanced accuracy
        with open(output.results, 'w') as f:
            f.write(json.dumps(confusion_matrix))

rule hocort_map_Kraken2Minimap2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/kraken2/",
        "indexes/minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Minimap2_{wildcards.tech_brand} && "
        "hocort map kraken2minimap2 -t {workflow.cores} -k indexes/kraken2/human -m indexes/minimap2_{wildcards.tech_brand}/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} -p {wildcards.tech_brand} && "
        "rm {params.R2}"

rule hocort_map_Kraken2Minimap2_unpaired:
    input:
        "indexes/kraken2/",
        "indexes/minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Minimap2_{tech_brand}/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Minimap2_{wildcards.tech_brand} && "
        "hocort map kraken2minimap2 -t {workflow.cores} -k indexes/kraken2/human -m indexes/minimap2_{wildcards.tech_brand}/human -i {input.R1} -o {output.R1} -p {wildcards.tech_brand}"

rule hocort_map_Kraken2Bowtie2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bowtie2/",
        "indexes/kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Bowtie2 && "
        "hocort map kraken2bowtie2 -t {workflow.cores} -b indexes/bowtie2/human -k indexes/kraken2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Kraken2Bowtie2_unpaired:
    input:
        "indexes/bowtie2/",
        "indexes/kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2Bowtie2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2Bowtie2 && "
        "hocort map kraken2bowtie2 -t {workflow.cores} -b indexes/bowtie2/human -k indexes/kraken2/human -i {input.R1} -o {output.R1}"

rule hocort_map_Kraken2HISAT2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/hisat2/",
        "indexes/kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2HISAT2 && "
        "hocort map kraken2hisat2 -t {workflow.cores} -s indexes/hisat2/human -k indexes/kraken2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Kraken2HISAT2_unpaired:
    input:
        "indexes/hisat2/",
        "indexes/kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2HISAT2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2HISAT2 && "
        "hocort map kraken2hisat2 -t {workflow.cores} -s indexes/hisat2/human -k indexes/kraken2/human -i {input.R1} -o {output.R1}"

rule hocort_map_Bowtie2HISAT2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Bowtie2HISAT2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bowtie2/",
        "indexes/hisat2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2HISAT2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2HISAT2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2HISAT2 && "
        "hocort map bowtie2hisat2 -t {workflow.cores} -s indexes/hisat2/human -b indexes/bowtie2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Bowtie2HISAT2_unpaired:
    input:
        "indexes/bowtie2/",
        "indexes/hisat2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2HISAT2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2HISAT2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2HISAT2 && "
        "hocort map bowtie2hisat2 -t {workflow.cores} -s indexes/hisat2/human -b indexes/bowtie2/human -i {input.R1} -o {output.R1}"

rule hocort_map_Bowtie2Bowtie2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Bowtie2Bowtie2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2Bowtie2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2Bowtie2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2Bowtie2 && "
        "hocort map bowtie2bowtie2 -t {workflow.cores} -x indexes/bowtie2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Bowtie2Bowtie2_unpaired:
    input:
        "indexes/bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2Bowtie2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2Bowtie2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2Bowtie2 && "
        "hocort map bowtie2bowtie2 -t {workflow.cores} -x indexes/bowtie2/human -i {input.R1} -o {output.R1}"

# Kraken2 does not compress output
rule hocort_map_Kraken2_paired:
    params:
        out="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean#_paired.fastq.gz",
        R2="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2 && "
        "hocort map kraken2 -t {workflow.cores} -x indexes/kraken2/human -i {input.R1} {input.R2} -o {params.out} && "
        "rm {params.R2}"

# Kraken2 does not compress output
rule hocort_map_Kraken2_unpaired:
    input:
        "indexes/kraken2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Kraken2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Kraken2 && "
        "hocort map kraken2 -t {workflow.cores} -x indexes/kraken2/human -i {input.R1} -o {output.R1}"

rule hocort_map_Minimap2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Minimap2_{wildcards.tech_brand} && "
        "hocort map minimap2 -t {workflow.cores} -x indexes/minimap2_{wildcards.tech_brand}/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} -p {wildcards.tech_brand} && "
        "rm {params.R2}"

rule hocort_map_Minimap2_unpaired:
    input:
        "indexes/minimap2_{tech_brand}/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Minimap2_{tech_brand}/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Minimap2_{wildcards.tech_brand} && "
        "hocort map minimap2 -t {workflow.cores} -x indexes/minimap2_{wildcards.tech_brand}/human -i {input.R1} -o {output.R1} -p {wildcards.tech_brand}"

rule hocort_map_HISAT2_paired:
    params:
        R2="{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/hisat2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/HISAT2 && "
        "hocort map hisat2 -t {workflow.cores} -x indexes/hisat2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_HISAT2_unpaired:
    input:
        "indexes/hisat2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/HISAT2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/HISAT2 && "
        "hocort map hisat2 -t {workflow.cores} -x indexes/hisat2/human -i {input.R1} -o {output.R1}"

rule hocort_map_BWA_MEM2_paired:
    params:
        R2="{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bwamem2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BWA_MEM2 && "
        "hocort map bwamem2 -t {workflow.cores} -x indexes/bwamem2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_BWA_MEM2_unpaired:
    input:
        "indexes/bwamem2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BWA_MEM2/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BWA_MEM2 && "
        "hocort map bwamem2 -t {workflow.cores} -x indexes/bwamem2/human -i {input.R1} -o {output.R1}"

rule hocort_map_BBMap_default_paired:
    params:
        R2="{tech}_reads/{seed}/results/BBMap_default/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bbmap/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BBMap_default/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BBMap_default/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BBMap_default && "
        'hocort map bbmap -t {workflow.cores} -x indexes/bbmap/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} -c="" && '
        "rm {params.R2}"

rule hocort_map_BBMap_fast_paired:
    params:
        R2="{tech}_reads/{seed}/results/BBMap_fast/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bbmap/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BBMap_fast/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BBMap_fast/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BBMap_fast && "
        'hocort map bbmap -t {workflow.cores} -x indexes/bbmap/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} -c="fast=t" && '
        "rm {params.R2}"

rule hocort_map_BBMap_unpaired:
    input:
        "indexes/bbmap/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/BBMap/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/BBMap/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/BBMap && "
        "hocort map bbmap -p nanopore -t {workflow.cores} -x indexes/bbmap/human -i {input.R1} -o {output.R1}"

rule hocort_map_Bowtie2_paired:
    params:
        R2="{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_clean_2_paired.fastq.gz"
    input:
        "indexes/bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_clean_1_paired.fastq.gz"
    wildcard_constraints:
        mode='|'.join(['end-to-end', 'local'])
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2_{wildcards.mode} && "
        "hocort map bowtie2 -p {wildcards.mode} -t {workflow.cores} -x indexes/bowtie2/human -i {input.R1} {input.R2} -o {output.R1} {params.R2} && "
        "rm {params.R2}"

rule hocort_map_Bowtie2_unpaired:
    input:
        "indexes/bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_clean_1_unpaired.fastq.gz"
    wildcard_constraints:
        mode='|'.join(['end-to-end', 'local'])
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2_{mode}/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2_{wildcards.mode} && "
        "hocort map bowtie2 -p {wildcards.mode} -t {workflow.cores} -x indexes/bowtie2/human -i {input.R1} -o {output.R1}"

rule map_Bowtie2_ad_hoc_un_conc_paired:
    input:
        "indexes/bowtie2/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_paired.fastq.gz",
        R2="{tech}_reads/{seed}/{type}_microbiome_R2_paired.fastq.gz"
    output:
        R1="{tech}_reads/{seed}/results/Bowtie2_ad_hoc_un_conc/{type}_microbiome_clean_1_paired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/Bowtie2_ad_hoc_un_conc/{type}_microbiome_benchmark_paired.txt"
    shell:
        "mkdir -p {wildcards.tech}_reads/{wildcards.seed}/results/Bowtie2_ad_hoc_un_conc && "
        "bowtie2 -x indexes/bowtie2/human -p {workflow.cores} -1 {input.R1} -2 {input.R2} --un-conc-gz bowtie2_clean_temp > /dev/null && "
        "mv bowtie2_clean_temp.1 {output.R1} && "
        "rm bowtie2_clean_temp.2"

rule download_deconseq:
    threads:
        workflow.cores
    output:
        dir=directory("deconseq-standalone-0.4.3/")
    shell:
        "wget https://sourceforge.net/projects/deconseq/files/standalone/deconseq-standalone-0.4.3.tar.gz/download && "
        "tar -xvzf download && "
        "rm download && "
        "chmod +x deconseq-standalone-0.4.3/bwa64"

rule deconseq_index_gen:
    input:
        dir=rules.download_deconseq.output.dir,
        deconseq="deconseq-standalone-0.4.3",
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        dummy_output=directory("indexes/deconseq")
    threads:
        workflow.cores
    shell:
        "mkdir -p indexes/deconseq && "
        "mkdir -p deconseq-standalone-0.4.3/db && "
        "{input.deconseq}/bwa64 index -p {input.deconseq}/db/hs_ref_GRCh37 -a bwtsw {input.genome}"

rule deconseq_map_unpaired:
    params:
        deconseq="deconseq-standalone-0.4.3",
        out_dir="{tech}_reads/{seed}/results/deconseq_{type}_unpaired/",
        input_R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq",
        output_R1="{tech}_reads/{seed}/results/deconseq_{type}_unpaired/123_clean.fq"
    input:
        "indexes/deconseq/",
        R1="{tech}_reads/{seed}/{type}_microbiome_R1_unpaired.fastq.gz"
    output:
        directory("{tech}_reads/{seed}/results/deconseq_{type}_unpaired/"),
        R1="{tech}_reads/{seed}/results/deconseq_{type}_unpaired/{type}_microbiome_clean_1_unpaired.fastq.gz"
    threads:
        workflow.cores
    benchmark:
        "{tech}_reads/{seed}/results/deconseq_{type}_unpaired/{type}_microbiome_benchmark_unpaired.txt"
    shell:
        "mkdir -p {params.out_dir} && "
        "gunzip -dk {input.R1} && "
        "cd {params.deconseq} && "
        "perl deconseq.pl --dbs hsref -f ../{params.input_R1} --out_dir ../{params.out_dir} --id 123 && "
        "cd .. && "
        "rm {params.input_R1} && "
        "mv {params.output_R1} {output.R1}"

rule hocort_index_gen_minimap2:
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("indexes/minimap2_{tech}/")
    threads:
        workflow.cores
    shell:
        "mkdir -p indexes/minimap2_{wildcards.tech}/ && "
        "hocort index minimap2 -i {input.genome} -o indexes/minimap2_{wildcards.tech}/human -p {wildcards.tech}"

rule hocort_index_gen:
    input:
        genome="human_genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
    output:
        directory("indexes/{mapper}/")
    wildcard_constraints:
        mapper='|'.join(['bowtie2', 'hisat2', 'kraken2', 'bbmap', 'bwamem2'])
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
