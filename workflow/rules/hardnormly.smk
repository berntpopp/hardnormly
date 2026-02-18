rule run_hardnormly_pipeline:
    input:
        vcf=lambda wildcards: JOBS[wildcards.vcf_basename],
    output:
        vcf=os.path.join(OUTPUT_DIR, "{vcf_basename}.hardnormly.vcf.gz"),
        idx=os.path.join(OUTPUT_DIR, "{vcf_basename}.hardnormly.vcf.gz.tbi"),
    log:
        os.path.join(LOG_DIR, "{vcf_basename}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "{vcf_basename}.tsv")
    params:
        hardnormly_script=get_hardnormly_script(),
        fasta=REF,
        include_beds=" ".join(f"-b {bed}" for bed in INCLUDE_BEDS),
        exclude_beds=" ".join(f"-e {bed}" for bed in EXCLUDE_BEDS),
        caller=f"--caller {CALLER}" if CALLER else "",
        filters_file=f"--filters-file {FILTERS_FILE}" if FILTERS_FILE else "",
        strip_annotations=f"--strip-annotations {STRIP_ANNOTATIONS}" if STRIP_ANNOTATIONS else "",
        genome_file=f"-g {GENOME_FILE}" if GENOME_FILE else "",
        slop=SLOP,
        generate_stats="--generate-stats" if GENERATE_STATS else "",
        auto_index="--auto-index" if AUTO_INDEX else "",
        only_pass="--only-pass" if ONLY_PASS else "",
        plot_stats="--plot-stats" if PLOT_STATS else "",
        plot_output_dir=lambda wildcards: f"--plot-output-dir {os.path.join(PLOT_OUTPUT_DIR, wildcards.vcf_basename)}" if PLOT_STATS else "",
    conda:
        "../envs/hardnormly.yaml"
    shell:
        """
        {params.hardnormly_script} \
            -v {input.vcf} \
            -f {params.fasta} \
            {params.include_beds} \
            {params.exclude_beds} \
            {params.caller} \
            {params.filters_file} \
            {params.strip_annotations} \
            {params.genome_file} \
            --slop {params.slop} \
            {params.generate_stats} \
            {params.auto_index} \
            {params.only_pass} \
            {params.plot_stats} \
            {params.plot_output_dir} \
            -o {output.vcf} \
            2> {log}
        """
