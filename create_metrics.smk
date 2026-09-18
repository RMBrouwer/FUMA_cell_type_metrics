import os
import anndata

levels_per_dataset=dict()
all_levels = []

for i in config["scrnas"]:
   adata = anndata.read_h5ad(config["input"] + "/" + i + "/" + i + ".h5ad",backed='r')
   levels_per_dataset[i]=[col[10:] for col in adata.obs.columns if 'cell_type_level' in col]

def mem_from_input(multiplier=3, minimum=2000):
    def _mem(wildcards, input):
        size_mb = os.path.getsize(input.h5ad) / 1024**2
        return max(minimum, int(size_mb * multiplier))
    return _mem

rule all:
    input:
        [ os.path.join(config["metrics"], f"cellex/{id}/{lvl}/{id}.{lvl}.esmu_fmt.tsv") for id in config["scrnas"] for lvl in levels_per_dataset[id]] +
        [ os.path.join(config["metrics"], f"ewce/{id}/{lvl}/{id}.{lvl}_spec_ewce_top10.txt") for id in config["scrnas"] for lvl in levels_per_dataset[id]] +
        [ os.path.join(config["metrics"], f"cepo/{id}/{lvl}/{id}.{lvl}_cepo_top10.txt") for id in config["scrnas"] for lvl in levels_per_dataset[id]] +  
        [ os.path.join(config["metrics"], f"fuma/{id}/{lvl}/{id}.{lvl}_means_cell_log_counts_pM.tsv") for id in config["scrnas"] for lvl in levels_per_dataset[id]]

checkpoint fix_index:
    input: 
        h5ad = os.path.join(config["input"],"{id}/{id}.h5ad")
    output:
        temp("input_temp/{id}_temp_{lvl}.h5ad")  
    params:
        scriptdir = config["scriptdir"]
    resources:
        mem_mb = mem_from_input(multiplier=2, minimum=4000)
    shell:
        """
        python {params.scriptdir}/fix_index.py --level {wildcards.lvl} --h5ad {input.h5ad} --output {output}
        """    

    
rule run_cellex:
    input:
        h5ad = os.path.join("input_temp","{id}_temp_{lvl}.h5ad")
    output:
        os.path.join(config["metrics"],"cellex/{id}/{lvl}/{id}.{lvl}.esmu_fmt.tsv")
    params:
        outdir = os.path.join(config["metrics"],"cellex/{id}/{lvl}"),
        scriptdir = config["scriptdir"]
    resources:
        mem_mb=mem_from_input(multiplier=8, minimum=4000)
    conda:
       "cellex"
    shell:
        """
        python {params.scriptdir}/process_cellex.py --h5ad {input.h5ad} --level {wildcards.lvl} --outdir {params.outdir} --cellex_out_prefix {wildcards.id} 
        """

rule run_pre_ewce:
    input:
        h5ad = os.path.join("input_temp","{id}_temp_{lvl}.h5ad")
    output:
        temp(os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}.{lvl}.filtered_protein_coding_aov.h5ad"))
    params:
        scriptdir = config["scriptdir"]
    resources:
        mem_mb=mem_from_input(multiplier=8, minimum=4000)
    shell:
        """
        python {params.scriptdir}/filter_h5ad_protein_coding_anova_2.py --h5ad {input.h5ad} --level {wildcards.lvl} --output {output}
        """

rule run_ewce:
    input:
        h5ad = os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}.{lvl}.filtered_protein_coding_aov.h5ad")
    output:
        os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}_ensemble_{lvl}.rda")
    params:
        outdir = os.path.join(config["metrics"],"ewce/{id}/{lvl}"),
        scriptdir = config["scriptdir"],
        file_prefix = "{id}_ensemble"
    resources:
        mem_mb=mem_from_input(multiplier=8, minimum=4000)
    conda:
       "R_with_python"
    shell:
        """
        Rscript {params.scriptdir}/process_scrna_ewce.R --h5ad {input.h5ad} --level {wildcards.lvl} --outdir {params.outdir} --file_prefix {params.file_prefix} 
        """


rule post_ewce:
   input: 
      os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}_ensemble_{lvl}.rda")
   output:
      os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}.{lvl}_spec_ewce_linear.tsv"),
      os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}.{lvl}_mean_ewce_linear.tsv")
   params:
      scriptdir = config["scriptdir"],
      outdir = os.path.join(config["metrics"],"ewce/{id}/{lvl}")
   conda:
      "R_with_python"
   shell:
      """
      Rscript {params.scriptdir}/process_ewce.R --infile {input} --outdir {params.outdir} --level {wildcards.lvl} --dataset {wildcards.id} 
      """   

rule create_spec_top10:
    input:
        os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}.{lvl}_spec_ewce_linear.tsv")
    output:
        os.path.join(config["metrics"],"ewce/{id}/{lvl}/{id}.{lvl}_spec_ewce_top10.txt")
    params:
      	scriptdir = config["scriptdir"],
      	outdir_magma = os.path.join(config["metrics"],"ewce/{id}/{lvl}"),
    conda:
       "R_with_python"
    shell:
       """
       mkdir -p {params.outdir_magma} 
       python {params.scriptdir}/post_process_spec_ewce.py --spec_linear {input} --output {output} --outdir_magma {params.outdir_magma}
       """

rule compute_cepo:
    input: 
        h5ad = os.path.join("input_temp","{id}_temp_{lvl}.h5ad")
    output:
        os.path.join(config["metrics"],"cepo/{id}/{lvl}/{id}.{lvl}_cepo_filtered.tsv")
    params:
        outdir = os.path.join(config["metrics"],"cepo/{id}/{lvl}"),
        scriptdir = config["scriptdir"],
    resources:
        mem_mb=mem_from_input(multiplier=12, minimum=4000)
    conda:
       "R_CEPO"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.scriptdir}/compute_cepo.R --infile {input.h5ad} --dataset {wildcards.id} --level {wildcards.lvl} --outdir {params.outdir} 
        """

rule create_cepo_top10:
    input:
        os.path.join(config["metrics"],"cepo/{id}/{lvl}/{id}.{lvl}_cepo_filtered.tsv")
    output:
        os.path.join(config["metrics"],"cepo/{id}/{lvl}/{id}.{lvl}_cepo_top10.txt")
    params:
        scriptdir = config["scriptdir"],
        outdir_magma = os.path.join(config["metrics"],"cepo/{id}/{lvl}"),
    conda:
       "R_with_python"
    shell:
       """
       mkdir -p {params.outdir_magma}
       python {params.scriptdir}/post_process_cepo.py --spec_linear {input} --output {output} --outdir_magma {params.outdir_magma} 
       """

rule compute_sumstats_FUMA_cell_typing:
    input:
        h5ad = os.path.join("input_temp","{id}_temp_{lvl}.h5ad")
    output:
        os.path.join(config["metrics"],"fuma/{id}/{lvl}/{id}.{lvl}_means_cell_log_counts_pM.tsv")
    params:
        outdir = os.path.join(config["metrics"], "fuma/{id}/{lvl}"),
	script = os.path.join(config["scriptdir"],"compute_sumstat_magma.py")
    resources:
        mem_mb=mem_from_input(multiplier=4, minimum=4000)
    shell:
        """
        python {params.script} --h5ad {input} --id {wildcards.id} --level {wildcards.lvl} --out {params.outdir};
        """

