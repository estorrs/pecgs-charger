arguments:
- position: 0
  prefix: --out-dir
  valueFrom: output
- position: 0
  prefix: --format-vcf-script
  valueFrom: /pecgs-charger/src/scripts/format_vcf_for_CharGer.py
- position: 0
  prefix: --post-charger-script
  valueFrom: /pecgs-charger/src/scripts/post_CharGer.py
- position: 0
  prefix: --filter-charger-script
  valueFrom: /pecgs-charger/src/scripts/filter_CharGer.py
baseCommand:
- python
- /pecgs-charger/src/charger.py
class: CommandLineTool
cwlVersion: v1.0
id: charger
inputs:
- id: vcf
  inputBinding:
    position: '1'
  type: File
- id: inheritance_gene_list
  inputBinding:
    position: '0'
    prefix: --inheritance-gene-list
  type: File
- id: pp2_gene_list
  inputBinding:
    position: '0'
    prefix: --pp2-gene-list
  type: File
- id: pathogenic_variants
  inputBinding:
    position: '0'
    prefix: --pathogenic-variants
  type: File
- id: hotspot3d_clusters
  inputBinding:
    position: '0'
    prefix: --hotspot3d-clusters
  type: File
- id: clinvar_alleles
  inputBinding:
    position: '0'
    prefix: --clinvar-alleles
  type: File
- default: sample
  id: sample
  inputBinding:
    position: '0'
    prefix: --sample
  type: string?
- default: '0.0005'
  id: rare_threshold
  inputBinding:
    position: '0'
    prefix: --rare-threshold
  type: string?
- default: /miniconda/envs/charger/bin:$PATH
  id: environ_PATH
  type: string?
label: charger
outputs:
- id: filtered_tsv
  outputBinding:
    glob: output/4.filter_charger/$(inputs.sample).charged2vcf.filtered.tsv
  type: File
- id: rare_threshold_filtered_tsv
  outputBinding:
    glob: output/4.filter_charger/$(inputs.sample).charged2vcf.filtered.af$(inputs.rare_threshold).tsv
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/pecgs-charger:0.0.1
- class: ResourceRequirement
  coresMin: 1
  ramMin: 29000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)
