## genomic_neighborhood.py
![image](https://i.imgur.com/T6I0JTa.png)

### Description:
Versatile command-line tool designed to create publication quality genomic neighborhood plots built on top of [DnaFeaturesViewer](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/). This implementation streamlines usage and can handle excel, various delimiters, gff3, gtf, and various compressed files.  

### Requirements:
```
# Python packages:
pandas
numpy
matplotlib
dna_features_viewer
soothsayer_utils

# Optional:
conda
```

### Installation:
#### Download latest version
```
git clone https://github.com/jolespin/genomic_neighborhood
cd genomic_neighborhood
```

#### Conda
`conda env create --name genomic_neighborhood_env --file environment.yaml`

#### PyPi (PIP)
`pip install numpy pandas matplotlib soothsayer_utils genopype dna_feature_viewer`
________

### Arguments:
```
usage: genomic_neighborhood.py -f <feature_table> -a <annotation_table> -o <output_directory>

    Running: genomic_neighborhood.py v2020.03.03 via Python v3.7.3 | /Users/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Feature table arguments:
  -f FEATURE_TABLE, --feature_table FEATURE_TABLE
                        path/to/feature_table.[gff3,gtf][.gz,.bz2,.zip] (e.g. feature_table.gff3[.gz]) {gff3,gtf}
  --field FIELD         Query feature.  Note: locus_tag is the only feature accepted with current version. [Default: locus_tag]
  --feature_format FEATURE_FORMAT
                        Feature format. [Default: infer] {gff3,gtf}

Annotation table arguments:
  -a ANNOTATION_TABLE, --annotation_table ANNOTATION_TABLE
                        path/to/annotation_table.[ext][.compression] (e.g. annotation_table.tsv[.gz]).  Usable columns are [group, <--field>, color ].  Required column is [<--field>] (e.g. [locus_tag]) {tsv,csv,xlsx}
  --excel EXCEL         Input table is excel format {true, false, infer} [Default: infer]
  --sep SEP             Separator for input table [Default: '\t']
  --sheet_name SHEET_NAME
                        Sheetname if using excel

Image arguments:
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output direcotyr [Default: /Users/jespinoz/Google/Manuscripts/In_progress/NASA__alteromonas/genomic_neighborhood]
  --feature_color FEATURE_COLOR
                        Feature color as a hexcode (#929591) or named color (gray).  Note this is the default color and will be overrided if a 'Color' column is provided for `--anotation_table`
                        Reference: https://matplotlib.org/3.1.0/gallery/color/named_colors.html
  --feature_opacity FEATURE_OPACITY
                        Feature color opacity. [Default: 0.85] [0.0,..,1.0]
  --image_format IMAGE_FORMAT
                        Image format [Default: svg] {svg,png,pdf}
  --show_sequence_record SHOW_SEQUENCE_RECORD
                        Add sequence record title. [Default: false]
  --figure_width FIGURE_WIDTH
                        Width of figures [Default: 20]
  --draw_reference_line DRAW_REFERENCE_LINE
                        Draw reference line for features [Default: true]

Utility arguments:
  -v, --version         show program's version number and exit
  --citation            If you use this software, please cite the following sources:
                        =========
                        Citations
                        =========
                        DNA Features Viewer, a sequence annotations formatting and plotting library for Python
                        Valentin Zulkower et al. bioRxiv 2020.01.09.900589; doi: https://doi.org/10.1101/2020.01.09.900589

                        Transcriptomic study of substrate-specific transport mechanisms for iron and carbon in the marine copiotroph Alteromonas macleodii
                        Lauren E. Manck et al. 2020. mSystems

Copyright 2020 Josh L. Espinoza (jespinoz@jcvi.org) [BSD-3 License]
```
____________

### Usage tutorial:
#### `--annotation_table ./Data/annotation_table.xlsx`
Note: Columns are not case sensitive

```
Group	Locus_Tag	Color
9	MASE_00925	#FFFFFF
9	MASE_00930	#FFFFFF
9	MASE_00935	#FF0000
9	MASE_00940	#FFFFFF
```

#### `--feature_table ./Data/CP003841.1.gff3`
```
##sequence-region CP003841.1 1 4653851
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=529120
CP003841.1	Genbank	region	1	4653851	.	+	.	ID=CP003841.1:1..4653851;Dbxref=taxon:529120;Is_circular=true;Name=ANONYMOUS;country=Pacific Ocean: near Hawaii;gbkey=Src;genome=chromosome;isolation-source=seawater surface;mol_type=genomic DNA;old-lineage=Bacteria%3B Proteobacteria%3B Gammaproteobacteria%3B Alteromonadales%3B Alteromonadaceae%3B Alteromonas;strain=ATCC 27126
CP003841.1	Genbank	gene	473	2065	.	+	.	ID=gene-MASE_00005;Name=MASE_00005;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MASE_00005
CP003841.1	Genbank	CDS	473	2065	.	+	0	ID=cds-AFS35578.1;Parent=gene-MASE_00005;Dbxref=NCBI_GP:AFS35578.1;Name=AFS35578.1;Note=COG0593 ATPase involved in DNA replication initiation;gbkey=CDS;locus_tag=MASE_00005;product=chromosomal replication initiator protein dnaA;protein_id=AFS35578.1;transl_table=11
CP003841.1	Genbank	gene	2098	3198	.	+	.	ID=gene-MASE_00010;Name=MASE_00010;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MASE_00010
CP003841.1	Genbank	CDS	2098	3198	.	+	0	ID=cds-AFS35579.1;Parent=gene-MASE_00010;Dbxref=NCBI_GP:AFS35579.1;Name=AFS35579.1;Note=COG0592 DNA polymerase sliding clamp subunit (PCNA homolog);gbkey=CDS;locus_tag=MASE_00010;product=DNA polymerase III subunit beta;protein_id=AFS35579.1;transl_table=11
CP003841.1	Genbank	gene	3324	4412	.	+	.	ID=gene-MASE_00015;Name=MASE_00015;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MASE_00015
CP003841.1	Genbank	CDS	3324	4412	.	+	0	ID=cds-AFS35580.1;Parent=gene-MASE_00015;Dbxref=NCBI_GP:AFS35580.1;Name=AFS35580.1;Note=COG1195 Recombinational DNA repair ATPase (RecF pathway);gbkey=CDS;locus_tag=MASE_00015;product=recombinational DNA repair ATPase;protein_id=AFS35580.1;transl_table=11
CP003841.1	Genbank	gene	4421	6841	.	+	.	ID=gene-MASE_00020;Name=MASE_00020;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MASE_00020
```

#### Commands:
```
# Fe Responsive
./genomic_neighborhood.py -f ./Data/CP003841.1.gff3 -a ./Data/annotation_table.xlsx -o fe_responsive_output --sheet_name "Fe Responsive"

# C Responsive
./genomic_neighborhood.py -f ./Data/CP003841.1.gff3 -a ./Data/annotation_table.xlsx -o c_responsive_output --sheet_name "C Responsive"
```
_________
### Citation:
```
DNA Features Viewer, a sequence annotations formatting and plotting library for Python
Valentin Zulkower et al. bioRxiv 2020.01.09.900589; doi: https://doi.org/10.1101/2020.01.09.900589

Transcriptomic study of substrate-specific transport mechanisms for iron and carbon in the marine copiotroph Alteromonas macleodii
Lauren E. Manck et al. 2020. mSystems
                        ```