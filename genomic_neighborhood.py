#! /usr/bin/env python
from __future__ import print_function
import os, sys, argparse, warnings, gzip, bz2, zipfile, time
from io import TextIOWrapper
from collections import defaultdict, OrderedDict
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from soothsayer_utils import infer_compression, format_path, format_duration, format_header, boolean, assert_acceptable_arguments, pv
from dna_features_viewer import GraphicFeature, GraphicRecord


pd.options.display.max_colwidth = 100
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2020.03.03"
__cite__ = "{}\nDNA Features Viewer, a sequence annotations formatting and plotting library for Python\nValentin Zulkower et al. bioRxiv 2020.01.09.900589; doi: https://doi.org/10.1101/2020.01.09.900589\n\nTranscriptomic study of substrate-specific transport mechanisms for iron and carbon in the marine copiotroph Alteromonas macleodii\nLauren E. Manck et al. 2020. mSystems".format(format_header("Citations"))

# Get file object
def get_file_object(path, mode="infer", compression="infer", safe_mode="infer"):
    # Adapted from Soothsayer
    # https://github.com/jolespin/soothsayer

    # Init
    f = None
    file_will_be_written = False

    # Paths
    path = format_path(path)
    path_exists = os.path.exists(path)

    if compression == "infer":
        compression = infer_compression(path)


    # Inferring mode
    if mode == "infer": # Create new function for this? infer_filemode?
        if path_exists:
            mode = "read"
        else:
            mode = "write"
    assert mode != "infer", "The mode should be inferred by this point.  Please specify mode manually."
    assert compression != "infer", "The compression should be inferred by this point.  Please specify compression manually."

    # Generic read write
    if mode in ["write", "read"]:
        if mode == "write":
            mode = "w"
        if mode == "read":
            mode = "r"
        if compression in ["gzip", "bz2"]:
            mode = mode + "b"


    # Will a file be written?
    if "w" in mode:
        file_will_be_written = True

    # Ensure zip is not being written
    if file_will_be_written:
        assert compression != "zip", "Currently cannot handle writing zip files.  Please use gzip, bz2, or None."
        # Future do this:
        # https://docs.python.org/3/library/zipfile.html#zipfile.ZipFile.open

    # Safe mode
    if safe_mode == "infer":
        if file_will_be_written:
            safe_mode = True
        else:
            safe_mode = False
    assert safe_mode in {True,False}, "Please choose either True or False for `safe_mode`"
    if safe_mode:
        if all([file_will_be_written, path_exists]):
            raise Exception("Safe Mode: Please explicitly provide a writeable mode ['w', 'wb', or 'write'] because `{}` already exists and will be rewritten.".format(path))

    # GZIP compression
    if compression == "gzip":
        f = gzip.open(path, mode)
    # BZ2 compression
    if compression == "bz2":
        f = bz2.open(path, mode)
    if compression == "zip":
        filename, ext = os.path.splitext(os.path.split(path)[-1])
        f = zipfile.ZipFile(path,mode).open(filename)
    # No compression
    if f is None:
        return open(path, mode)

    # Reading and writing compressed files
    else:
        return TextIOWrapper(f, encoding="utf-8")

# GTF/GFF
def _read_gtf_gff_base(path, compression, record_type):
    # Adapted from Soothsayer
    # https://github.com/jolespin/soothsayer

    # Read the gff3 file
    with get_file_object(path, mode="read", compression=compression, safe_mode=False) as f:
        iterable_lines = f.readlines()
        data= list()
        if record_type is None:
            for line in iterable_lines:
                if not line.startswith("#"):
                    line = line.strip("\n")
                    if bool(line):
                        base_fields = line.split("\t")
                        data.append(base_fields)
        else:
            for line in iterable_lines:
                if not line.startswith("#"):
                    if "{}_id".format(record_type) in line:
                        line = line.strip("\n")
                        base_fields = line.split("\t")
                        data.append(base_fields)
    # Generate table
    df_base = pd.DataFrame(data)
    df_base.columns = ["seq_record", "source", "seq_type", "pos_start", "pos_end", ".1", "sense", ".2", "data_fields"]
    return df_base

def read_gff3(path, compression="infer", record_type=None, reset_index=False, name=True):
    def f(x):
        fields = x.split(";")
        data = dict()
        for item in fields:
            k, v = item.split("=")
            data[k] = v
        return data
    # Adapted from Soothsayer
    # https://github.com/jolespin/soothsayer

    path = format_path(path)

    accepted_recordtypes = {"exon", "gene", "transcript", "protein", None}
    assert record_type in accepted_recordtypes, "Unrecognized record_type.  Please choose from the following: {}".format(accepted_recordtypes)

    # Read the gff3 file
    df_base = _read_gtf_gff_base(path, compression, record_type)

    try:
        df_fields = pd.DataFrame(df_base["data_fields"].map(f).to_dict()).T
        df_gff3 = pd.concat([df_base[["seq_record", "source", "seq_type", "pos_start", "pos_end", ".1", "sense", ".2"]], df_fields], axis=1)
         # Contains identifier
        if "ID" in df_gff3.columns:
            df_gff3["seq_id"] = df_gff3["ID"].map(lambda x: "-".join(x.split("-")[1:]) if not pd.isnull(x) else x)
        if reset_index:
            df_gff3 = df_gff3.reset_index(drop=True)
    except IndexError:
        warnings.warn("Could not successfully parse file: {}".format(path))
        df_gff3 = df_base

    # Path index
    if name is not None:
        if name == True:
            name = path
        df_gff3.index.name = name

    return df_gff3

#GTF
def read_gtf(path, compression="infer", record_type=None, reset_index=False, name=True):
    # Adapted from Soothsayer
    # https://github.com/jolespin/soothsayer

    path = format_path(path)

    accepted_recordtypes = {"exon", "gene", "transcript", "protein", None}
    assert record_type in accepted_recordtypes, "Unrecognized record_type.  Please choose from the following: {}".format(accepted_recordtypes)

    # Read the gtf file
    df_base = _read_gtf_gff_base(path, compression, record_type)

    # Splitting fields
    iterable_fields =  df_base.iloc[:,-1].iteritems()
    # Parse fields
    dataframe_build = dict()
    for i, gtf_data in iterable_fields:
        gtf_data = gtf_data.replace("''","").replace('"',"")
        fields = filter(bool, map(lambda field:field.strip(), gtf_data.split(";")))
        fields = map(lambda field: field.split(" "), fields)
        dataframe_build[i] = dict(fields)
    df_fields = pd.DataFrame(dataframe_build).T
    df_gtf = pd.concat([df_base.iloc[:,:7], df_fields], axis=1)

    # Reset index
    if reset_index:
        df_gtf = df_gtf.reset_index(drop=True)

    # Path index
    if name is not None:
        if name == True:
            name = path
        df_gtf.index.name = name
    return df_gtf

# Read Dataframe
def read_dataframe(path, sep="infer", index_col=0, header=0, compression="infer", pickled="infer", engine="c",  excel="infer", sheet_name=None,  **args):
    # Adapted from Soothsayer
    # https://github.com/jolespin/soothsayer

    start_time = time.time()
    path = format_path(path, str)
    dir , ext = os.path.splitext(path)
    ext = ext.lower()

    if excel == "infer":
        if ext in {".xlsx", ".xls"}:
            excel = True
        else:
            excel = False

    if excel:
        assert sheet_name is not None, "Please specify `sheet_name` if using Excel files.  Alternatively, you can give it a text table as tab or comma-separated."
        df = pd.read_excel(path, sheet_name=sheet_name, index_col=index_col, header=header, **args)

    else:
        # Seperator
        if any(list(map(lambda x:ext.endswith(x),[".csv", "csv.gz", "csv.zip"]))):
            sep = ","
        else:
            sep = "\t"

        # Serialization
        if pickled == "infer":
            if ext in {".pkl", ".pgz", ".pbz2"}:
                pickled = True
            else:
                pickled = False
        # Compression
        if compression == "infer":
            if pickled:
                if ext == ".pkl":
                    compression = None
                if ext == ".pgz":
                    compression = "gzip"
                if ext == ".pbz2":
                    compression = "bz2"
            else:
                if ext == ".gz":
                    compression = "gzip"
                if ext == ".bz2":
                    compression = "bz2"
                if ext == ".zip":
                    compression = "zip"

        if pickled:
            df = pd.read_pickle(path, compression=compression)
        else:
            df = pd.read_csv(path, sep=sep, index_col=index_col, header=header, compression=compression, engine=engine, **args)



    return df

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = sys.argv[0].split("/")[-1]
    description = """
    Running: {} v{} via Python v{} | {}\n
    Versatile command-line tool designed to create publication quality genomic neighborhood plots built on top of DnaFeaturesViewer

    [Example: --annotation_table ./Data/annotation_table.xlsx]
    # Note: Columns are not case sensitive
    Group	Locus_Tag	Color
    9	MASE_00925	#FFFFFF
    9	MASE_00930	#FFFFFF
    9	MASE_00935	#FF0000
    9	MASE_00940	#FFFFFF

    [Example: --feature_table ./Data/CP003841.1.gff3]
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

    [Example commands:]
    ./genomic_neighborhood.py -f ./Data/CP003841.1.gff3 -a ./Data/annotation_table.xlsx -o genomic_neighborhood_output --sheet_name "Fe Responsive"
    """.format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -f <feature_table> -a <annotation_table> -o <output_directory>".format(__program__)
    epilog = "Copyright 2020 Josh L. Espinoza (jespinoz@jcvi.org) [BSD-3 License]"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Features Table
    parser_features = parser.add_argument_group('Feature table arguments')
    parser_features.add_argument("-f", "--feature_table", type=str, help = "path/to/feature_table.[gff3,gtf][.gz,.bz2,.zip] (e.g. feature_table.gff3[.gz]) {gff3,gtf}")
    parser_features.add_argument( "--field", type=str, default="locus_tag", help = "Query feature.  Note: locus_tag is the only feature accepted with current version. [Default: locus_tag]")
    parser_features.add_argument( "--feature_format", type=str, default="infer", help = "Feature format. [Default: infer] {gff3,gtf}")


    # Annotation Table
    parser_annotations = parser.add_argument_group('Annotation table arguments')
    parser_annotations.add_argument("-a", "--annotation_table", type=str, help = "path/to/annotation_table.[ext][.compression] (e.g. annotation_table.tsv[.gz]).  Usable columns are [group, <--field>, color ].  Required column is [<--field>] (e.g. [locus_tag]) {tsv,csv,xlsx}")
    parser_annotations.add_argument("--excel", type=str, default="infer", help = "Input table is excel format {true, false, infer} [Default: infer]")
    parser_annotations.add_argument("--sep", type=str, default="\t", help = "Separator for input table [Default: '\\t']")
    parser_annotations.add_argument("--sheet_name", type=str, help = "Sheetname if using excel")

    # Image
    parser_images = parser.add_argument_group('Image arguments')
    parser_images.add_argument( "-o","--output_directory", type=str, help = "Output direcotyr [Default: {}]".format(os.getcwd()))
    parser_images.add_argument( "--feature_color", type=str, default="gray", help = "Feature color as a hexcode (#929591) or named color (gray).  Note this is the default color and will be overrided if a 'Color' column is provided for `--anotation_table`\nReference: https://matplotlib.org/3.1.0/gallery/color/named_colors.html")
    parser_images.add_argument( "--feature_opacity", type=str, default=0.85, help = "Feature color opacity. [Default: 0.85] [0.0,..,1.0]")
    parser_images.add_argument( "--image_format", type=str, default="svg", help = "Image format [Default: svg] {svg,png,pdf}")
    parser_images.add_argument( "--show_sequence_record", type=str, default="f", help = "Add sequence record title. [Default: false]")
    parser_images.add_argument("--figure_width", type=float, default=20.0, help="Width of figures [Default: 20]")
    parser_images.add_argument("--draw_reference_line", type=str, default="t", help="Draw reference line for features [Default: true]")

    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))
    parser_utility.add_argument("--citation", action='store_true', help="If you use this software, please cite the following sources:\n{}".format(__cite__))

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    if opts.citation:
        print(__cite__,file=sys.stderr)
        sys.exit(0)


    print(format_header(__program__), file=sys.stdout)

    # Read in annotations
    if opts.sep in {"comma", "csv"}:
        opts.sep = ","
    if opts.sep in {"tab", "tsv", "t"}:
        opts.sep = "\t"
    if opts.sep in {"\s", "space"}:
        opts.sep = " "

    opts.field = opts.field.lower()
    assert opts.field == "locus_tag", "Currently only `locus_tag` is supported for --field"

    df_annotations = read_dataframe(
        path = opts.annotation_table,
        sep = opts.sep,
        excel = opts.excel,
        sheet_name = opts.sheet_name,
        index_col=None,
    )
    df_annotations.columns = df_annotations.columns.map(lambda x:x.strip().lower())
    print("Reading annotation table: {} | {}".format(opts.annotation_table, df_annotations.shape), file=sys.stderr)

    assert opts.field in df_annotations.columns, "--field ({}) not --feature_table columns".format(opts.field)
    if "group" not in df_annotations.columns:
        df_annotations["group"] = None
    if "color" not in df_annotations.columns:
        df_annotations["color"] = opts.feature_color

    # Read in feature table
    opts.feature_format = opts.feature_format.lower()
    if opts.feature_format == "infer":
        if any(x in opts.feature_table for x in {".gff", ".gff3"}):
            opts.feature_format = "gff3"
        if ".gtf" in opts.feature_table:
            opts.feature_format = "gtf"
    assert opts.feature_format != "infer", "Could not infer `feature_format`.  Please specify either {gff3, gtf}"
    assert_acceptable_arguments(opts.feature_format, {"gff", "gff3", "gtf"})
    if opts.feature_format in {"gff", "gff3"}:
        df_features = read_gff3(opts.feature_table)
    if opts.feature_format in {"gtf"}:
        df_features = read_gtf(opts.feature_table)
    print("Reading features table: {} | {}".format(opts.feature_table, df_features.shape), file=sys.stderr)

    df = df_features.query("seq_type == 'region'")
    d_seq_length = dict(zip(df["seq_record"], df["pos_end"].astype(int)))
    df_features = df_features.loc[df_features["seq_type"][lambda x:x == "CDS"].index,:].reset_index()

    # Output directory
    os.makedirs(opts.output_directory, exist_ok=True)

    # Genomic neighborhoods
    d_id_group = dict(zip(df_annotations[opts.field], df_annotations["group"]))
    d_id_color = dict(zip(df_annotations[opts.field], df_annotations["color"]))

    d_group_features = defaultdict(list)
    d_group_positions = defaultdict(list)
    positions = list()
    n_peripheral = 10

    for i, data in df_features.iterrows():
        seq_record = data["seq_record"]
        sequence_length = d_seq_length[seq_record]

        start = int(data["pos_start"])
        end = int(data["pos_end"])
        strand = {"+":+1, "-":-1}[data["sense"]]
        label = data[opts.field]
        if label in d_id_color:
            group = d_id_group[label]
            feature = GraphicFeature(start=start, end=end, strand=strand, label=label, color=(*to_rgb(d_id_color[label]), opts.feature_opacity))
            d_group_features[group].append(feature)
            d_group_positions[group] += [start, end]

    lengths = dict()
    limits = dict()
    for group, positions in d_group_positions.items():
        lengths[group] = max(positions) - min(positions)
        limits[group] = (min(positions), max(positions))

    max_length = max(lengths.values())
    pad = max_length//2 + 10
    for group in d_group_features:
        features = d_group_features[group]
        positions = d_group_positions[group]
        midpoint = np.mean(limits[group])
        record = GraphicRecord(sequence_length=sequence_length, features=features).crop((max(midpoint-pad,0), midpoint+pad))

        ax, results = record.plot(figure_width=opts.figure_width, draw_line=boolean(opts.draw_reference_line))
        if boolean(opts.show_sequence_record):
            ax.set_title("{} [{}..{}]".format(seq_record,limits[group][0], limits[group[1]]), fontsize=15, fontweight="bold")
        if group is not None:
            output_filepath = os.path.join(opts.output_directory, "{}__{}.{}".format(seq_record, group, opts.image_format))
        else:
            output_filepath = os.path.join(opts.output_directory, "{}.{}".format(seq_record,opts.image_format))

        ax.figure.savefig(output_filepath, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
