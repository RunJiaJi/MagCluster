#args.py handles arguments

def get_magcluster_parser():
    import argparse
    
    parser = argparse.ArgumentParser(
             prog="magcluster", 
             description='Magnetosome Gene Clusters Analyzer', 
             usage='%(prog)s [options]', 
             formatter_class=argparse.RawTextHelpFormatter, 
             epilog=
            "General usage\n-------------\n"
            "Genome annotation:\n"
            "  $ magcluster prokka genome.fasta\n\n"
            "MGCs screening:\n"
            "  $ magcluster mgc_screen genbank.gbk\n\n"
            "MGCs alignment and mapping:\n"
            "  $ magcluster clinker genbank_mgc.gbk\n\n"
            
            "Runjia, 2021"
            )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.9', help='show magcluster version number and exit')
    # parser.add_argument('-ascm', action='store_true', help='directly analyse from genome file to genecluster mapping.')
    #构建子命令
    subparsers = parser.add_subparsers(title='Options', dest="subparser_name")
    #构建maga子命令
    parser_maga = subparsers.add_parser('prokka', help='Genome annotation with Prokka')
    parser_maga.add_argument('fafile', type=str, help='Genome files need to be annotated', nargs='+')

    General = parser_maga.add_argument_group('General')
    General.add_argument('--version', help='Print version and exit', action="store_true")
    General.add_argument('--docs', help='Show full manual/documentation', action="store_true")
    General.add_argument('--citation', help='Print citation for referencing Prokka', action="store_true")
    General.add_argument('--quiet', help='No screen output (default OFF)', action="store_true")
    General.add_argument('--debug', help='Debug mode: keep all temporary files (default OFF)', action="store_true")

    Setup = parser_maga.add_argument_group('Setup')
    Setup.add_argument('--listdb', help='List all configured databases', action="store_true")
    Setup.add_argument('--setupdb', help='Index all installed databases', action="store_true")
    Setup.add_argument('--cleandb', help='Remove all database indices', action="store_true")
    Setup.add_argument('--depends', help='List all software dependencies', action="store_true")

    Outputs = parser_maga.add_argument_group('Outputs')
    Outputs.add_argument('--outdir', type=str, help="Output folder [auto] (default 'XXX_annotation')")
    Outputs.add_argument('--prefix', type=str, help='Filename output prefix [auto] (default "XXX_")')
    Outputs.add_argument('--force', help='Force overwriting existing output folder (default OFF)', action="store_true")
    Outputs.add_argument('--addgenes', help="Add 'gene' features for each 'CDS' feature (default OFF)", action="store_true")
    Outputs.add_argument('--addmrna', help="Add 'mRNA' features for each 'CDS' feature (default OFF)", action="store_true")
    Outputs.add_argument('--locustag', type=str, help="Locus tag prefix (default 'PROKKA')")
    Outputs.add_argument('--increment', type=int, help="Locus tag counter increment (default '1')")
    Outputs.add_argument('--gffver', type=int, help="GFF version (default '3')")
    Outputs.add_argument('--compliant', help='Force Genbank/ENA/DDJB compliance: --genes --mincontiglen 200 --centre XXX (default OFF)', action="store_true")

    somthing = parser_maga.add_argument_group('XXX (default OFF)')
    somthing.add_argument('--centre', type=str, help="Sequencing centre ID. (default '')")
    somthing.add_argument('--accver', type=int, help="Version to put in Genbank file (default '1')")

    Organism_details = parser_maga.add_argument_group('Organism details')
    Organism_details.add_argument('--genus', type=str, help="Genus name (default 'Genus')")
    Organism_details.add_argument('--species', type=str, help="Species name (default 'species')")
    Organism_details.add_argument('--strain', type=str, help="Strain name (default 'strain')")
    Organism_details.add_argument('--plasmid', type=str, help="Plasmid name or identifier (default '')")

    Annotations = parser_maga.add_argument_group('Annotations')
    Annotations.add_argument('--kingdom', type=str, help="Annotation mode: Archaea|Bacteria|Viruses (default 'Bacteria')")
    Annotations.add_argument('--gcode', type=int, help="Genetic code / Translation table (set if --kingdom is set) (default '0')")
    Annotations.add_argument('--gram', type=str, help="Gram: -/neg +/pos (default '')")
    Annotations.add_argument('--usegenus', help='Use genus-specific BLAST databases (needs --genus) (default OFF)', action="store_true")
    Annotations.add_argument('--proteins', type=str, help='Fasta file of trusted proteins to first annotate from (default "Magnetosome_protein_data.fasta")')
    Annotations.add_argument('--hmms', type=str, help="Trusted HMM to first annotate from (default '')")
    Annotations.add_argument('--metagenome', help='Improve gene predictions for highly fragmented genomes (default OFF)', action="store_true")
    Annotations.add_argument('--rawproduct', help='Do not clean up /product annotation (default OFF)', action="store_true")
    Annotations.add_argument('--cdsrnaolap', help="Allow [tr]RNA to overlap CDS (default OFF)", action="store_true")

    Computation = parser_maga.add_argument_group('Computation')
    Computation.add_argument('--cpus', type=int, help="Number of CPUs to use [0=all] (default '8')")
    Computation.add_argument('--fast', help='Fast mode - skip CDS /product searching (default OFF)', action="store_true")
    Computation.add_argument('--noanno', help='For CDS just set /product="unannotated protein" (default OFF)', action="store_true")
    Computation.add_argument('--mincontiglen', type=int, help="Minimum contig size [NCBI needs 200] (default '1')")
    Computation.add_argument('--evalue', type=float, help="Similarity e-value cut-off (default '1e-06')")
    Computation.add_argument('--rfam', help="Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default '0')", action="store_true")
    Computation.add_argument('--norrna', help="Don't run rRNA search (default OFF)", action="store_true")
    Computation.add_argument('--notrna', help="Don't run tRNA search (default OFF)", action="store_true")
    Computation.add_argument('--rnammer', help="Prefer RNAmmer over Barrnap for rRNA prediction (default OFF)", action="store_true")
    
    #构建magsc子命令
    parser_magsc = subparsers.add_parser('mgc_screen', help='Magnetosome gene cluster screening')
    parser_magsc.add_argument('gbkfile', type=str, help='.gbk/.gbf files to be analyzed. Multiple files or files-containing folder is acceptable.', nargs='+')
    parser_magsc.add_argument('-l', '--contiglength', type=int, default=2000, help="The minimum length of contigs to be considered (default '2,000 bp')")
    parser_magsc.add_argument('-w', '--windowsize', type=int, default=10000, help="The length of MGCs screening window (default '10,000 bp')")
    parser_magsc.add_argument('-th','--threshold', type=int, default=3, help="The minimum number of magnetosome genes in a given contig and a given length of screening window (default '3')")
    parser_magsc.add_argument('-o','--outdir', type=str, help="Output folder (default 'mgc_screen')")
    
    #构建magm子命令
    parser_magm = subparsers.add_parser('clinker', help='Magnetosome gene cluster mapping with Clinker')
    inputs = parser_magm.add_argument_group("Input options")
    inputs.add_argument('gbkfiles', help="Gene cluster GenBank files", nargs="*")
    inputs.add_argument("-r", "--ranges", 
           help="Scaffold extraction ranges. If a range is specified, only features within the range will be extracted from the scaffold. Ranges should be formatted like: scaffold:start-end"
           " (e.g. scaffold_1:15000-40000)", nargs="+",
    )

    alignment = parser_magm.add_argument_group("Alignment options")
    alignment.add_argument(
        "-na",
        "--no_align",
        help="Do not align clusters",
        action="store_true",
    )
    alignment.add_argument(
        "-i",
        "--identity",
        help="Minimum alignment sequence identity [default: 0.3]",
        type=float,
        default=0.3
    )
    alignment.add_argument(
        "-j",
        "--jobs",
        help="Number of alignments to run in parallel (0 to use the number of CPUs) [default: 0]",
        type=int,
        default=0,
    )

    output = parser_magm.add_argument_group("Output options")
    output.add_argument("-s", "--session", help="Path to clinker session")
    output.add_argument("-ji", "--json_indent", type=int, help="Number of spaces to indent JSON [default: none]")
    output.add_argument("-f", "--force", help="Overwrite previous output file", action="store_true")
    output.add_argument("-o", "--output", help="Save alignments to file")
    output.add_argument(
        "-p",
        "--plot",
        nargs="?",
        const=True,
        default=False,
        help="Plot cluster alignments using clustermap.js. If a path is given,"
        " clinker will generate a portable HTML file at that path. Otherwise,"
        " the plot will be served dynamically using Python's HTTP server."
    )
    output.add_argument("-dl", "--delimiter", help="Character to delimit output by [default: human readable]")
    output.add_argument("-dc", "--decimals", help="Number of decimal places in output [default: 2]", default=2)
    output.add_argument(
        "-hl",
        "--hide_link_headers",
        help="Hide alignment column headers",
        action="store_true",
    )
    output.add_argument(
        "-ha",
        "--hide_aln_headers",
        help="Hide alignment cluster name headers",
        action="store_true",
    )

    viz = parser_magm.add_argument_group("Visualisation options")
    viz.add_argument(
        "-ufo",
        "--use_file_order",
        action="store_true",
        help="Display clusters in order of input files"
    )

    args = parser.parse_args()
    return args

