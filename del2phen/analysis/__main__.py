import argparse
import sys

from del2phen.analysis.analyze import analyze
from del2phen.analysis.patient_comparison import (predict_phenotypes_for_cnv_strings,
                                                  _convert_cnv_str_to_cnv)
from del2phen.analysis.data_objects import Patient


class CustomHelp(argparse.HelpFormatter):
    """Custom help formatter_class that only displays metavar once."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar
        parts = []
        if action.nargs == 0:
            parts.extend(action.option_strings)
        else:
            default = self._get_default_metavar_for_optional(action)
            args_string = self._format_args(action, default)
            for option_string in action.option_strings:
                parts.append("%s" % (option_string))
            parts[-1] += " %s " % args_string
        return ", ".join(parts)

    def _format_action(self, action):
        parts = super()._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts


def _setup_argparser():
    parser = argparse.ArgumentParser()
    parser.prog = "Del2Phen"
    parser.formatter_class = CustomHelp
    parser.fromfile_prefix_chars = "@"

    # Data.
    data_args = parser.add_argument_group(title="Input Data Arguments")
    data_args.add_argument(
        "-g", "--genotypes",
        help="TSV file with 5 columns, labeled as follows with the "
             "relevant CNV information:"
             "\n  -id: patient ID"
             "\n  -chromosome: CNV contig name as listed in the relevant GTF file"
             "\n  -start: 1-indexed first affected base of the CNV"
             "\n  -stop: 1-indexed last affected base of the CNV (inclusive)"
             "\n  -copy_number: integer copy number (e.g. 1 is a single "
             "deletion, 3 is a single duplication)"
             "\nMultiple entries are allowed per patient ID."
        )
    data_args.add_argument(
        "-p", "--phenotypes",
        help="TSV file in which the first column contains patient IDs and the "
             "remaining columns contain 'T' or 'F' to indicate presence or absence of "
             "a phenotype. Column headers must be present, where column 1 is 'id', and "
             "the remaining column headers are HPO term IDs (formatted e.g. "
             "HPO:000001) or custom term IDs defined in --custom-phenotypes."
        )
    data_args.add_argument(
        "-cp", "--custom-phenotypes",
        dest="custom_phenotype_file",
        help="TSV file containing IDs and labels of custom non-HPO phenotypes included "
             "in the --phenotypes file header."
        )

    # data_args.add_argument(
    #     "--expand-hpos", dest="expand_hpo_terms",
    #     help="Recursively add parent HPO terms to each case based on presence of child "
    #          "HPO terms."
    #     )

    # Filtering.
    filtering_args = parser.add_argument_group(title="Filtering")
    filtering_args.add_argument(
        "-d", "--drop-list", dest="drop_list_file",
        help="Text file listing case IDs to skip, formatted with one ID per line. Empty "
             "lines are allowed and lines beginning with '#' are skipped. Cases listed "
             "will not be included in the set of cases analyzed."
        )
    filtering_args.add_argument(
        "-c", "--included-contigs", type=str,
        dest="chromosomes", nargs="+",
        help="Space-separated list of contig names to be included in analysis. Only "
             "CNVs on these contigs will be analyzed. Any cases without CNVs on any of "
             "these contigs are removed before analysis, unless --keep-uncompared is "
             "used. Default behavior is to analyze all contigs."
        )
    filtering_args.add_argument(
        "-cn", "--included-copy-numbers", type=int,
        dest="copy_numbers", nargs="+",
        help="Space-separated list of CNV copy number integers to be included in "
             "analysis. Only CNVs with one of the listed copy numbers will be "
             "analyzed. Any cases without any CNVs of the listed copy number are "
             "removed before analysis, unless --keep-uncompared is used. Default "
             "behavior is to analyze all CNV copy numbers."
        )
    filtering_args.add_argument(
        "-py", "--phenotype-termset",
        dest="phenotype_termset_yaml",
        help="YAML file specifying which phenotypes to use. See full documentation for "
             "specifications on how to build a YAML phenotype file, including custom "
             "combination phenotypes."
        )

    filtering_args.add_argument(
        "--keep-ungenotyped", action="store_true",
        help="Include cases that have no listed CNVs, but who are still listed in the "
             "phenotype file. Default behavior is to ignore any cases that do not have "
             "both CNV and phenotype information."
        )
    filtering_args.add_argument(
        "--keep-unphenotyped", action="store_true",
        help="Include cases that have no phenotype information, but who have CNV "
             "information. Default behavior is to ignore any cases that do not have "
             "both CNV and phenotype information."
        )
    filtering_args.add_argument(
        "--keep-uncompared", action="store_true",
        help="Include cases that will not be compared due to not having a CNV of a "
             "listed CNV type on a listed contig. Cases will not be included in "
             "comparative analysis, but will be included in summary statistics."
        )

    # Gene settings.
    gene_args = parser.add_argument_group(title="Gene Arguments")
    gene_args.add_argument(
        "-pli", "--pli-threshold", dest="pLI_threshold", default=0.9,
        help="Set minimum pLI threshold for defining HI genes. Default: 0.9"
        )
    gene_args.add_argument(
        "-hi", "--hi-threshold", dest="HI_threshold", default=10,
        help="Set maximum HI Score threshold for defining HI genes. Default: 10"
        )
    gene_args.add_argument(
        "-phaplo", "--phaplo-threshold", default=0.86,
        help="Set minimum pHaplo score for defining HI genes. Default: 0.86"
        )
    gene_args.add_argument(
        "-m", "--hi-mode", dest="mode", default="confirm",
        help="Set HI gene definition mode from one of the following: 'any': HI genes "
             "meet the threshold of at least one HI score. '2': HI genes meet the "
             "threshold at at least 2 HI scores. 'all': HI genes meet the threshold of "
             "all 3 HI scores. 'confirm': HI genes meet the threshold of at least 2/3 "
             "HI scores when at least 2 scores are known, otherwise HI genes meet the "
             "threshold of 1 HI score when only 1 is known. Default: confirm"
        )
    gene_args.add_argument(
        "-de", "--dominant-effect-genes", dest="dominant_gene_list", nargs="+",
        help="Space-separated list of gene IDs to treat as fully penetrant. By default, "
             "a full match of these genes in CNV gene content is required for group "
             "membership in the comparative analysis. This behavior can be disabled "
             "using --allow-de-gene-mismatch. Overriden by --de-file."
        )
    gene_args.add_argument(
        "--de-file", dest="dominant_gene_file",
        help="Text file listing gene IDs to treat as fully penetrant. By default, "
             "a full match of these genes in CNV gene content is required for group "
             "membership in the comparative analysis. This behavior can be disabled "
             "using --allow-de-gene-mismatch. Overrides the -de argument."
        )

    # Comparison threshold arguments.
    sim_thresholds = parser.add_argument_group(
        title="Similarity Thresholds",
        description="CNV similarity thresholds. Each metric is used to filter patient "
                    "subgroups based on Jaccard index similarity score between 0 and 1."
        )
    sim_thresholds.add_argument(
        "--length-sim", dest="length_similarity", default=0, type=float,
        help="CNV similarity by combined raw length of affected loci."
        )
    sim_thresholds.add_argument(
        "--loci-sim", dest="loci_similarity", default=0, type=float,
        help="CNV similarity by total overlap of affected loci."
        )
    sim_thresholds.add_argument(
        "--gene-sim", dest="gene_similarity", default=0, type=float,
        help="CNV similarity by all genes overlapping affected loci."
        )
    sim_thresholds.add_argument(
        "--hi-gene-sim", dest="hi_gene_similarity", default=0, type=float,
        help="CNV similarity by all haploinsufficient genes overlapping affected loci."
        )
    sim_thresholds.add_argument(
        "--allow-de-gene-mismatch", dest="dom_gene_match", action="store_false",
        help="Allow patients to form groups when dominant effects genes do not match. "
             "Default behavior is to only group patients when they have matching "
             "affected dominant effect genes (including when they each have none)."
        )

    prediction_args = parser.add_argument_group(title="Prediction Arguments")
    prediction_args.add_argument(
        "-abs", "--absolute-threshold", default=2, dest="abs_threshold",
        help="Minimum absolute occurrence of a phenotype within a patient group required "
             "for the phenotype to be considered 'predicted'. Default=2"
        )
    prediction_args.add_argument(
        "-rel", "--relative-threshold", default=0.2, dest="rel_threshold",
        help="Minimum relative frequency of a phenotype within a patient group required "
             "for the phenotype to be considered 'predicted'. Default=0.2"
        )
    prediction_args.add_argument(
        "--ignore-nas", action="store_true", dest="use_adjusted_frequency",
        help="When calculating relative frequency, ignore NA phenotype responses when "
             "calculating the denominator, i.e., only consider True/(True+False) "
             "instead of True/(True+False+NA)."
        )
    prediction_args.add_argument(
        "--group-size", default=5, dest="group_size_threshold",
        help="Minimum patient group size required to predict phenotypes. Default=5"
        )

    # Output arguments.
    output_args = parser.add_argument_group(title="Output Arguments")
    output_args.add_argument(
        "-op", "--output-predictions",
        help="Directory where predictions will be written. Must be specified to "
             "produce prediction files."
        )
    output_args.add_argument(
        "-os", "--output-stats",
        help="TSV file where prediction accuracy metrics will be written. Must be "
             "specified to produce stats file. Ignored and no stats will be produced "
             " if --cnv-query-only is also specified."
        )
    output_args.add_argument(
        "-cnv", "--cnv-predict", nargs="+", action="extend",
        help="Include phenotype predictions for the specified CNV. CNV string format "
             "is chr:start-stop:copy_number. Multiple CNVs can be specified "
             "sequentially separated by spaces or by passing this argument multiples "
             "times. Note that predictions can also be made by listing CNVs in the "
             "--genotypes file, specifying --keep-unphenotyped, and by omitting the "
             "--cnv-query-only option."
        )
    output_args.add_argument(
        "-x", "--cnv-query-only", action="store_true",
        help="Output predictions only for CNVs specified using --cnv-predict, and not "
             "patients/CNVs present in --genotypes."
        )

    # Reference.
    reference_args = parser.add_argument_group(
        title="Reference Arguments",
        description="Arguments to override the default reference files used for "
                    "analysis. NOTE: Using any of the following arguments will disable "
                    "the defaults for all four reference files to prevent accidental "
                    "application of incorrect values. See README for more information.")
    reference_args.add_argument("--gtf-file",
                                help="File path to a GTF format file containing gene "
                                     "information to construct the gene set for "
                                     "analysis.")
    reference_args.add_argument("--pli-file",
                                help="File path to a gnomad file containing pLI gene "
                                     "info.")
    reference_args.add_argument("--hi-file",
                                help="File path to a file containing HI gene info.")
    reference_args.add_argument("--phaplo-file",
                                help="File path to a file containing pHaplo gene info.")

    return parser


def main():
    parser = _setup_argparser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Compare provided patient data.
    comparison, geneset, ontology, termset = analyze(**vars(args))

    # If no output requested, finish.
    if not (args.output_predictions or args.output_stats):
        return comparison, geneset, ontology, termset

    # If cnv_query_only, make predictions for the command line CNVs, output a table,
    # and finish.
    if args.cnv_query_only and args.cnv_predict:
        predictions = predict_phenotypes_for_cnv_strings(
            cnv_strings=args.cnv_predict,
            geneset=geneset, comparison_database=comparison, tabulate=True,
            **vars(args))
        predictions.to_csv(f"{args.output_predictions}/cnv_query_predictions.tsv",
                           sep="\t")
        return comparison, geneset, ontology, termset, predictions

    # Otherwise add the CNVs as a patient to the comparisons.
    elif args.cnv_predict:
        cnvs = [_convert_cnv_str_to_cnv(cnv) for cnv in args.cnv_predict]
        for cnv in cnvs:
            cnv.set_genes(geneset)
        patient = Patient("cnv_predict")
        patient.cnvs = cnvs
        comparison.compare_patient_vs_others(patient, save_results=True, **vars(args))

    predictions = (comparison.test_all_patient_pheno_predictions(**vars(args),
                                                                 skip_no_hpos=False,
                                                                 filter_unknowns=False))
    if args.output_stats:
        stats_table = predictions.make_individual_precision_table(termset=termset,
                                                                  **vars(args))
        stats_table.to_csv(args.output_stats, sep="\t")
    if args.output_predictions:
        predictions.write_all_predictions(args.output_predictions, termset=termset,
                                          **vars(args))
    return comparison, geneset, ontology, termset, predictions


if __name__ == '__main__':
    results = main()
