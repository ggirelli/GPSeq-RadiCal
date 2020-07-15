#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# License: MIT - Copyright (c) 2020 Gabriele Girelli
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

depends_on = c("argparser", "data.table", "logging",
    "outliers", "pbapply", "rtracklayer")
tmp = lapply(depends_on, function(x) {
    expr = x %in% names(installed.packages()[, "Package"])
    if (!expr) {
        cat(sprintf(" Missing '%s' package.\n", x),
            sprintf("Install it with 'install.packages(\"%s\")'.\n", x))
    }
    stopifnot(expr)
})
pbapply::pboptions(type="timer")

# FUNCTIONS ====================================================================

chrom_to_chrom_id = function(chrom, nchrom = 24, hetero = c("X", "Y")) {
    # Convert a chromosome name (e.g., "chr1", "chrX") to a numerical ID.
    if (grepl(":", chrom)) {
        return(floor(as.numeric(gsub(":", ".",
            substr(chrom, 4, nchar(chrom))))))
    } else {
        chrom_id = substr(chrom, 4, nchar(chrom))
        if (chrom_id %in% hetero)
            chrom_id = nchrom - which(rev(hetero) == chrom_id) + 1
        chrom_id = as.numeric(chrom_id)
        stopifnot(!is.na(chrom_id))
        return(chrom_id)
    }
}

add_chrom_id = function(data, key = "chrom") {
    # Add chromosome ID to a data.table.
    # key should be the name of the column with chromosome names.
    stopifnot(data.table::is.data.table(data))
    stopifnot("chrom" %in% colnames(data))

    cid_table = data.table::data.table(
        chrom = as.character(unique(data$chrom)))
    cid_table$chrom_id = unlist(lapply(
        cid_table$chrom, FUN = chrom_to_chrom_id))
    data.table::setkeyv(cid_table, "chrom")
    data.table::setkeyv(data, key)

    return(data[cid_table, , nomatch = 0])
}

assert = function(expr, msg) {
    if (!expr) logging::logerror(msg)
    stopifnot(expr)
}

btag2extremes = function(btag) {
    if (grepl(":", btag)) {
        extremes = as.numeric(unlist(strsplit(btag, ":")))
    } else {
        extremes = rep(as.numeric(btag), 2)
    }
    extremes = data.table::as.data.table(t(extremes))
    colnames(extremes) = c("size", "step")
    extremes$tag = btag
    return(extremes)
}

bstring2specs = function(bstring) {
    taglist = unlist(strsplit(bstring, ","))
    binspecs = data.table::rbindlist(lapply(taglist, btag2extremes))
    return(binspecs)
}

mkbins = function(brid, bspecs, cinfo, elongate_ter_bin=FALSE) {
    bins = cinfo[, .(start=seq(start, end, by=bspecs[brid, step]),
        size=end), by=chrom]
    bins[, end := start + bspecs[brid, size] - 1]
    if (elongate_ter_bin) {
        bins = data.table::rbindlist(list(bins[end<=size, .(chrom, start, end)],
            bins[end > size, .(start=min(start), end=min(end)), by=chrom]))
    } else {
        bins = data.table::rbindlist(list(bins[end<=size, .(chrom, start, end)],
            bins[end > size, .(start=min(start), end=size[1]), by=chrom]))
    }
    bins = add_chrom_id(bins)
    bins[, chrom := reorder(chrom, chrom_id)]
    bins[, chrom_id := NULL]
    bins = bins[order(chrom, start)]
    bins$tag = bspecs[brid, sprintf("%.0e:%.0e", size, step)]
    return(bins)
}

import_gpseq_bed = function(brid, bmeta) {
    assert("fpath" %in% colnames(bmeta),
        "Missing 'fpath' column in metadata file.")
    brmeta = data.table::copy(bmeta[brid])
    o = data.table::as.data.table(rtracklayer::import.bed(brmeta$fpath))
    data.table::setnames(o, "seqnames", "chrom")
    o[, c("width", "strand", "name") := NULL]
    brmeta[, fname := NULL]
    o[, condition := sprintf("cid_%d", brid)]
}

parse_bed_meta = function(bbmeta, args) {
    logging::loginfo(sprintf("Reading %d bed files.", nrow(bbmeta)))
    bd = data.table::rbindlist(pbapply::pblapply(seq_len(nrow(bbmeta)),
        import_gpseq_bed, bbmeta, cl=args$threads))
    bd[, end := start]
    logging::loginfo(sprintf("Dcasting bed data."))
    bd = data.table::dcast(bd, chrom+start+end~condition,
        value.var="score", fill=0)
    if (3 <= args$export_level) {
        logging::loginfo("Exporting dcasted input bed...")
        saveRDS(bd, file.path(args$exp_output_folder, "input_bed.rds"))
    }
    return(bd)
}

otag2specs = function(otag) {
    outlier_specs = unlist(strsplit(otag, ":"))
    outlier_specs = outlier_specs[1:min(2, length(outlier_specs))]
    outlier_specs = data.table::as.data.table(t(outlier_specs))
    colnames(outlier_specs) = c("method", "threshold")
    outlier_specs[, method := tolower(method)]
    outlier_specs[, threshold := as.numeric(threshold)]
    assert(outlier_specs$method %in% outlier_methods,
        sprintf("Unrecognized outlier method '%s'.", outlier_specs$method))
    outlier_specs$alpha = .5  # Default
    outlier_specs$lim = 1.5   # Default
    if (outlier_specs$method == "iqr") {
        outlier_specs[, lim := threshold]
    } else {
        outlier_specs[, alpha := threshold]
    }
    outlier_specs[, threshold := NULL]
    return(outlier_specs)
}

export_output = function(odata, odir, format, suffix, rm_tag=FALSE) {
    assert("tag" %in% colnames(odata), "Cannot find required 'tag' column.")
    assert(format %in% c("tsv", "tsv.gz", "csv", "csv.gz", "rds"),
        sprintf("Unrecognized output format '%s'.", format))
    opath_base = sprintf("%s.bins_%s", suffix, gsub(":", "_", odata[1, tag]))
    columns_to_export = colnames(odata)
    if (rm_tag) {
        columns_to_export = columns_to_export["tag" != columns_to_export]
    }
    if ("rds" == format) {
        saveRDS(odata[, .SD, .SDcols=columns_to_export],
            file.path(odir, sprintf("%s.rds", opath_base)))
    }
    if (format %in% c("tsv", "tsv.gz")) {
        data.table::fwrite(odata[, .SD, .SDcols=columns_to_export],
            file.path(odir, sprintf("%s.%s", opath_base, format)),
            sep="\t", na="NA", quote=F)
    }
    if (format %in% c("csv", "csv.gz")) {
        data.table::fwrite(odata[, .SD, .SDcols=columns_to_export],
            file.path(odir, sprintf("%s.%s", opath_base, format)),
            sep=",", na="NA", quote=F)
    }
}

get_condition_outliers_stats = function(x) {
    x = x[0 != x]
    is_outlier = outliers::scores(x,
        type=bed_outlier_specs$method,
        prob=1-bed_outlier_specs$alpha,
        lim=bed_outlier_specs$lim)
    ostats = as.numeric(summary(x[is_outlier]))
    ostats = c(ostats, sum(is_outlier), length(x))
    return(ostats)
}

rm_condition_outliers = function(x) {
    is_outlier = outliers::scores(x[0 != x],
        type=bed_outlier_specs$method,
        prob=1-bed_outlier_specs$alpha,
        lim=bed_outlier_specs$lim)
    x[which(0 != x)[is_outlier]] = 0
    return(x)
}

calc_bed_outliers_stats = function(bd, args) {
    logging::loginfo(sprintf(
        "Calculating outlier stats. [%s]", args$bed_outlier_tag))
    bed_outlier_specs = otag2specs(args$bed_outlier_tag)
    outlier_stats = bd[, lapply(.SD,
        get_condition_outliers_stats), .SDcols=args$cond_cols]
    outlier_stats = data.table::data.table(t(outlier_stats))
    colnames(outlier_stats) = c(names(summary(1)), "nout", "ntot")
    outlier_stats_opath = file.path(args$exp_output_folder, "outlier_stats.tsv")
    logging::loginfo(sprintf(
        "Exporting outlier stats to '%s'.", outlier_stats_opath))
    data.table::fwrite(outlier_stats, outlier_stats_opath, sep="\t")
}

rm_bed_outliers = function(bd, args) {
    logging::loginfo(sprintf(
        "Removing outliers. [%s]", args$bed_outlier_tag))
    bd[, c(args$cond_cols) := lapply(.SD,
        rm_condition_outliers), .SDcols=args$cond_cols]
    bd = bd[0 != apply(bd[, .SD, .SDcols=args$cond_cols], MARGIN=1, FUN=sum)]
    if (3 <= args$export_level) {
        logging::loginfo(
            "Exporting dcasted input bed after outlier removal...")
        saveRDS(bd, file.path(args$exp_output_folder, "clean_bed.rds"))
    }
    return(bd)
}

bin_bed_data = function(bbins, cond_cols, bd, site_domain, site_universe=NULL) {
    logging::loginfo(sprintf("Binning... [%s]", bbins[1, tag]))
    binned = data.table::rbindlist(pbapply::pblapply(split(bd, bd$chrom),
        bin_chromosome, cond_cols, bbins, site_domain, site_universe
        ))[order(tag, chrom, start, cid)]
    return(binned)
}

bin_chromosome = function(bbd,
        cond_cols, bins, site_domain, site_universe=NULL) {
    data.table::setkeyv(bins, bed3_colnames)

    selected_chromosome = bbd[1, chrom]
    ovlps = data.table::foverlaps(bbd, bins)

    nreads = ovlps[, lapply(.SD, sum), by=c(bed3_colnames, "tag"),
        .SDcols=cond_cols][order(tag, chrom, start)]
    nreads = data.table::melt(nreads, id.vars=c(bed3_colnames, "tag"))
    data.table::setnames(nreads, c("variable", "value"), c("cid", "nreads"))
    nreads[, cid := match(cid, cond_cols)]
    data.table::setkeyv(nreads, c(bed3_colnames, "tag", "cid"))

    if ("universe" == site_domain) {
        assert(!is.null(site_universe),
            "Missing site universe data with site domain 'universe'.")
        nsites = data.table::foverlaps(
            site_universe, bins[chrom==selected_chromosome]
            )[!is.na(start), .(
                tag=bins[1, tag], cid=seq_len(cond_cols), nsites=.N
            ), by=bed3_colnames]
    } else {
        if ("union" == site_domain) {
            nsites = ovlps[, lapply(.SD, function(x) length(x)),
                by=c(bed3_colnames, "tag"), .SDcols=cond_cols
                ][order(tag, chrom, start)]
        } else {
            nsites = ovlps[, lapply(.SD, function(x) sum(0 != x)),
                by=c(bed3_colnames, "tag"), .SDcols=cond_cols
                ][order(tag, chrom, start)]
        }
        nsites = data.table::melt(nsites, id.vars=c(bed3_colnames, "tag"))
        data.table::setnames(nsites, c("variable", "value"), c("cid", "nsites"))
        nsites[, cid := match(cid, cond_cols)]
    }
    data.table::setkeyv(nsites, c(bed3_colnames, "tag", "cid"))

    combined = nreads[nsites]
    combined[, lib_nreads := as.numeric(total_lib_nreads)[cid]]
    combined[, chr_nreads := as.numeric(total_chr_nreads[
        selected_chromosome==chrom, .SD, .SDcols=cond_cols])[cid]]

    return(combined)
}

export_binned_bed_data = function(binned, odir, format="rds") {
    export_output(binned, odir, format, "binned")
}

mask_track = function(bbins) {
    if ("chrom:wide" == bbins[1, tag]) return(bbins)
    data.table::setkeyv(bbins, bed3_colnames)

    masked = unique(data.table::foverlaps(bbins, mask)[,
        .(tag, nreads, nsites, lib_nreads, chr_nreads,
            mask_overlaps=.N, mask_overlapped=!is.na(start)),
        by=c(bed3_colnames[1], paste0("i.", bed3_colnames[2:3]),  "cid")])
    data.table::setnames(masked,
        paste0("i.", bed3_colnames[2:3]), bed3_colnames[2:3])

    masked[!(mask_overlapped), mask_overlaps := 0]
    masked[, mask_overlapped := NULL]
    masked[0 < mask_overlaps, c("nreads", "nsites") := .(0, 0)]
    masked[, mask_overlaps := NULL]

    return(masked)
}

export_masked_data = function(masked, odir, format="rds") {
    export_output(masked, odir, format, "binned_masked")
}

estimate_centrality = function(bbind, normalize_by) {
    bbind = bbind[,
        .(pres=nreads/nsites/get(sprintf("%s_nreads", normalize_by))),
        by=c(bed3_colnames, "tag", "cid")]
    centr = bbind[order(chrom, start, cid), .(
            score=sum(pres[2:.N]/pres[1:(.N-1)])
        ), by=c(bed3_colnames, "tag")]
    return(centr)
}

export_estimated_centrality = function(centr, odir, format="rds") {
    export_output(centr, odir, format, "estimated")
}

calc_rescale_factor_by_lib = function(estmd, specs) {
    estmd[!is.na(score), outlier_status := outliers::scores(
        estmd[!is.na(score), score], type=score_outlier_specs$method,
        prob=1-score_outlier_specs$alpha, lim=score_outlier_specs$lim)]
    lowest = estmd[FALSE == outlier_status, min(score, na.rm = T)]
    highest = estmd[FALSE == outlier_status, max(score-lowest, na.rm = T)]
    estmd[, outlier_status := NULL]
    return(c(lowest, highest))
}

rescale_by_lib = function(estmd, specs) {
    rescaling_factors = calc_rescale_factor_by_lib(estmd, specs)
    estmd[, score := .(score - rescaling_factors[1])]
    estmd[, score := .(score / rescaling_factors[2])]
    estmd[, score := 2**score]
    return(estmd)
}

rescale_by_chr = function(estmd, specs) {
    data.table::rbindlist(by(estmd, estmd$chrom, rescale_by_lib, specs))
}

export_rescaled_centrality = function(rscld, odir, format="tsv.gz") {
    export_output(rscld, odir, format, "rescaled", rm_tag=TRUE)
}

apply_intersection_site_domain = function(bd, args) {
    logging::loginfo("Retaining only sites shared across conditions...")
    n_condition_empty = rowSums(0 == bd[, .SD, .SDcols=args$cond_cols])
    bd = bd[0 == n_condition_empty]
    logging::loginfo(sprintf("Retained %d/%d (%.1f%%) common sites.",
        nrow(bd), length(n_condition_empty),
        nrow(bd)/length(n_condition_empty)*100))
    if (3 <= args$export_level) {
        logging::loginfo(
            "Exporting dcasted input bed after site intersection...")
        saveRDS(bd, file.path(args$exp_output_folder,
            "clean_bed.intersected.rds"))
    }
    return(bd)
}

mask_binned = function(binned, args) {
    assert(file.exists(args$mask_bed),
        sprintf("Cannot find mask bed file '%s'.", args$mask_bed))
    mask = data.table::as.data.table(
        rtracklayer::import.bed(args$mask_bed))[, .(chrom=seqnames, start, end)]
    data.table::setkeyv(mask, bed3_colnames)
    if (args$chromosome_wide) {
        logwarn("Skipped masking for chromosome-wide bins.")
    }
    binned = pbapply::pblapply(binned, mask_track, cl=args$threads)
    if (1 <= args$export_level) {
        logging::loginfo(sprintf("Exporting binned bed data..."))
        tmp = lapply(binned, export_masked_data, args$exp_output_folder)
    }
    return(binned)
}

rescale_estimated = function(estimated, args) {
    logging::loginfo(sprintf("Rescaling estimates... [%s]", args$normalize_by))
    score_outlier_specs = otag2specs(args$score_outlier_tag)
    if ("chr" == args$normalize_by) {
        if (args$chromosome_wide) logwarn(
            "Skipped rescaling by chromosome for chromosome-wide bins.")
        rescaled = pbapply::pblapply(estimated, function(estmd) {
            if ("chrom:wide" == estmd[1, tag]) return(estmd)
            estmd = rescale_by_chr(estmd, score_outlier_specs)
        }, cl=args$threads)
    }
    if ("lib" == args$normalize_by) {
        rescaled = pbapply::pblapply(estimated, function(estmd) {
            estmd = rescale_by_lib(estmd, score_outlier_specs)
        }, cl=args$threads)
    }
    logging::loginfo(sprintf("Exporting rescaled centrality..."))
    tmp = lapply(rescaled, export_rescaled_centrality, args$exp_output_folder)
    return(rescaed)
}

process_experiment = function(bbmeta, bins, cinfo, args) {
    exid = bbmeta[1, exid]
    logging::loginfo(sprintf("Processing experiment '%s'.", exid))
    args$exp_output_folder = file.path(args$output_folder, exid)
    dir.create(args$exp_output_folder)

    logging::loginfo("Storing metadata.")
    data.table::fwrite(bbmeta,
        file.path(args$exp_output_folder, "bed.metadata.tsv"), sep="\t")

    args$cond_cols = sprintf("cid_%d", seq_len(nrow(bbmeta)))

    # Read bed files -----------------------------------------------------------

        bd = parse_bed_meta(bbmeta, args)

    # Get outlier stats --------------------------------------------------------

        if (0 != nchar(args$bed_outlier_tag)) {
            calc_bed_outliers_stats(bd, args)
        } else {
            logging::loginfo(sprintf("Skipped outlier removal."))
        }

    # Clean bed outliers -------------------------------------------------------

        if (0 != nchar(args$bed_outlier_tag)) bd = rm_bed_outliers(bd, args)

    # Apply intersection site domain -------------------------------------------

        if ("intersection" == args$site_domain) {
            bd = apply_intersection_site_domain(bd, args)
        } else if ("universe" == args$site_domain) {
            logging::loginfo(sprintf(
                "Reading site bed file '%s'", args$site_bed))
            site_universe = data.table::as.data.table(
                rtracklayer::import.bed(args$site_bed))[,
                    .(chrom=seqnames, start, end=start)]
        }

    # Calculate normalization factors ------------------------------------------

        logging::loginfo(sprintf("Calculating normalization factors."))
        total_lib_nreads = bd[, lapply(.SD, sum), .SDcols=args$cond_cols]
        total_chr_nreads = bd[, lapply(.SD, sum),
            by=chrom, .SDcols=args$cond_cols]

    # Assign to bins -----------------------------------------------------------

        bin_tags = bins[, unique(tag)]
        binned = by(bins, bins$tag, bin_bed_data,
            args$cond_cols, bd, args$site_domain, site_universe)
        if (1 <= args$export_level) {
            logging::loginfo(sprintf("Exporting binned bed data..."))
            tmp = lapply(binned, export_binned_bed_data, args$exp_output_folder)
        }

    # Masking track ------------------------------------------------------------

        if (!is.na(args$mask_bed)) binned = mask_binned(binned, args)

    # Calculate centrality -----------------------------------------------------

        logging::loginfo(sprintf("Estimating centrality..."))
        estimated = pbapply::pblapply(
            binned, estimate_centrality, args$normalize_by, cl=args$threads)
        if (1 <= args$export_level) {
            logging::loginfo(sprintf("Exporting estimated centrality..."))
            tmp = lapply(estimated,
                export_estimated_centrality, args$exp_output_folder)
        }

    # Rescale estimates --------------------------------------------------------

        if (0 == nchar(args$score_outlier_tag)) {
            logging::loginfo(sprintf("Skipped rescaling."))
            logging::loginfo(sprintf("Exporting estimated centrality..."))
            tmp = lapply(estimated, export_estimated_centrality,
                args$exp_output_folder, format="tsv.gz")
        } else {
            rescaled = rescale_estimated(estimated, args)
        }
}

# COMMON PARAMETERS ============================================================

bed3_colnames = c("chrom", "start", "end")
outlier_methods = c("z", "t", "chisq", "iqr", "mad")

# INPUT ========================================================================

parser = argparser::arg_parser("
Provide the path to a 4-columns tabulation-separated metadata file, containing
the sequencing run ID, a condition description, the library ID, and full
absolute path to each BED file (possibly gzipped) from the GPSeq experiment.
The bed files should be reported in order of condition strength, with the top
rows being *weaker* than bottom ones. Also, the first line should contain the
column headers: exid, cond, libid, and path. An example metadata file would be
is available at
https://github.com/ggirelli/GPSeq-RadiCal/blob/master/example_meta.tsv.

Fore more details, use the '--more-help' option.
", name="gpseq-radical.R")

parser = argparser::add_argument(parser, arg="bmeta_path",
    help="Path to bed metadata tsv file.")
parser = argparser::add_argument(parser, arg="output_folder",
    help="Path to output directory. Stops if it already exists.")

parser = argparser::add_argument(parser, arg="--cinfo-path",
    help="Path to bed file with chromsome sizes. Queries UCSC if not provided.",
    default=NULL)
parser = argparser::add_argument(parser, arg="--ref-genome", help=paste0(
        "Used when --cinfo-path is not provided to query UCSC. ",
        "If the provided reference genome is not found, reverts to 'hg38'."),
    default="hg19")

parser = argparser::add_argument(parser, arg="--bin-tags",
    help="Comma-separated bin tags. Use --more-help for more details.",
    default="1e6:1e5,1e5")
parser = argparser::add_argument(parser, arg="--bed-outlier-tag",
    help="Method:threshold for input bed outlier removal.",
    default="chisq:0.01")
parser = argparser::add_argument(parser, arg="--score-outlier-tag",
    help="Method:threshold for output score outlier removal (rescaling).",
    default="iqr:1.5")

parser = argparser::add_argument(parser, arg="--normalize-by", help=paste0(
        "Whether rescaling should be performed library-wise ('lib') ",
        "or chromosome-wise ('chr')."),
    default="lib")

parser = argparser::add_argument(parser, arg="--site-domain",
    help="Site domain method. Use --more-help for more detail.",
    default="separate")
parser = argparser::add_argument(parser, arg="--site-bed", help=paste0(
        "Path to bed file with recognition site locations for the enzyme ",
        "in use. This is required when using '--site-domain universe'."),
    default=NULL)

parser = argparser::add_argument(parser, arg="--mask-bed",
    help="Path to bed file with regions to be masked.",
    default=NULL)

parser = argparser::add_argument(parser, arg="--threads", help=paste0(
        "Number of threads for parallelization. ",
        "Optimal when using at least one core per bed file."),
    default=1)
parser = argparser::add_argument(parser, arg="--export-level",
    help="Limits the amount of output. Use --more-help for more details",
    default=0)

parser = argparser::add_argument(parser, arg="--chromosome-wide", flag=TRUE,
    help="Use this option to calculate also on chromosome-wide bins.")
parser = argparser::add_argument(parser, arg="--elongate-ter-bin", flag=TRUE,
    help=paste0("Use this option to elongate chromosome-terminal bins ",
        "and have equally-sized bins."))
parser = argparser::add_argument(parser, arg="--more-help", flag=TRUE,
    help="Show extended help page and exit")

if ("--more-help" %in% commandArgs(trailingOnly=TRUE)) {
    cat("
Provide the path to a 4-columns tabulation-separated metadata file, containing
the sequencing run ID, a condition description, the library ID, and full
absolute path to each BED file (possibly gzipped) from the GPSeq experiment.
The bed files should be reported in order of condition strength, with the top
rows being *weaker* than bottom ones. Also, the first line should contain the
column headers: exid, cond, libid, and path. An example metadata file would be
is available at
https://github.com/ggirelli/GPSeq-RadiCal/blob/master/example_meta.tsv.

Binning can be specified by providing comma-separated bin 'tags', in the format
of 'bin_size:bin_step' or 'bin_size'. If 'bin_step' is not specified, the bins
are generated in a non-overlapping manner. Use the '--chromosome-wide' option to
estimated centrality also on chromosome-wide bins. Use the '--elongate-ter-bin'
option to elongate chromosome terminal bins to the specified bin size, the
default behaviour is for terminal bins' end to coincide with the chromosome end.

Outlier removal method can be specified with the 'method:threshold' format.
Available methods: z, t, chisq, iqr, and mad. The specified threshold is an
alpha threshold on the outlier score p-value for the z, t, chisq, and mad
methods, and should thus be 0 < threshold < 1. In the case of the iqr method,
outliers are identified by comparing their distance to the closest quartile (Q1
or Q3) and the product threshold*iqr. More details are available here:
https://www.rdocumentation.org/packages/outliers/versions/0.14/topics/scores

Outliers can be removed from input bed files by using the '--bed-outlier-tag'
option. To skip outlier removal use '--bed-outlier-tag \"\"'.

Moreover, centrality outliers are detected and used to rescale the final
estimates when using the '--score-outlier-tag' option. To skip score rescaling
use '--score-outlier-tag \"\"'. Rescaling can be performed in a
chromosome-wise or library-wise manner, by using the '--normalize-by chr' or
'--normalize-by lib' options. The '--normalize-by' is also used for centrality
calculation, even if rescaling is skipped.

The '--site-domain' option can be used to specify which recognition sites to
consider when estimating centrality. The default 'separate' domain considers all
(non-outlier) sites with at least one mapped read in a condition. The 'union'
domain considers all (non-outlier) sites with at least one mapped read in at
least one condition. The 'intersection' domain retains only recognition sites
with at least on mapped read in all conditions. The 'universe' domain considers
all known recognition site locations in the reference genome. When using
'--site-domain universe', a bed file with known recognition site locations must
be provided using the '--site-bed' option.

The output centrality estimation can be masked (masking performed after
binning) to removed known problematic regions. Use the '--mask-bed' option to
provide the path to a bed file with regions to be masked. A bin is masked if it
overlaps with any regions in the provided mask bed file.

The '--export-level' option is useful to limit the amount of output files
produced by the script. The default export level is set to 0.
 0: export only final results and statistics.
 1: export intermediate results.
 2: export supplementary (bins) files.
 3: export dcasted input (bed, clean bed, domained bed).

Finally, use the '--threads' option to specify how many threads to use for
parallelization. An optimal number of threads coincides with the number of
bed files in the input metadata file.
    \n")
    quit()
}

args = argparser::parse_args(parser)
assert(args$normalize_by %in% c("chr", "lib"),
    sprintf("Unrecognized 'normalize_by' value: '%s'", args$normalize_by))
assert(args$site_domain %in% c("separate", "union", "intersection", "universe"),
    sprintf("Unrecognized 'site_domain' value: '%s'", args$site_domain))
if ("universal" == args$site_domain) {
    assert(!is.na(args$site_bed),
        "A 'site_bed' path must be specified when 'site_domain' is 'universe'.")
    assert(file.exists(args$site_bed),
        sprintf("Cannot find 'site_bed' file '%s'.", args$site_bed))
}

# RUN ==========================================================================

# Prelude ----------------------------------------------------------------------

    assert(!dir.exists(args$output_folder), paste0(
        sprintf("Output folder '%s' already exists. ", args$output_folder),
        "This can lead to mixed settings output and should be avoided. ",
        "Stopping."))
    dir.create(args$output_folder)

    logging::basicConfig()
    log_path = file.path(args$output_folder, "gpseq-radicalc.log")
    logging::loginfo(sprintf("This log will be stored at '%s'.", log_path))
    logging::addHandler(logging::writeToFile, file=log_path)

    logging::loginfo(sprintf("Created output folder '%s'.", args$output_folder))
    settings_path = file.path(args$output_folder, "gpseq-radicalc.opts.rds")
    saveRDS(args, settings_path)
    logging::loginfo(sprintf(
        "Exported input parameters to '%s'.", settings_path))

    data.table::setDTthreads(args$threads)
    if (args$threads != data.table::getDTthreads()) {
        logwarn(sprintf("Changed from the requested %d to %d threads.",
            args$threads, data.table::getDTthreads()))
        args$threads = data.table::getDTthreads()
    }

# Parse metadata ---------------------------------------------------------------

    logging::loginfo(sprintf("Parsing metadata from '%s'.", args$bmeta_path))
    bmeta = data.table::fread(args$bmeta_path)
    assert(2 < nrow(bmeta), "Provide at least two bed files.")
    logging::loginfo("Storing metadata.")
    data.table::fwrite(bmeta,
        file.path(args$output_folder, "bed.metadata.tsv"), sep="\t")

# Read chromosome info bed -----------------------------------------------------

    cinfo = NULL
    if (is.na(args$cinfo_path)) {
        logging::loginfo("Opening UCSC browser session...")
        ucsc = rtracklayer::browserSession("UCSC")
        rtracklayer::genome(ucsc) = args$ref_genome
        logging::loginfo(sprintf("Querying UCSC for '%s' chromosome info...",
            rtracklayer::genome(ucsc)))
        cinfo = data.table::data.table(rtracklayer::getTable(
            rtracklayer::ucscTableQuery(ucsc,
                track="Chromosome Band", table="cytoBand")))
        cinfo = cinfo[, .(start=1, end=max(chromEnd)), by=chrom]
    } else {
        assert(file.exists(args$cinfo_path),
            sprintf("Cannot find chromosome info bed file '%s'.",
                args$cinfo_path))
        logging::loginfo(sprintf(
            "Reading chromosome info from '%s'.", args$cinfo_path))
        cinfo = data.table::as.data.table(
            rtracklayer::import.bed(args$cinfo_path))
        data.table::setnames(cinfo, "seqnames", "chrom")
        cinfo[, c("width", "strand") := NULL]
    }
    assert(!is.null(cinfo), "Failed to build or retrieve chromosome info.")

# Build bins -------------------------------------------------------------------

    logging::loginfo(sprintf("Building bins."))
    bspecs = bstring2specs(args$bin_tags)
    if (0 == nrow(bspecs)) {
        bins = data.table()
    } else {
        bins = data.table::rbindlist(pbapply::pblapply(seq_len(nrow(bspecs)),
            mkbins, bspecs, cinfo, args$elongate_ter_bin, cl=args$threads))
    }
    if (args$chromosome_wide) {
        args$chromosome_wide_bins = data.table::copy(cinfo)
        args$chromosome_wide_bins[, tag := "chrom:wide"]
        bins = data.table::rbindlist(list(args$chromosome_wide_bins, bins))
    }
    data.table::setkeyv(bins, bed3_colnames)
    assert(0 < nrow(bins), "No bins built. Stopping.")
    if (2 <= args$export_level) {
        logging::loginfo("Exporting bins...")
        saveRDS(bins, file.path(args$output_folder, "bins.rds"))
    }

# Process one experiment at a time ---------------------------------------------

    assert("exid" %in% colnames(bmeta), "Missing 'exid' column from metadata.")
    tmp = by(bmeta, bmeta$exid, process_experiment, bins, cinfo, args)

# ------------------------------------------------------------------------------

logging::loginfo(sprintf("Done."))

################################################################################
