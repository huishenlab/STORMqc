require(platetools)
require(ggplot2)
require(viridisLite)
require(yaml)

make_plot <- function(df, column, title=NULL, limits=NULL) {
    p <- raw_map(data=df[column], well=df$well, plate=384) +
        ggtitle(title) +
        theme_dark() +
        scale_fill_viridis_c(limits=limits)
    
    p
}

create_plot <- function(data_file, limits_yaml, out_file, log_file) {
    cat("Loading data\n", file=log_file)
    df <- read.delim(data_file)

    cat("Loading limits configuration\n", file=log_file, append=TRUE)
    limits <- read_yaml(limits_yaml)

    cat("Creating plots\n", file=log_file, append=TRUE)
    pdf(out_file, width=8, height=4)
    print(
        make_plot(
            df,
            "n_reads_mil_R1",
            title="Number of reads (Millions)",
            limits=c(limits$n_reads$lo, limits$n_reads$hi)
        )
    )
    print(
        make_plot(
            df,
            "qual_q30_R1",
            title="Percent of Reads with Avg Qual >= 30 (read 1)",
            limits=c(limits$qual_q30_r1$lo, limits$qual_q30_r1$hi)
        )
    )
    print(
        make_plot(
            df,
            "qual_q30_R2",
            title="Percent of Reads with Avg Qual >= 30 (read 2)",
            limits=c(limits$qual_q30_r2$lo, limits$qual_q30_r2$hi)
        )
    )
    print(
        make_plot(
            df,
            "qual_avg_R1",
            title="Average Quality Score (read 1)",
            limits=c(limits$qual_avg_r1$lo, limits$qual_avg_r1$hi)
        )
    )
    print(
        make_plot(
            df,
            "qual_avg_R2",
            title="Average Quality Score (read 2)",
            limits=c(limits$qual_avg_r2$lo, limits$qual_avg_r2$hi)
        )
    )
    print(
        make_plot(
            df,
            "percent_uniq_map",
            title="Percent of uniquely mapped reads",
            limits=c(limits$percent_uniq_map$lo, limits$percent_uniq_map$hi)
        )
    )
    print(
        make_plot(
            df,
            "rrna",
            title="Percent of ribosomal RNA reads",
            limits=c(limits$rrna$lo, limits$rrna$hi)
        )
    )
    print(
        make_plot(
            df,
            "mito",
            title="Percent of mitochondrial RNA reads",
            limits=c(limits$mito$lo, limits$mito$hi)
        )
    )
    print(
        make_plot(
            df,
            "ercc",
            title="Percent of ERCC spike-in reads",
            limits=c(limits$ercc$lo, limits$ercc$hi)
        )
    )
    print(
        make_plot(
            df,
            "non_annot",
            title="Percent of reads in non-annotated space",
            limits=c(limits$non_annotated$lo, limits$non_annotated$hi)
        )
    )
    dev.off()
}

create_plot(
    snakemake@input[["pt_data"]],
    snakemake@params[["limits"]],
    snakemake@output[["pdf"]],
    snakemake@log[["log_file"]]
)
