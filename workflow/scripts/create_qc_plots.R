require(platetools)
require(ggplot2)
require(viridisLite)

make_plot <- function(df, column, title=NULL, limits=NULL) {
    p <- raw_map(data=df[column], well=df$well, plate=384) +
        ggtitle(title) +
        theme_dark() +
        scale_fill_viridis_c(limits=limits)
    
    p
}

create_plot <- function(data_file, out_file, log_file) {
    cat("Loading data\n", file=log_file)
    df <- read.delim(data_file)

    cat("Creating plots\n", file=log_file, append=TRUE)
    pdf(out_file, width=8, height=4)
    print(make_plot(df, "n_reads_mil_R1", title="Number of reads (Millions)", limits=c(0,5)))
    print(make_plot(df, "qual_q30_R1", title="Percent of Reads with Avg Qual >= 30 (read 1)", limits=c(0, 100)))
    print(make_plot(df, "qual_q30_R2", title="Percent of Reads with Avg Qual >= 30 (read 2)", limits=c(0, 100)))
    print(make_plot(df, "qual_avg_R1", title="Average Quality Score (read 1)", limits=c(0, 40)))
    print(make_plot(df, "qual_avg_R2", title="Average Quality Score (read 2)", limits=c(0, 40)))
    print(make_plot(df, "percent_uniq_map", title="Percent of uniquely mapped reads", limits=c(0, 100)))
    print(make_plot(df, "rrna", title="Percent of ribosomal RNA reads", limits=c(0, 100)))
    print(make_plot(df, "mito", title="Percent of mitochondrial RNA reads", limits=c(0, 100)))
    print(make_plot(df, "ercc", title="Percent of ERCC spike-in reads", limits=c(0, 100)))
    print(make_plot(df, "non_annot", title="Percent of reads in non-annotated space", limits=c(0, 100)))
    dev.off()
}

create_plot(
    snakemake@input[["pt_data"]],
    snakemake@output[["pdf"]],
    snakemake@log[["log_file"]]
)
