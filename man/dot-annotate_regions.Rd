% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clipseq_related_functions.R
\name{.annotate_regions}
\alias{.annotate_regions}
\title{Single function to annotate by utr3, utr5, cds and introns}
\usage{
.annotate_regions(x, utr3_file, utr5_file, intron_file, coding_exon_file)
}
\arguments{
\item{x}{an object of class GRanges}
}
\value{

}
\description{
Single function to annotate by utr3, utr5, cds and introns
}
\examples{

intron_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_introns.bed", package = "SpongeAnalysis")

x <- rtracklayer::import(intron_bed_file)
utr3_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_3UTR.bed", package = "SpongeAnalysis")
utr5_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_5UTR.bed", package = "SpongeAnalysis")
cds_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_coding_exons.bed", package = "SpongeAnalysis")

.annotate_regions(x , utr3_file,utr5_file,intron_file,coding_exon_file)

}
\keyword{internal}
