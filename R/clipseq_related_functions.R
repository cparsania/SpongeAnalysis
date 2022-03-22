# given the list of peaks in bed file format
#
# 1) calculate average depth for each peak
#
# 2) annotate peaks by different genomic regions
#
# 3) assign genes to each peak
#
# 4) assign IGV search text
#
# 5) assign peak length, peak seq, GC
#
#



#' Annotate clipseq peaks
#'
#' @param x an object of class GRanges where each row denotes clipseq peaks. Object must have column 'peak_id'.
#'
#' @return an object of class GRanges
#' @export
#'
#' @examples
#'
#' peak_file <- system.file("extdata", "boudreau_et_al_clipseq_peaks.bed", package = "SpongeAnalysis")
#' peak_data <- readr::read_delim(peak_file,delim = "\t")
#' peak_data_gr <- peak_data %>% plyranges::as_granges()
#'
#' intron_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_introns.bed", package = "SpongeAnalysis")
#' utr3_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_3UTR.bed", package = "SpongeAnalysis")
#' utr5_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_5UTR.bed", package = "SpongeAnalysis")
#' cds_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_coding_exons.bed", package = "SpongeAnalysis")
#'
#' clipseq_annotate_peaks(x = peak_data_gr,coding_exon_file = cds_bed_file , utr3_file = utr3_bed_file, utr5_file = utr5_bed_file, intron_file = intron_bed_file)
#'
clipseq_annotate_peaks <- function(x, coding_exon_file , utr3_file, utr5_file, intron_file){

  # validate x class
  stopifnot("x must be the object of class 'GRanges' " = is(x , "GRanges"))

  # x must have column peak_id
  check_peak_id <- "peak_id" %in% (x %>% tibble::as_tibble() %>% colnames())
  stopifnot("column 'peak_id' must present in the x" = check_peak_id)

  # values in the column peak_id must be unique.
  check_peak_id_unique <- x %>% tibble::as_tibble() %>% dplyr::pull(peak_id) %>% duplicated() %>% any()
  stopifnot("values in the column peak_id must be unique" = !check_peak_id_unique)

  x <- x %>% dplyr::select(peak_id)

  # prepare reference regions
  ucsc_introns <- .annotate_introns(intron_file = intron_file,utr3_file = utr3_file, utr5_file = utr5_file,coding_exon_file = cds_bed_file)
  ucsc_coding_exon <- .get_ucsc_coding_exon(coding_exon_file = coding_exon_file)
  ucsc_utr5 <- .get_ucsc_utr5(utr5_file = utr5_file)
  ucsc_utr3 <- .get_ucsc_utr3(utr3_file = utr3_file)

  # combine them in a single object
  comb_annot <-  plyranges::bind_ranges(ucsc_introns,ucsc_coding_exon,ucsc_utr5, ucsc_utr3) %>% plyranges::as_granges()

  # map annotations to clip peaks
  with_annot <- x %>%
    plyranges::as_granges() %>%
    plyranges::join_overlap_left_within_directed(comb_annot, suffix = c("_peak", "_reference"))


  #  annotate peaks based on following priority
  # 3'UTR > CDS > 5'UTR > Introns >

  oo <- with_annot %>%
    tibble::as_tibble()  %>%
    dplyr::select(peak_id, type,name) %>%
    dplyr::group_by(peak_id) %>%
    dplyr::mutate(type_final = dplyr::case_when(any(type == "ucsc_utr3") ~ "ucsc_utr3",
                                                any(type == "ucsc_coding_exon") ~ "ucsc_coding_exon",
                                                any(type == "ucsc_utr5") ~ "ucsc_utr5",
                                                any(type == "intron_utr3") ~ "intron_utr3",
                                                any(type == "intron_utr5") ~ "intron_utr5",
                                                any(type == "non_utr_introns") ~ "non_utr_introns",
                                                TRUE ~ "other"))
  # get final peak locations
  peak_loc <- oo %>%
    dplyr::ungroup() %>%
    dplyr::select(peak_id, type_final) %>% unique() %>%
    dplyr::rename(reference_annotation = type_final)

  x <- x %>% tibble::as_tibble() %>% dplyr::left_join(peak_loc) %>% plyranges::as_granges()
  return(x)
}




#' Get annotation of coding exons
#'
#' @param coding_exon_file a character string pointing coding regions in a bed file format
#'
#' @return an object of GRanges
#' @export
#'
#' @examples
#'
#' cds_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_coding_exons.bed", package = "SpongeAnalysis")
#' .get_ucsc_coding_exon(cds_bed_file)
#'
#' @keywords internal
.get_ucsc_coding_exon <- function(coding_exon_file){
  x <- rtracklayer::import(coding_exon_file)
  x %<>% dplyr::mutate(type = "ucsc_coding_exon")
  return(x)
}

#' Get annotation of utr5 regions
#'
#' @param utr5_file a character string pointing utr5 regions in a bed file format
#'
#' @return
#' @export
#'
#' @examples
#'
#' utr5_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_5UTR.bed", package = "SpongeAnalysis")
#' .get_ucsc_utr5(utr5_bed_file)
#'
#' @keywords internal
#'
.get_ucsc_utr5 <- function(utr5_file){

  x <- rtracklayer::import(utr5_file)
  x %<>% dplyr::mutate(type = "ucsc_utr5")
  return(x)
}

#' Get annotation of utr3 regions
#'
#' @param utr3_file a character string pointing utr3 regions in a bed file format
#'
#' @return
#' @export
#'
#' @examples
#' utr3_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_3UTR.bed", package = "SpongeAnalysis")
#' .get_ucsc_utr3(utr3_bed_file)
#'
#' @keywords internal
#'
.get_ucsc_utr3 <- function(utr3_file){
  x <- rtracklayer::import(utr3_file)
  x %<>% dplyr::mutate(type = "ucsc_utr3")
  return(x)
}


#' Get annotation of intron regions
#'
#' @param intron_file character string pointing intron regions in a bed file format
#'
#' @return
#' @export
#'
#' @examples
#' intron_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_introns.bed", package = "SpongeAnalysis")
#' .get_ucsc_introns(intron_bed_file)
#'
#' @keywords internal
.get_ucsc_introns <- function(intron_file){

  x <- rtracklayer::import(intron_file)
  x %<>% dplyr::mutate(type = "ucsc_intron")
  return(x)
}

#' categorize introns into introns (coding + non-coding), intron_utr3 and intron_utr5
#'
#' @param intron_file character string pointing intron regions in a bed file format
#'
#' @return
#' @export
#'
#' @examples
#' intron_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_introns.bed", package = "SpongeAnalysis")
#' utr3_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_3UTR.bed", package = "SpongeAnalysis")
#' utr5_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_5UTR.bed", package = "SpongeAnalysis")
#' cds_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_coding_exons.bed", package = "SpongeAnalysis")
#' .annotate_introns(intron_file = intron_bed_file,utr3_file = utr3_bed_file ,utr5_file = utr5_bed_file,coding_exon_file = cds_bed_file)
#' @keywords internal
.annotate_introns <- function(intron_file, utr3_file,utr5_file,coding_exon_file){

  x <- rtracklayer::import(intron_file)

  y <- .annotate_regions(x = x ,intron_file = intron_file,utr3_file = utr3_file, utr5_file = utr5_file , coding_exon_file = coding_exon_file)

  y %<>% tibble::as_tibble()
  y %>% dplyr::group_by(type) %>% dplyr::tally()

  # categorize introns into introns (coding + non-coding), intron_utr3 and intron_utr5

  # 1. get non utr introns
  y_grpd <- y %>% dplyr::filter(type != "ucsc_coding_exon") %>%  dplyr::group_by(name) %>% dplyr::add_count()

  type_non_utr_introns <- y_grpd %>%
    dplyr::ungroup() %>%
    dplyr::filter(type == "ucsc_intron" & n ==1 ) %>%
    dplyr::mutate(type = "non_utr_introns")

  type_intron_utr3 <- y_grpd %>%
    dplyr::ungroup() %>%
    dplyr::filter(type == "ucsc_utr3") %>%
    dplyr::mutate(type = "intron_utr3")

  type_intron_utr5 <- y_grpd %>%
    dplyr::ungroup() %>%
    dplyr::filter(type == "ucsc_utr5") %>%
    dplyr::mutate(type = "intron_utr5")

  out <- dplyr::bind_rows(type_non_utr_introns, type_intron_utr3, type_intron_utr5) %>% dplyr::select(-n)

  return( out %>% plyranges::as_granges()  )
}


#' Single function to annotate by utr3, utr5, cds and introns
#'
#' @param x an object of class GRanges
#'
#' @return
#' @export
#'
#' @examples
#'
#' intron_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_introns.bed", package = "SpongeAnalysis")
#'
#' x <- rtracklayer::import(intron_bed_file)
#' utr3_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_3UTR.bed", package = "SpongeAnalysis")
#' utr5_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_5UTR.bed", package = "SpongeAnalysis")
#' cds_bed_file <- system.file("annot_grch38_genecode_v36", "grch38_genecode_v36_coding_exons.bed", package = "SpongeAnalysis")
#'
#' .annotate_regions(x , utr3_file,utr5_file,intron_file,coding_exon_file)
#'
#' @keywords internal
.annotate_regions <- function(x , utr3_file,utr5_file,intron_file,coding_exon_file){

  ucsc_utr3 <- .get_ucsc_utr3(utr3_file = utr3_file)
  ucsc_utr5 <- .get_ucsc_utr5(utr5_file = utr5_file)
  ucsc_introns <- .get_ucsc_introns(intron_file = intron_file)
  ucsc_coding_exon <- .get_ucsc_coding_exon(coding_exon_file = coding_exon_file)

  all_regions <- plyranges::bind_ranges(ucsc_utr3, ucsc_utr5,ucsc_introns, ucsc_coding_exon) %>% dplyr::select(type)

  mapped <- x %>% plyranges::join_overlap_left_within_directed(all_regions) %>%
    tibble::as_tibble()
  gr <- mapped %>% dplyr::distinct() %>% plyranges::as_granges()
  return(gr)

}

