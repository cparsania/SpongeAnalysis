# IR analysis functions

retained_introns_files  <-
  list.files(path  = "000_data_files/ir_finder_s_output/",
             pattern =  "*nondir_SRR*", full.names = T,recursive = T)

names(retained_introns_files) <- stringr::str_extract(retained_introns_files,
                                                      pattern = "SRR([^_])+")

# read irfinders output

#' Read IR-finderS output
#'
#' @param files a character string denoting irfinderS output file ending with *IR-nondir.
#' @param add_prefix_chr
#' @param remove_prefix_chr
#' @param select_columns
#'
#' @return
#' @export
#'
#' @examples
read_irfinderS_output <- function(files, add_prefix_chr = TRUE, remove_prefix_chr = FALSE, select_columns = TRUE){

  # check if each values in files are named, not null or not NA.
  stopifnot("files must be named vector" = !(is.null(files) %>% all()) | !(is.na(files) %>% any()))

  stopifnot("'add_prefix_chr' must be logical" = is.logical(add_prefix_chr))

  stopifnot("'remove_prefix_chr' must be logical" = is.logical(remove_prefix_chr))


  retained_introns_list <- purrr::map(retained_introns_files, ~ ..1 %>%

                                        # read
                                        readr::read_delim(delim = "\t",show_col_types = FALSE) %>%

                                        # rename to lower case
                                        dplyr::rename_all(~tolower(.)))

  # add prefix chr
  if(add_prefix_chr){
    retained_introns_list <- retained_introns_list %>%
      purrr::map( ~ ..1 %>% dplyr::mutate(chr = stringr::str_c("chr", chr, sep ="")))
  }


  # remove prefix chr
  if(remove_prefix_chr){
    retained_introns_list <- retained_introns_list %>%
      purrr::map( ~ ..1 %>% dplyr::mutate(chr = stringr::str_replace(chr, pattern = "chr",replacement = "")))
  }

  # add class
  class(retained_introns_list) <- c("irfinders_output", class(retained_introns_list))

  # select columns
  if(select_columns){
    retained_introns_list <- select_cols_irfinderS_output(retained_introns_list)
  }

  return(retained_introns_list)
}

# select columns from irfinder-s output

select_cols_irfinderS_output <- function(x, keep_columns = c("chr",
                                                             "start",
                                                             "end",
                                                             "name",
                                                             "null",
                                                             "strand",
                                                             "coverage",
                                                             "introndepth",
                                                             "spliceleft",
                                                             "spliceright",
                                                             "spliceexact",
                                                             "irratio",
                                                             "warnings")){

  .validate_irfinders_object(x)

  x <- purrr::map(x , ~..1 %>% dplyr::select(keep_columns))

  return(x)


}


# get intron master list
get_intron_master_list <- function(x, add_meta_data = T, bs_genome_object = NULL){

  .validate_irfinders_object(x)
  stopifnot("'add_meta_data' must be logical." = is.logical(add_meta_data))

  # if x has more than one elements, select first to prepare intron list.
  # As all elements have same rows taking any one to prepare intron list is ok.

  x <- x[[1]]

  x <- x  %>%

    # select cols
    dplyr::select(chr, start, end, name, null, strand) %>%

    # prepare gene symbol and gene id
    dplyr::mutate(gene_name = stringr::str_replace(name,"/.*" ,""),
                  gene_id = stringr::str_match(string = name,pattern = "/(.*)/")[,2],
                  intron_type = stringr::str_replace(name,".*/" ,"") ) %>%

    # prepare intron id
    dplyr::mutate(intron_id = stringr::str_c("intron_",dplyr::row_number()))

  # add meta data
  if(add_meta_data){
    stopifnot("when add_meta_data = T, bs_genome_object must be an object of class BSgenome" = is(bs_genome_object, "BSgenome"))

    x <- x %>% dplyr::rename("seqnames" = "chr")  %>% plyranges::as_granges()

    x <- x %>%
      # add sequence for each intron
      dplyr::mutate(seq =  BSgenome::getSeq(bs_genome_object,x)) %>%

      # add GC for each intron
      dplyr::mutate(GC =  BSgenome::letterFrequency(seq, letters = "GC", as.prob = T) %>%  as.numeric()) %>%

      # add length for each intron
      dplyr::mutate(length =  BSgenome::width(seq))

    # convert back into tibble
    x <- x %>% tibble::as_tibble()

  }


  return(x)

}


# add intron metadata to irfinder s output

map_intron_meta_data <- function(x , bs_genome_object = BSgenome.Hsapiens.UCSC.hg38){
  .validate_irfinders_object(x)
  meta_data <- get_intron_master_list(x, add_meta_data = T, bs_genome_object = bs_genome_object) %>% dplyr::select(8:dplyr::last_col())
  x_mapped <- purrr::map(x, ~ ..1 %>% dplyr::bind_cols(meta_data))
  class(x_mapped) <- class(x)
  return(x_mapped)
}

# filter IR results

filter_irfinderS_output <- function(x , min_intron_cov = 0.95,
                                    min_intron_depth = 5,
                                    minimum_splice_exact = 5,
                                    min_irratio = 0.0001){
  .validate_irfinders_object(x)

  apply_intron_coverage_cutoff = T

  # filter by coverage
  if(apply_intron_coverage_cutoff) {
    x_filt <- x %>%
      purrr::map(~ ..1 %>% dplyr::filter(coverage >= min_intron_cov ))
  }


  # filter by depth
  apply_depth_cutoff <- TRUE

  if(apply_depth_cutoff){
    x_filt <- x_filt %>%
      purrr::map(~ ..1 %>% dplyr::filter(introndepth >= min_intron_depth ))
  }


  # filter by splice exact
  filter_by_splice_exact <- TRUE

  if(filter_by_splice_exact){
    x_filt <- x_filt %>%
      purrr::map(~ ..1 %>% dplyr::filter(spliceexact >= minimum_splice_exact))
  }

  # filter by IR ratio
  apply_ir_ratio_cuoff <- TRUE
  if(apply_ir_ratio_cuoff){
    x_filt <- x_filt %>%
      purrr::map(~ ..1 %>%
                   dplyr::filter(irratio >= min_irratio ))
  }

  class(x_filt) <- class(x)
  return(x_filt)

}



# internals

.validate_irfinders_object <- function(x){
  stopifnot("x must be an object of class irfinders_output" = is(x, "irfinders_output") )
}



