# IR analysis functions

retained_introns_files  <-
  list.files(path  = "000_data_files/ir_finder_s_output/",
             pattern =  "*nondir_SRR*", full.names = T,recursive = T)

names(retained_introns_files) <- stringr::str_extract(retained_introns_files,
                                                      pattern = "SRR([^_])+")

# read irfinders output

#' Read IR-finderS output
#'
#' @param files a character vector denoting irfinderS output file(s) ending with suffix "IR-nondir".
#' @param add_prefix_chr logical, whether to add prefix 'chr' in the column seqnames
#' @param remove_prefix_chr logical, whether to remove prefix 'chr' from the column seqnames
#' @param select_columns logical, whether to subset columns. If TRUE below columns will be subsetted
#'  + chr
#'  + start
#'  + end
#'  + name
#'  + null
#'  + strand
#'  + coverage
#'  + introndepth
#'  + spliceleft
#'  + spliceright
#'  + spliceexact
#'  + irratio
#'  + warnings
#' @return
#' @export
#'
#' @examples
#' example_dir <- "~/Documents/Projects/15_SpongeAnalysisRpkg/SpongeAnalysis/inst/extdata"
#' example_files <- fs::dir_ls(example_dir, regexp = "exp*")
#' names(example_files) <- c("exp1" , "exp2")
#' read_irfinderS_output(files = example_files,  add_prefix_chr = F)
read_irfinderS_output <- function(files,
                                  add_prefix_chr = TRUE,
                                  remove_prefix_chr = FALSE){

  # check if each values in files are named, not null or not NA.
  stopifnot("files must be named vector" = !(is.null(files) %>% all()) | !(is.na(files) %>% any()))

  stopifnot("'add_prefix_chr' must be logical" = is.logical(add_prefix_chr))

  stopifnot("'remove_prefix_chr' must be logical" = is.logical(remove_prefix_chr))


  retained_introns_list <- purrr::map(files, ~ ..1 %>%

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
  retained_introns_list <- .assign_class_spongeAnalysis(retained_introns_list)

  return(retained_introns_list)
}

# select columns from irfinder-s output
#' Subset columns from IRfinder-S output
#'
#' @param x an object of class spongeAnalysis
#' @param keep_columns a character vector denoting column names to keep in the output dataframe
#'  + chr
#'  + start
#'  + end
#'  + name
#'  + null
#'  + strand
#'  + coverage
#'  + introndepth
#'  + spliceleft
#'  + spliceright
#'  + spliceexact
#'  + irratio
#'  + warnings
#' @return
#' @export
#'
#' @examples
#' example_dir <- "~/Documents/Projects/15_SpongeAnalysisRpkg/SpongeAnalysis/inst/extdata"
#' example_files <- fs::dir_ls(example_dir, regexp = "exp*")
#' names(example_files) <- c("exp1" , "exp2")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = F)
#' select_cols_irfinderS_output(x)
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

  class(x) <- c("spongeAnalysis", class(x))
  return(x)


}


#' Prepare a master list of introns
#'
#' @param f a character string denoting irfinderS output file ending with suffix "IR-nondir".
#' @param add_meta_data logical, whether to map meta data or not
#' @param bs_genome_object an object of class BSgenome
#' @param add_prefix_chr logical, whether to add prefix 'chr' in the column seqnames
#' @param remove_prefix_chr logical, whether to remove prefix 'chr' from the column seqnames
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' example_dir <- "~/Documents/Projects/15_SpongeAnalysisRpkg/SpongeAnalysis/inst/extdata"
#' example_files <- fs::dir_ls(example_dir, regexp = "exp*")
#' get_intron_master_list(f = example_files[1],add_prefix_chr = TRUE,bs_genome_object = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
get_intron_master_list <- function(f, add_meta_data = T,
                                   add_prefix_chr = FALSE,
                                   remove_prefix_chr = FALSE,
                                   bs_genome_object = NULL){


  stopifnot("'add_meta_data' must be logical." = is.logical(add_meta_data))

  names(f) <- "file"
  x <- read_irfinderS_output(files = f,
                        add_prefix_chr = add_prefix_chr,
                        remove_prefix_chr = remove_prefix_chr)

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


# add intron metadata to irfinder-s output

#' Map metadata irfinder-s output
#' @description This function allows mapping DNA sequence, GC content and intron length to irfinder-s output
#' @param x an object of class spongeAnalysis.
#' @param bs_genome_object an object of class BSgenome
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' example_dir <- "~/Documents/Projects/15_SpongeAnalysisRpkg/SpongeAnalysis/inst/extdata"
#' example_files <- fs::dir_ls(example_dir, regexp = "exp*")
#' names(example_files) <- c("exp1" , "exp2")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = T, select_columns = T)
#' map_intron_meta_data(x = x,  bs_genome_object = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
map_intron_meta_data <- function(x , bs_genome_object = BSgenome.Hsapiens.UCSC.hg38){
  .validate_irfinders_object(x)
  meta_data <- get_intron_master_list(x, add_meta_data = T, bs_genome_object = bs_genome_object) %>% dplyr::select(8:dplyr::last_col())
  x_mapped <- purrr::map(x, ~ ..1 %>% dplyr::bind_cols(meta_data))
  class(x_mapped) <- class(x)
  return(x_mapped)
}

# filter IR results

#' Filter irfinder-s results
#'
#' @param x an object of class sponge analysis
#' @param min_intron_cov a numeric value between 0 and 1, default 0.95, denoting minimum value for coverage cutoff
#' @param min_intron_depth a numeric value, default 5, denoting average sequence depth for each intron.
#' @param minimum_splice_exact a numeric value, default 5, denoting number of reads supporting intron splicing.
#' @param min_irratio a numeric value, default 0.0001, denoting ir-ratio.
#'
#' @return an object of class spongeAnalysis
#' @export
#'
#' @examples
#' example_dir <- "~/Documents/Projects/15_SpongeAnalysisRpkg/SpongeAnalysis/inst/extdata"
#' example_files <- fs::dir_ls(example_dir, regexp = "exp*")
#' names(example_files) <- c("exp1" , "exp2")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = T, select_columns = T)
#' filter_irfinderS_output(x)
#'
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

  return(x_filt)

}





#' Check if the object belongs to class spongeAnalysis.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @keywords internal
.validate_irfinders_object <- function(x){

  stopifnot("x must be an object of class spongeAnalysis" = is(x, "spongeAnalysis") )
}


#' Assign class spongeAnalysis
#' @description This function does all mandatory checks before it assigns class spongeAnalysis.
#' @param x a list or dataframe to which class spongeAnalysis to assign.
#'
#' @return an object of class spongeAnalysis.
#' @export
#' @keywords internal
.assign_class_spongeAnalysis <- function(x){

  # x can be a dataframe or list of dataframes
  # if dataframe it must have mandatory columns
  # if a list it must have mandatory columns in each dataframe and same number of rows in each dataframe

  if(is(x , "data.frame")){
    .check_mendate_columns(x)
    x <- list(x)
    class(x) <- c("spongeAnalysis", class(x))
  }

  if(is(x , "list")){
    purrr::walk(x, ~ .check_mendate_columns(..1))

    # all elems of the list have same number of rows
    x_nrows <- purrr::map_int(x, ~..1 %>% nrow())

    if(!all(x_nrows==x_nrows[1])){
      stop("All elements in x must have same number of rows.")
    }

    class(x) <- c("spongeAnalysis", class(x))
  }

  return(x)

}


#' Check mandatory columns in a dataframe
#' @description This function helps to create mandatory column in a dataframe
#' @param x dataframe
#' @param mandat_columns a character vector denoting mandatory columns in x.
#'
#' @return TRUE or ERROR
#' @export
#' @keywords internal
.check_mendate_columns <- function(x,
                                   mandat_columns = c("chr",
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

  stopifnot("x must be a dataframe" = is(x , "data.frame"))

  x_cols <- colnames(x)

  col_not_found <- mandat_columns[is.na(match(mandat_columns, x_cols))]
  if(length(col_not_found) > 0){
    stop(cli::format_error(c("x" ="Column{?s} {col_not_found} not found")))
  }

  return(TRUE)

}





