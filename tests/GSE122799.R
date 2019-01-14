library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(plyr)
library(rrbsData)
library(biomaRt)

library(doParallel)
registerDoParallel(cores=32)

rdata_path = '~/projs/RData'
cov_path = '~/projs/samples'
samples_url = 'http://www.epigenetic.ru/projects/rrbs-bc/GSE122799.tsv'
xmaI_lib_ranges_url = 'http://tools.epigenetic.ru/LibraRe/hg19/CCCGGG/90-220'

if (!dir.exists(rdata_path))
  stop(sprintf("Path rdata_path='%s' does not exist.", rdata_path))

if (!dir.exists(cov_path))
  stop(sprintf("Path cov_path='%s' does not exist.", cov_path))

tss_gr_rdata_fn <- file.path(rdata_path, 'tss_gr.RData')
if (!file.exists(tss_gr_rdata_fn)){
  # define biomart object
  bioMart_hg19 <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host="grch37.ensembl.org");

  hg19_ann <- with(list(
    ann=getBM(
      attributes = c("hgnc_symbol",
                     "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id",
                     "chromosome_name", "strand", "transcript_start", "transcript_end"),
      mart = bioMart_hg19)),
    cbind(
      ann,
      data.frame(
        chr = ifelse(grepl("^\\d+$", ann$chromosome_name) | ann$chromosome_name == "X",
                     sprintf("chr%s", ann$chromosome_name),
                     ann$chromosome_name),
        tss_start = ifelse(ann$strand > 0, #"+",
                           ann$transcript_start,
                           ann$transcript_end) - 2000,
        tss_end = ifelse(ann$strand > 0,
                         ann$transcript_start,
                         ann$transcript_end) + 2000
      ))
  );

  system.time(
    tss_gr <- GRanges(
      seqnames = hg19_ann$chr,
      ranges = IRanges(hg19_ann$tss_start, hg19_ann$tss_end),
      strand = ifelse(hg19_ann$strand > 0, "+", "-"),
      geneName = hg19_ann$hgnc_symbol,
      transcriptId = hg19_ann$ensembl_transcript_id,
      proteinId = hg19_ann$ensembl_peptide_id
    ));

  save(tss_gr, file=tss_gr_rdata_fn);
} else{
  load(file=tss_gr_rdata_fn, verbose=TRUE)
}

# Get all CpG-pairs of reference genome belonging to XmaI-RRBS library only
cpg_rdata_fn <- file.path(rdata_path, 'cpg.RData')
if (!file.exists(cpg_rdata_fn)) {
  cpg <- subsetByOverlaps(
    load_g_cpg(BSgenome.Hsapiens.UCSC.hg19::Hsapiens), # все CpG пары генома
    {
      # Loading XmaI-RRBS library genomic elements layout
      lib <- import.bed(xmaI_lib_ranges_url)
      # Correcting elements to capture CpG's in restrictase recognition site (overlapping?)
      start(lib) <- start(lib) - 4; end(lib) <- end(lib) + 4;
      lib;
    });
  save(cpg, file=cpg_rdata_fn);
} else {
  load(file=cpg_rdata_fn, verbose=TRUE);
}

rrbs_data_rdata_fn <- file.path(rdata_path, 'rrbs_data.RData')
if (!file.exists(rrbs_data_rdata_fn)){
  gse <- GEOquery::getGEO('GSE122799')[[1]]
  gse_phen <- gse@phenoData@data

  # Select subsample
  gse_phen_subs <- gse_phen[
    grepl('Cell line', gse_phen$characteristics_ch1) | grepl('Normal Breast', gse_phen$characteristics_ch1.1),
    c('title','geo_accession','characteristics_ch1','characteristics_ch1.1','supplementary_file_1') ]
  gse_phen_subs$supplementary_file_1 = as.character(gse_phen_subs$supplementary_file_1)
  print(gse_phen_subs$supplementary_file_1)

  file_df <- ddply(
    gse_phen, c('geo_accession', 'supplementary_file_1'),
    function(row) {
      # row <- gse_phen[4,]
      url <- as.character(row[['supplementary_file_1']])
      dst_fn <- file.path(cov_path, basename(url));

      head_resp <- curlGetHeaders(url)
      url_content_length = with(
        list(grep_list= str_match(head_resp, regex("Content-length: (\\d+)", ignore_case = T))),
        as.numeric(grep_list[!is.na(grep_list[,1]), 2])
      )[1]

      if (!file.exists(dst_fn)) { # || file.size(dst_fn) != url_content_length) {
        cat("Downloading ...", '\n', ' url: ', url, '\n   to:', dst_fn, '\n')
        download.file(url,destfile = dst_fn, method = 'wget')
      } else{
        cat(sprintf("File %s of size %d exists already.\n", dst_fn, url_content_length));
      }

      c('Cov_FileName' = dst_fn, 'Sample_Name' = row[['geo_accession']]) # return
    })

  clinic_df <- read.table(samples_url, quote = '\"', header=TRUE, stringsAsFactors = FALSE)
  sample_df <- merge(x=file_df, y=clinic_df, by.x='geo_accession', by.y='Accession')

  # sample_list ####
  rrbs_samples <- rrbsSampleList(
    {
      sample_list <- lapply(
        1:nrow(sample_df),
        function(i) {
          # i <- 1
          s <- lapply(sample_df, "[", i);
          new("rrbsSampleInfo", s);
        });
      names(sample_list) <- lapply(sample_list, '[[', "geo_accession")

      sample_list;
    },
    groups = list(
      "EType" = c(
        "NORM" = "green",
        "CL" = "black",
        "LumA" = "blue",
        "LumB" = "cyan",
        "HER2" = "magenta",
        "TN" = "red",
        "UN" = "lightgray"
      ),
      "MTypeHE" = c(
        "CL" = "black",
        "hiHER" = "darkmagenta",
        "hiHET" = "gray",
        "hiLumAB" = "darkturquoise",
        "NORM" = "green",
        "modLumB" = "cyan",
        "modTNBC" = "red",
        "modHER" = "magenta"
      )));

  rrbs_data <- loadRrbsData(rrbs_samples, cpg, "Breast cancer")
  save(rrbs_data, file=rrbs_data_rdata_fn)
} else {
  load(file=rrbs_data_rdata_fn, verbose=TRUE)
}

m_rdata_path <- file.path(rdata_path, 'm_flt,m_he.RData')
if (!file.exists(m_rdata_path)) {
  # Filtering CpG in each specimen for minimal coverage depth and samples for minimal coverage width
  rrbs_data_flt <- filterRrbsData(rrbs_data, depth_threshold = 10, width_threshold = 10000)

  # Create object of 'rrbsMatrix' class with minimal quality filtered data
  rrbs_m <- rrbsMatrix(rrbs_data_flt, tss_gr);
  dim(rrbs_m)

  # удалим образцы, которые имеют мало общих цитозинов с другими образцами
  m_dj <- filter_disjoint(rrbs_m, min_shared_cpg = 10000);

  # удалим цитозины, чтобы в общем наборе данных было не более 20% пропущенных
  m_flt <- filter_miss(m_dj, max_miss = 0.2);

  ent_c <- apply(m_flt, 1, function(v) { entropy(v[!is.na(v)]);} );
  m_he <- with(list(
    x=m_flt[ent_c >= 3.5 & !is.na(ent_c),]),
    new(
      "rrbsMatrix",
      x@.Data,
      coverage = x@coverage,
      samples = x@samples,
      coords = x@coords,
      genes = x@genes,
      description=sprintf("%s,\n CpG entropy >= 3.5", description(x))));

  save(m_flt, m_he, file=m_rdata_path)
} else {
  load(file=m_rdata_path, verbose=TRUE)
}


