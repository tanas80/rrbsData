
# featuresData class ----
#' Class to store features' properties to separate sample for groups
#' .Data - feature values fo specimens (to plot?)
#' @slot samples markup sample for groups
#' @slot groups mean (avg?) value of feature for groups to know direction of change B-value
#' @slot features significance (rownames - id, p.value, adj.p.value, q.value)
#'                class data.frame for convinient access to columns
#'                data.frame unlike matrix allows columns of different value types
#' @slot adj.p.value_threshold 
#' @slot meta 
setClass(
  # Claas
  'featuresData',
  contains='matrix',
  slots = c(
    'samples' = 'character', 
    'groups' = 'matrix', 
    'features' = 'data.frame',
    'adj.p.value_threshold' = 'numeric',
    'meta' = 'list'
  )
)

#' Constructor object of 'featuresData' class with wilcox.test (U-Mann-Whitney)
#' base code taken from function ftrs_cpgs_get.w.R() of 2015 year
#' @param x matrix (markers x samples) values
#' @param samples_group character vector with group labels of samples
features.w <- function(x, samples_group, meta = list(stat_threshold = 0.05)) {
  
  if (!inherits(x, 'matrix')){
    stop("Argument 'x' must be of class 'matrix'.");
  }
  
  if (length(samples_group) != ncol(x))
    stop("Argument 'ftrs_smpls' length must agree data count of samples");
  
  group_names <- with(list(v=names(sort(table(samples_group)))), { v[v != "" & !is.na(v)]; });
  
  # list of vectors (indicies of specimen)
  group_list_samples <- {
    r <- lapply(group_names, function(group_name) { which(samples_group == group_name);});
    names(r) <- group_names;
    r; }
  if (length(Reduce(intersect, group_list_samples)) > 0)
    stop("Same samples in different groups (it is impossible).")
  
  # For now only two group comparison 
  if (length(group_list_samples) != 2)
    stop("Two groups are mandatory.");
  
  # Calculate mean (avg?) value of feature for every group
  groups <- aaply(
    x, 1,
    function(x_row){
      sapply(group_list_samples, function(group_samples) { mean(x_row[group_samples], na.rm=TRUE); });
    });
  
  system.time(
    features <- ldply(
      1:nrow(x),
      function(i_row, g1, g2, val, rows){
        
        v1 <- g1[i_row,];
        v2 <- g2[i_row,];
        
        p.value = if ((sum(!is.na(v1)) >= 4) & (sum(!is.na(v2)) >= 4))
          tryCatch(wilcox.test(v1,v2)$p.value, error=NA) # U-Mann-Whitney-Wilcoxon
          else NA; 
        
        c(val[i_row,], p.value=p.value);
      },
      g1 = x[, group_list_samples[[1]]],
      g2 = x[, group_list_samples[[2]]],
      val=groups[match(rownames(groups), rownames(x)),],
      rows=rownames(x),
      .parallel=FALSE));
  rownames(features) <- rownames(x);
  features$adj.p.value <- p.adjust(features$p.value, method="fdr");
  
  # ---- q.value ----
  # Calculate qvalue adjusting significance with qvalue package (bioconductor) by Storey
  # Features having calculated p.value
  pv_ixs <- which(!is.na(features$p.value));
  features$q.value <- NA;
  features$q.value[pv_ixs] <- qvalue(features$p.value[pv_ixs])$qvalue;
  
  new (
    "featuresData",
    x,
    samples = samples_group,
    groups = groups,
    features = features,
    meta = meta
  );
}

setGeneric(
  "group_names",
  function(x) standardGeneric("group_names"));

setMethod(
  "group_names",
  signature("featuresData"),
  function(x) {
    with(list(v=names(sort(table(x@samples)))), { v[v != "" & !is.na(v)]; })
  })