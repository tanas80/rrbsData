merge_pairwise <- function(dl, merge_fun, .parallel=FALSE){

  cat("\nmerge_pairwise:\n");
  cat(paste(sapply(dl, function(d) { paste(names(d), collapse=", ");}), collapse="\n"), "\n");

  if (length(dl)==1){
    dl[[1]];
  }
  else if (length(dl)==2) {
    merge_fun(dl[[1]], dl[[2]]);
  }
  else {
    pairwise_dl <- llply(
      1:(length(dl)/2),
      function(i){
        merge_fun(dl[[i*2-1]], dl[[i*2]])
      },
      .parallel=.parallel);

    #cat(sapply(pairwise_dl, names));

    merge_pairwise(
      if ((length(dl) %% 2) == 0)
        pairwise_dl
      else
        c(pairwise_dl, list(dl[[length(dl)]])),
      merge_fun,
      .parallel=.parallel);
  }
};
