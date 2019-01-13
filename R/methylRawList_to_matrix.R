#' @import plyr
#' @include merge_pairwise.R

methylRawList_to_matrix <- function(x, .parallel=TRUE )
{
  if (
    is(x,"methylRawList") ||
    (is(x, "list") && (sum(sapply(x, function(item) { !is(item, "methylRaw") })) == 0))
  ) {

    # x_dl - List of data.frames of samples' methylation B-value by coord defined as "chr1:10000-10001"
    x_dl <- with(
      list(d=x),
      {
        dl <- llply(
          1:(length(d)),
          function(i, d)
          {
            sample_name <- d[[i]]@sample.id;
            mr <- getData(d[[i]]); # mr - данные по одному образцу methylRaw

            cat(sprintf("Sample %d - %s\n", i, sample_name));

            coord <- if (nrow(mr) > 0) paste(mr$chr, ":",mr$start, "-",mr$end, "/", mr$strand, sep="") else rep(1,0);

            df <- data.frame(
              coord = coord,
              B.val = mr$numCs/mr$coverage,
              cvrg = mr$coverage
            );
            names(df) <- c(
              "coord",
              paste(sample_name, "bval", sep="."),
              paste(sample_name, "cvrg", sep="."));
            df;
          },
          d,.parallel=.parallel);
        names(dl) <- sapply(d, attr, "sample.id");
        dl;
      });

    coord_all <- data.frame(
      coord=unique(do.call(c, lapply(x_dl, function(d){ as.character(d$coord);})))
    );
    cat(sprintf("unique all coord: %d\n", nrow(coord_all)));

    x_dl_exp <- llply(
      x_dl,
      function(sample_data){
        # expand sample data for missed coords
        d_missed <- coord_all$coord[!(coord_all$coord %in% sample_data$coord)];
        d_supp <- data.frame(coord=d_missed, rep(NA, length(d_missed)), rep(NA, length(d_missed)));
        names(d_supp) <- names(sample_data);
        d_exp <- rbind(sample_data, d_supp);
        d_exp;
      },
      .parallel=.parallel);

    x_dm_bval <- merge_pairwise(
      lapply(1:length(x_dl_exp), function(i) { r <- x_dl_exp[[i]][, c(1,2)]; colnames(r) <- c("coord" , x[[i]]@sample.id); r; }),
      merge_fun=function(d1,d2)
      {
        # cat("merge: \nd1(",names(d1),"), \nd2(", names(d2),")\n");
        merge(d1, d2, by.x="coord", by.y="coord");
      },
      .parallel=.parallel); # параллельность помогает - проверено

    x_dm_cvrg <- merge_pairwise(
      lapply(1:length(x_dl_exp), function(i) { r <- x_dl_exp[[i]][, c(1,3)]; colnames(r) <- c("coord" , x[[i]]@sample.id); r; }),
      merge_fun=function(d1,d2)
      {
        # cat("merge: \nd1(",names(d1),"), \nd2(", names(d2),")\n");
        merge(d1, d2, by.x="coord", by.y="coord");
      },
      .parallel=.parallel); # параллельность помогает - проверено

    list(
      b.value = {
        m <- as.matrix(x_dm_bval[, names(x_dm_bval) != "coord"])
        mode(m) <- "numeric";
        rownames(m) <- x_dm_bval$coord;
        m;
      },
      coverage = {
        m <- as.matrix(x_dm_cvrg[, names(x_dm_cvrg) != "coord"])
        mode(m) <- "numeric";
        rownames(m) <- x_dm_cvrg$coord;
        m;
      });
  }
  else
    stop("Invalid class of x (methylRawList or list of methylRaw objs expected).");
}
