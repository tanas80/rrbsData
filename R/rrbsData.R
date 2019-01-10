#' @import plyr
#' @import dplyr
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import doParallel
#' @import methylKit
#' @import Biostrings
#' @importFrom grDevices rainbow
#' @importFrom graphics legend plot
#' @importFrom methods is new slotNames
#' @importFrom stats as.dendrogram cor dendrapply dist fisher.test hclust is.leaf p.adjust quantile sd
#' @importFrom utils read.table


registerDoParallel();


####Support Functions####
dist_pairwise <- function(x, method="manhattan"){

  if (!is.element(method, c("manhattan", "pearson")))
    stop(sprintf("Invalid distance method '%s", method));

  res <- if (
    is(x, "methylRawList") ||
    (is(x, "list") && sum(sapply(x, function(item) { !is(item, "methylRaw")})) == 0) ||
    is(x, "rrbsDataList")
  ) {
    ix <- 1:length(x);
    n <- length(ix);

    dists <- ldply(
      1:(n-1),
      function(i)
      {
        # построим список каждый-с-каждым (но только один раз)
        ddply(
          data.frame(
            "i"=rep(i, n-i),
            "j"=(i+1):n
          ),
          c("i","j"),
          function(row){
            s1 <- x[[row[1,1]]];
            s2 <- x[[row[1,2]]];

            pd <- new("methylRawList",list(s1,s2));
            pm <- unite(pd, destrand=FALSE);

            if (method == "pearson") {
              1 - (1 + cor(pm$numCs1/pm$coverage1, pm$numCs2/pm$coverage2, method="pearson"))/2;
            } else if (method == "manhattan") {
              sum(abs(pm$numCs1/pm$coverage1 - pm$numCs2/pm$coverage2))/nrow(pm);
            } else {
              stop("Distance method error");
            };
          });
      },
      .parallel=TRUE);

    list(n=n, dists=dists, labels=sapply(x, attr, "sample.id"));
  } else if (class(x) == "matrix") {
    m <- x;

    ix <- 1:ncol(m);
    n <- length(ix);

    dists <- ldply(
      1:(n-1),
      function(i)
      {
        # построим список каждый-с-каждым (но только один раз)
        ddply(
          data.frame(
            "i"=rep(i, n-i),
            "j"=(i+1):n
          ),
          c("i","j"),
          function(row){
            v1 <- m[,row[1,1]];
            v2 <- m[,row[1,2]];

            if (method == "pearson") {
              1- (1 + cor(v1, v2, use="pairwise.complete.obs", method="pearson"))/2;
            } else if (method == "manhattan") {
              common_i <- !is.na(v1) & !is.na(v2);
              sum(abs(v1[common_i] - v2[common_i]))/sum(common_i);
            } else {
              stop("Distance method error");
            };
          });
      },
      .parallel=TRUE);

    list(n=n, dists=dists, labels=colnames(m));
  } else {
    stop("Invalid class of x ('matrix' or 'methylRawList' or 'list of methylRaw' expected).");
  };

  n <- as.numeric(res$n); dists <- res$dists;

  l <- (n/2)*(n-1);
  if (nrow(dists) != l)
    stop("Invalid operation");

  names(dists) <- c("i", "j", "dist");

  structure(
    as.numeric(dists$dist),
    Size = n, Labels = res$labels,
    Diag = FALSE, Upper = FALSE,
    method = "user", class = "dist");
};

.dist.cor=function(x, method="pearson", abs=TRUE, diag=FALSE, upper=FALSE)
{
  if (!is.na(pmatch(method, "pearson")))
    method <- "pearson"
  METHODS <- c("pearson")
  method <- pmatch(method, METHODS)
  if (is.na(method))
    stop("invalid correlation method")
  if (method == -1)
    stop("ambiguous correlation method")

  xcor = cor(t(x), method=METHODS[method])
  if(abs)
    xcor = 1-abs(xcor)
  else
    xcor = 1-xcor
  if(upper)
    d <- xcor[upper.tri(xcor,diag=diag)]
  else
    d <- xcor[lower.tri(xcor,diag=diag)]
  attr(d, "Size") <- nrow(x)
  attr(d, "Labels") <- dimnames(x)[[1L]]
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- METHODS[method]
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
}

.cluster=function(x, dist.method="correlation", hclust.method="ward", plot=TRUE,
                  treatment=treatment,sample.ids=sample.ids,context){
  DIST.METHODS <- c("correlation", "euclidean", "maximum", "manhattan", "canberra",
                    "binary", "minkowski")
  dist.method <- pmatch(dist.method, DIST.METHODS)

  HCLUST.METHODS <- c("ward", "single", "complete", "average", "mcquitty",
                      "median", "centroid")
  hclust.method <- pmatch(hclust.method, HCLUST.METHODS)
  if (is.na(hclust.method))
    stop("invalid clustering method")
  if (hclust.method == -1)
    stop("ambiguous clustering method")

  if(DIST.METHODS[dist.method] == "correlation")
    d = .dist.cor(t(x))
  else
    d=dist(scale(t(x)), method=DIST.METHODS[dist.method]);

  hc=hclust(d, HCLUST.METHODS[hclust.method]);

  if(plot){
    treatment=treatment
    sample.ids=sample.ids

    my.cols=rainbow(length(unique(treatment)), start=1, end=0.6)

    col.list=as.list(my.cols[treatment+1])
    names(col.list)=sample.ids

    colLab <- function(n, col.list)
    {
      if(is.leaf(n))
      {
        a <- attributes(n)

        attr(n, "nodePar") <- c(a$nodePar, list(
          lab.col = col.list[[a$label]],
          lab.cex=1,
          col=col.list[[a$label]],
          cex=1,
          pch=16 ));
      }

      n;
    }

    dend = as.dendrogram(hc)
    dend_colored <- dendrapply(dend, colLab, col.list)

    plot(dend_colored, main = paste(context, "methylation clustering"),
         sub = paste("Distance method: \"", DIST.METHODS[dist.method],
                     "\"; Clustering method: \"", HCLUST.METHODS[hclust.method],"\"",sep=""),
         xlab = "Samples", ylab = "Height"
    );
    # end of plot
  }
  return(hc)
}

rowSds <- function(x, center=NULL, ...) {
  n <- !is.na(x);
  n <- rowSums(n);
  n[n <= 1] <- NA;

  if (is.null(center)) {
    center <- rowMeans(x, ...);
  }

  x <- x - center;
  x <- x*x;
  x <- rowSums(x, ...);
  x <- x/(n-1);

  sqrt(x);
}


clusterSamples <- function(.Object, dist="correlation", method="ward",
                           sd.filter=TRUE,sd.threshold=0.5,
                           filterByQuantile=TRUE, plot=TRUE)
{
  mat      =getData(.Object)
  # remove rows containing NA values, they might be introduced at unite step
  mat      =mat[ rowSums(is.na(mat))==0, ]

  meth.mat = mat[, .Object@numCs.index]/
    (mat[,.Object@numCs.index] + mat[,.Object@numTs.index] )
  names(meth.mat)=.Object@sample.ids

  # if Std. Dev. filter is on remove rows with low variation
  if(sd.filter){
    if(filterByQuantile){
      sds=rowSds(as.matrix(meth.mat))
      cutoff=quantile(sds,sd.threshold)
      meth.mat=meth.mat[sds>cutoff,]
    }else{
      meth.mat=meth.mat[rowSds(as.matrix(meth.mat))>sd.threshold,]
    }
  }

  .cluster(meth.mat, dist.method=dist, hclust.method=method,
           plot=plot, treatment=.Object@treatment,
           sample.ids=.Object@sample.ids,
           context=.Object@context)

}

treatment <- c (0,1,10,10,0,1,11);

#' plot hierarhical clustering
#' @param hc hclust object
#' @param types character vector with names by hc$labels
#' @param colors character vector with names by types
my.plot.hclust <- function(hc, types, colors, main="CpG methylation clustering", cex=1.0, pch=20)
{
  label_colors <- {
    v <- as.character(colors[types[hc$labels]]);
    names(v) <- hc$labels;
    v;
  };

  colLab <- function(n, label_colors)
  {
    if( is.leaf(n))
    {
      a <- attributes(n);

      color <- label_colors[a$label];

      attr(n, "nodePar") <- c(
        a$nodePar,
        list(
          lab.col=color,
          lab.cex=cex,
          col = color,
          cex=cex,
          pch=pch
        )
      );
    }

    n;
  };

  dend <- as.dendrogram(hc);
  dend_colored <- dendrapply(dend, colLab, label_colors);

  plot(dend_colored, main=main,
       #sub=paste("Distance method: \"", hc$dist.method, "\"", sep=""),
       #xlab="Samples",
       ylab="Distance"
  );
  if (!is.null(colors)) {
    legend("topright",
           legend=names(colors),
           fill = as.character(colors),
           col = as.character(colors),
           text.col = as.character(colors),
           cex=cex);
  }
}

load_g_cpg <- function(g){
  do.call(c,
          llply(names(g), function(chr) {
            chr_ranges = matchPDict(PDict("CG"), g[[chr]])[[1]];

            #d <- data.frame(
            #  chr = rep(chr, length(chr_ranges)),
            #  start = start(chr_ranges), data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAkAAAAJCAYAAADgkQYQAAAAMElEQVR42mNgIAXY2Nj8x8cHC8AwMl9XVxe3QqwKcJmIVwFWhehW4LQSXQCnm3ABAHD6MDrmRgfrAAAAAElFTkSuQmCC
            #  end = end(chr_ranges),
            #  strand = rep("*", length(chr_ranges)),
            #  stringsAsFactors = FALSE
            #  );

            GRanges(seqnames=chr, ranges=chr_ranges, seqinfo=seqinfo(g));
          }, .parallel=TRUE));
};

load_g_sites <- function(g, site) {
  do.call(c,
          llply(
            names(g),
            function(chr){
              chr_ranges = matchPDict(PDict(site), g[[chr]])[[1]];

              GRanges(seqnames=rep(chr, length(chr_ranges)), ranges=chr_ranges, seqinfo = seqinfo(g));
            },
            .parallel = TRUE));
};

parse_sample_name <- function(rrbs_sample_name)
{
  parts <- strsplit(x=rrbs_sample_name, split="-", fixed=TRUE)[[1]];
  lab <- parts[1];
  year <- parts[2];
  project <- parts[3];
  sample <- parts[4];
  tissue <- parts[5];
  extra <- if (length(parts) > 5) parts[6] else "";

  c("lab"=lab,
    "year"=as.numeric(year),
    "project"=project,
    "sample"=sample,
    "tissue"=tissue,
    "extra"=extra);
}

####Main Functions####


rrbsSampleInfo <- setClass(
  "rrbsSampleInfo",
  contains = "list"
  );

setMethod(
  "show",
  signature = c("rrbsSampleInfo"),
  function(object){
    for (prop_name in names(object)){
      cat(paste("$", prop_name, "\n",sep=""));
      print(object[[prop_name]]);
      cat("\n");
    }
  });

rrbsSampleList <- setClass(
  "rrbsSampleList",
  contains = "list",
  slots = c(
    groups = "list" # list of "data.frame" objects
    ));

setMethod(
  "as.data.frame",
  signature =c("rrbsSampleList"),
  function(x, row.names=NULL){
    cols <- unique(lapply(x, names));
    ldply(
      1:length(x),
      function(i){
        data.frame(x[[i]]);
      }
    )
  });

setGeneric(
  "getfield",
  function(x, field) standardGeneric("getfield"))

# getfield <- function(x, field) {
#   sapply(x, function(smpl){ smpl[field] })
# }

setMethod(
  'getfield',
  signature = c("rrbsSampleList"),
  function(x, field) {
    sapply(x, function(smpl){ smpl[field] })
  }
)

setMethod(
  "show",
  signature = c("rrbsSampleList"),
  function(object){
    show(names(object));
  });

`[.rrbsSampleList` <- function(x, i){
  #cat('[.rrbsSampleList\n');
  #cat("i =", i, "\n");
  #cat("class(i) =", class(i), "\n");

  new(
    "rrbsSampleList",
    {
      indecies <- c();

      if (is.logical(i))
      {
        indicies <- (1:length(x))[i];
      } else if (is.numeric(i)) {
        indicies <- i;
      } else if (is.character(i)) {
        indicies <- match(i, names(x));
      }
      # cat("indicies = ", indicies, "\n");
      ss <- lapply(indicies, function(li) { x[[li]];});
      names(ss) <- lapply(indicies, function(li) { names(x)[[li]];});
      ss;
    },
    groups = x@groups
  );
};

`$.rrbsSampleList` <- function(x, attr_name){
  #cat("attr_name:", attr_name, "\n");
  sapply(x, "[[", attr_name);
};


loadRrbsSamples <- function(filename){

  # Ошибка "line NNN did not have NNN elements" бывает из-за того, что на листе были ошибки формул
  t <- read.table(filename, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", dec=".");

  readRrbsSamples(t);
}

readRrbsSamples <- function(tbl, lis){
  # tbl - данные из гуглотаблицы
  # lis - данные из БД ЛИС

  old_path = "/media/DATA-EPI/PROJECTS/RRBS2/samples"

  # из таблицы удалим строки, пустые в столбце Sample_Name
  tbl <- tbl[ !is.na(tbl$Sample_Name) & tbl$Sample_Name != "",]
  #rrbs_samples <- rrbs_samples[rrbs_samples$Chip %in% c("epi-2015-0012"), ]

  #rrbs_samples <- rrbs_samples[ rrbs_samples$EType=="NORM",]

  tbl$MethType <- ifelse(
    tbl$MethType %in%
      c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10"),
    tbl$MethType, NA);

  tbl$Sample_Name_Short <- ifelse(
    !is.na(tbl$Sample_Name_Short),
    tbl$Sample_Name_Short,
    sapply( tbl$Sample_Name, function(sn) { p <- parse_sample_name(sn); paste(p["sample"], p["extra"], sep="") })
  );

  tbl$Src <- "epi";

  tbl$Cov_FileName <- with(
    list(m=match(tbl$Sample_Name, lis$name)), {
      ifelse(
        !is.na(m),
        { # ЕСТЬ данные в ЛИС
          # lis$bam_filename

          # использовать папку данных направления (если в БД ЛИС станет правильно)
          file.path(lis$filepath, gsub("\\.bam$", ".bismark.cov.gz", lis$bam_filename))[m]
          #gsub("\\.bam$", ".bismark.cov.gz", lis$bam_filename)[m]
        },
        { # НЕТ данных в ЛИС

          file.path(
            old_path, tbl$Sample_Name, "/",
            paste(tbl$Sample_Name, ".bismark_bt2_cln.sorted.bismark.cov.gz", sep=""))
        }
      )
    }
  )

  # tbl$Cov_FileName <- paste(
  #   "../samples/", tbl$Sample_Name, "/",
  #   tbl$Sample_Name, ".bismark_bt2_cln.sorted.bismark.cov", sep="");

  # tbl$Bed_FileName <- paste(
  #   "../samples/", tbl$Sample_Name, "/",
  #   tbl$Sample_Name, ".bismark_bt2_cln.sorted.bam.bed", sep=""
  # );

  new(
    "rrbsSampleList", {
      sample_list <- lapply(
        1:nrow(tbl),
        function(i) {
          s <- lapply(tbl, "[", i);
          new("rrbsSampleInfo", s);
        });
      names(sample_list) <- lapply(sample_list, '[[', "Sample_Name_Short")
      sample_list;
    },
    #     groups = list(
    #       "EType" = data.frame(
    #         name=c("NORM", "CL", "LumA", "LumB", "HER2", "TN", "UN"),
    #         id=c(10, 11, 12, 13, 14, 15, 99),
    #         color=c("green", "black", "blue", "cyan", "magenta", "red", "lightgray")
    #       ),
    #       "GType" = data.frame(
    #         name=c("NORM", "CL", "duct", "lobul", "mix", "medul", "meta", "mpc", "neuro", "UN"),
    #         id=c(20, 21, 22, 23, 24, 25, 26, 27, 28, 99),
    #         color=c("green", "black", "blue", "red", "yellow", "gray", "gray", "gray", "gray", "lightgray")
    #       )));
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
      "GType" = c(
        "NORM" = "green",
        "CL" = "black",
        "duct" = "black",
        "lobul" = "red",
        "mix" = "yellow",
        "medul" = "gray",
        "meta" = "gray",
        "mpc" = "gray",
        "neuro" = "gray",
        "UN" = "lightgray"
      ),
      "MTypeHE" = c(
        "CL" = "black",
        "NORM" = "green",
        "modTNBC" = "red",
        "hiHET" = "gray",
        "modHER" = "magenta",
        "hiHER" = "darkmagenta",
        "hiLumAB" = "darkturquoise",
        "modLumB" = "cyan"
      ),
      "ACT_EffA" =c(
        "Full Regression" = "green",
        "Partial Regression" = "darkgreen",
        "Stabilization" = "orange",
        "Progression" = "red",
        "Residual Tumor" = "red",
        "pCR" = "green",
        "In Process" = "violet"
      ))
  );
}

`[.rrbsDataList` <- function(x, i){
  cat('[.rrbsDataList\n');
  #cat("i =", i, "\n");
  #cat("class(i) =", class(i), "\n");

  indicies <- if (is.logical(i))
  {
    which(i);
  } else if (is.numeric(i)) {
    i;
  } else if (is.character(i)) {
    match(i, names(x));
  }
  cat("indicies = ", indicies, "\n");

  new(
    "rrbsDataList",
    x@.Data[indicies],
    samples = x@samples[indicies],
    description = x@description);
};

rrbsDataList <- setClass(
  "rrbsDataList",
  contains = "list",
  slots = c(
    samples = "rrbsSampleList",
    description = "character"
    ));

setMethod(
  "show",
  signature=c("rrbsDataList"),
  function(object){
    show(new("methylRawList", object));
  });

setMethod(
  "names",
  signature = c("rrbsDataList"),
  function(x){
    res <- sapply(x, attr, "sample.id");
    names(res) <- NULL;
    res;
  });

setGeneric(
  "description",
  function(x) standardGeneric("description"))

setMethod(
  "description",
  signature = c("rrbsDataList"),
  function(x){
    x@description;
  });

setMethod(
  "as.matrix",
  signature("rrbsDataList"),
  function(x){

    m <- methylRawList_to_matrix(x);

    m_coords <- str_match(rownames(m$b.value), "([^:]+):(\\d+)-(\\d+)/(.)")

    m_coords_gr <- GRanges(
      seqnames = m_coords[,2],
      ranges = IRanges(as.numeric(m_coords[,3]), as.numeric(m_coords[,4])),
      strand = "*"
    );
    names(m_coords_gr) <- rownames(m);

    tss_gr_filename <- "RData/2015-12/tss_gr.RData";
    load(file=tss_gr_filename, verbose=T); # загружаем структуру данных в формате RData

    h <- findOverlaps(m_coords_gr, tss_gr)

    # c новой версией случилась замена названий слотов subjectHits и queryHits
    # объекта SortedByQueryHits с потерей обратной совместимости.
    h.to <- if ("subjectHits" %in% slotNames(h)) { h@subjectHits } else { h@to };
    h.from <- if ("queryHits" %in% slotNames(h)) { h@queryHits } else { h@from };

    genes <- llply(
      1:nrow(m$b.value),
      function(i){
        genes_str <- paste(names(table(tss_gr[h.to[h.from == i]]$geneName)), collapse=", ");
        genes_str;
      });

    res <- new(
      "rrbsMatrix",
      m$b.value,
      coverage = m$coverage,
      samples = x@samples,
      coords = m_coords_gr,
      genes = genes,
      depth_threshold = x@depth_threshold,
      width_threshold = x@width_threshold,
      description = description(x));

    res;
  });

#' Loads .cov files of bismark format to S4 list like object.
#'
#' @param ss A data.frame with sample properties (metadata).
#'           Column named 'Cov_FileName' is mandatory.
#' @param cpg A
loadRrbsData <- function(ss, cpg, description, ...){

  files <- sapply(ss, "[[", "Cov_FileName")
  if (sum(!file.exists(files)) > 0)
    stop(as.character(sprintf("Some files are missed:\n %s",paste(with(list(fns=files), fns[!file.exists(fns)]), collapse="\n"))));

  rrbs_data_list <- llply(
    1:length(ss),
    function(i, ss, cpg)
    {
      si <- ss[[i]];
      switch(
        as.character(si$Src),
        epi=load_epi2(
          filename=as.character(si$Cov_FileName),
          sample.id=as.character(si$Sample_Name_Short),
          cpg,
          assembly="hg19", context="CpG", resolution="base",
          ...
        ),
        ucsc="load_ucsc"
      );
    },
    ss, cpg, .parallel = TRUE);
  names(rrbs_data_list) <- names(ss);

  new(
    "rrbsDataList",
    rrbs_data_list,
    samples = ss,
    description = description
  );
};

rrbsDataFiltered <- setClass(
  "rrbsDataFiltered",
  contains = "rrbsDataList",
  slots = c(
    width_threshold = "numeric",
    depth_threshold = "numeric"
    ));

setMethod(
  "description",
  signature = c("rrbsDataFiltered"),
  function(x){
    #cat(class(x), "\n");
    #cat(x@width_threshold, "\n");
    sprintf(
      "%s, width: %d, depth: %d",
      x@description,
      x@width_threshold,
      x@depth_threshold);
  });

filterRrbsData <- function(rrbs_data, width_threshold, depth_threshold){

  if (class(rrbs_data) != "rrbsDataList")
    stop("Argument 'rrbs_data' class must be 'rrbsDataList'");

  mrl <- new(
    "methylRawList",
    lapply(
      rrbs_data,
      function(sample_data){
        sample_data[!is.na(sample_data$coverage)];
      }));

  # фильтруем не только по минимальному покрытие, но и по чрезмерному
  mrl_fltDpt <- filterByCoverage(
    mrl,
    lo.count=depth_threshold, lo.perc=NULL,
    hi.count=NULL, hi.perc=99.95
  );

  # фильтр по ширине (числу покрытых CpG в образце),
  # необходимо применять ДО нормализации глубины покрытия
  mrl_fltDpt_fltWidth <- with(
    list(width_flt=sapply(mrl_fltDpt, nrow) > width_threshold),
    {
      # так как мы удаляем образцы, то объект MethylRawList нужно пересоздать
      mrl_fltDpt_fltWdt <- new("methylRawList", mrl_fltDpt[width_flt]);

      # нормализация ПОСЛЕ фильтра по ширине
      mrl_fltDpt_fltWdt_nrm <- normalizeCoverage(mrl_fltDpt_fltWdt);

      new(
        "rrbsDataFiltered",
        new(
          "rrbsDataList",
          mrl_fltDpt_fltWdt_nrm,
          samples=rrbs_data@samples[width_flt],
          description = rrbs_data@description
        ),
        width_threshold=width_threshold,
        depth_threshold=depth_threshold);
    });

  mrl_fltDpt_fltWidth;
};

# уровень метилирования CpG в образцах
setClass(
  "rrbsMatrix",
  contains = "matrix",
  slots = c(
    "coverage" = "matrix",
    "samples" = "rrbsSampleList",
    "coords" = "GRanges",
    "genes" = "list",
    "depth_threshold" = "numeric",
    "width_threshold" = "numeric",
    "description" = "character"
  ),
  validity = function(object){
    if (is.null(object@.Data)) {
      ".Data CAN NOT be null"
    }
    else if (is.null(object@coverage)) {
      "coverage CAN NOT be null"
    }
    else if (is.null(object@samples)) {
      "samples CAN NOT be null"
    }
    else if (is.null(object@coords)){
      "coords CAN NOT be null"
    }
    else if (is.null(object@genes)) {
      "genes CAN NOT be NULL"
    }
    else if (sum(c(dim(object@.Data) != dim(object@coverage))) != 0) {
      "'.Data' and 'coverage' dimensions and size MUST be equal."
    }
    else if (nrow(object@.Data) != length(object@coords)){
      "length(coords) MUST be equal to nrow(.Data)"
    }
    else if (nrow(object@.Data) != length(object@genes)) {
      "length(genes) MUST be equal to nrow(.Data)"
    }
    else if (ncol(object@.Data) != length(object@samples)) {
      "length(samples) MUST be equal to ncol(.Data)"
    }
    else TRUE;
  });

setMethod(
  "description",
  signature = c("rrbsMatrix"),
  function(x){
    cat("Override: ", class(x), "\n");
    x@description;
  });

`[.rrbsMatrix` <- function(x, i, j){
  # cat('[.rrbsMatrix\n');
  # if (missing(i)) {
  #   cat("i is missing\n")
  # }
  # else {
  #   cat("i =", i, "\n")
  # }
  #
  #   if (missing(j)) {
  #     cat("j is missing\n")
  #   } else {
  #     cat("j =", j, "\n")
  #   }

  if (!missing(i) && !missing(j)){
    new("rrbsMatrix",
        x@.Data[i,j, drop=FALSE],
        coverage = x@coverage[i,j, drop=FALSE],
        samples = x@samples[j],
        coords = x@coords[i],
        genes = x@genes[i],
        depth_threshold = x@depth_threshold,
        width_threshold = x@width_threshold,
        description=x@description)
  } else if( missing(i) && !missing(j)) {
    new("rrbsMatrix",
        x@.Data[,j, drop=FALSE],
        coverage = x@coverage[,j, drop=FALSE],
        samples = x@samples[j],
        coords = x@coords,
        genes = x@genes,
        depth_threshold = x@depth_threshold,
        width_threshold = x@width_threshold,
        description=x@description)
  } else if (!missing(i) && missing(j)) {
    new("rrbsMatrix",
        x@.Data[i,, drop=FALSE],
        coverage = x@coverage[i,, drop=FALSE],
        samples = x@samples,
        coords = x@coords[i],
        genes = x@genes[i],
        depth_threshold = x@depth_threshold,
        width_threshold = x@width_threshold,
        description=x@description)
  } else
    stop("Not implemented");
};

#' Удаляем образцы, по которым нельзя посчитать расстояния
#' с другими образцами. Начинаем с самых несвязанных
#' @param  x rrbsMatrix object to filter
#' @param min threshold for min shared CpG count
filter_disjoint <- function(x, min_shared_cpg = 1){
  if (!is(x, "rrbsMatrix"))
    stop("Argument 'x' class 'rrbsMatrix' only allowed");

  count_shared_cpg <- function(x){
    with(
      list(xx = ifelse(is.na(x@.Data), 0, 1), n = ncol(x)),
      ldply(
        1:(n-1),
        function(i)
        {
          # построим список каждый-с-каждым
          ddply(
            data.frame("i"=rep(i, n-i), "j"=(i+1):n),
            c("i","j"),
            function(row){
              i <- row[1,1]; j <- row[1,2];
              sum(xx[,i] * xx[,j]);
            });
        },
        .parallel=TRUE));
  }

  # count shared CpG pairwise
  shared_cpg <- count_shared_cpg(x);

  while (sum(shared_cpg[,3] < min_shared_cpg) > 0){
    #stop("Not implemented");
    out_i <- which.min(sapply(1:ncol(x), function(i){ sum(shared_cpg[shared_cpg[,1] == i | shared_cpg[,2] == i, 3]); }));
    cat("out_i =", out_i, "\n");
    # x <- x[, -out_i];
    x <- x[, !1:ncol(x) %in% c(out_i)];
    shared_cpg <- count_shared_cpg(x);
  }

  new(
    "rrbsMatrix",
    x@.Data,
    coverage = x@coverage,
    samples = x@samples,
    coords = x@coords,
    genes = x@genes,
    depth_threshold = x@depth_threshold,
    width_threshold = x@width_threshold,
    description = sprintf("%s, min shared CpGs: %d", description(x), min_shared_cpg)
    );
}


filter_miss <- function(x, max_miss = 0.2){
  if (!is(x, "rrbsMatrix"))
    stop("Argument 'x' class 'rrbsMatrix' only allowed");

  miss_x <- miss(x);

  while (
    #(sum(cvrg_smpl_flt) + sum(sd_smpl_flt) + sum(cvrg_c_flt) + sum(sd_c_flt)) > 0 | miss > 0.25
    miss_x > max_miss
  ) {

    xx <- ifelse(is.na(x@.Data), 0, 1);

    # в скольких образцах есть данные по каждому локусу
    cvrg_c <- apply(xx, 1, sum);
    cvrg_c_mean <- mean(cvrg_c);
    cvrg_c_sd <- sd(cvrg_c);
    filter_out_plan <- nrow(x)*(miss_x - max_miss);

    cvrg_c_thrsh <- cvrg_c[tail(head(order(cvrg_c), filter_out_plan) ,1)];
    cvrg_c_flt <- cvrg_c <= cvrg_c_thrsh;

    x <- x[!cvrg_c_flt, ];

    miss_x <- miss(x);

    cat(sprintf(
      "Filters out: #plan=%f, #fact=%d CpGs  ==>  Dim: %d x %d (samples x CpGs); miss(x) = %.4f \n",
      filter_out_plan, sum(cvrg_c_flt), dim(x)[2], dim(x)[1], miss_x ));
  }

  # TODO: теперь удалим CpG-пары с NA-расстояниями до других CpG-пар
  # критика: если CpG-пар очень много, то матрица расстояний будет слишком большой


  new(
    "rrbsMatrix",
    x@.Data,
    coverage = x@coverage,
    samples = x@samples,
    coords = x@coords,
    genes = x@genes,
    depth_threshold = x@depth_threshold,
    width_threshold = x@width_threshold,
    description=sprintf("%s, max missed: %f", description(x), max_miss));
}


setGeneric(
  "miss",
  function(x) standardGeneric("miss"));

setMethod(
  "miss",
  signature = c("rrbsMatrix"),
  function(x){
    sum(is.na(x))/(sum(!is.na(x))+sum(is.na(x)));
  });

setGeneric(
  "dist_pairwise",
  function(x, method="manhattan") standardGeneric("dist_pairwise"));

setMethod(
  "dist_pairwise",
  signature=c("rrbsDataList"),
  function(x, method="manhattan"){

    ix <- 1:length(x);
    n <- length(ix);

    dists <- ldply(
      1:(n-1),
      function(i)
      {
        # построим список каждый-с-каждым (но только один раз)
        ddply(
          data.frame(
            "i"=rep(i, n-i),
            "j"=(i+1):n
          ),
          c("i","j"),
          function(row){
            s1 <- with(list(v=x[[row[1,1]]]), v[!is.na(v$coverage)]);
            s2 <- with(list(v=x[[row[1,2]]]), v[!is.na(v$coverage)]);

            pd <- new("methylRawList",list(s1,s2));
            pm <- unite(pd, destrand=FALSE);

            if (method == "pearson") {
              1 - (1 + cor(pm$numCs1/pm$coverage1, pm$numCs2/pm$coverage2, method="pearson"))/2;
            } else if (method == "manhattan") {
              sum(abs(pm$numCs1/pm$coverage1 - pm$numCs2/pm$coverage2))/nrow(pm);
            } else {
              stop("Distance method error");
            };
          });
      },
      .parallel=TRUE);

    l <- (n/2)*(n-1);
    if (nrow(dists) != l)
      stop("Invalid operation");

    names(dists) <- c("i", "j", "dist");

    structure(
      as.numeric(dists$dist),
      Size = n, Labels = names(x),
      Diag = FALSE, Upper = FALSE,
      method = "user", class = "dist");
  });

setMethod(
  "dist_pairwise",
  signature=c("rrbsMatrix"),
  function(x, method="manhattan"){
    m <- x;

    ix <- 1:ncol(m);
    n <- length(ix);

    dists <- ldply(
      1:(n-1),
      function(i)
      {
        # построим список каждый-с-каждым (но только один раз)
        ddply(
          data.frame(
            "i"=rep(i, n-i),
            "j"=(i+1):n
            ),
          c("i","j"),
          function(row){
            v1 <- m[,row[1,1]];
            v2 <- m[,row[1,2]];

            if (method == "pearson") {
              1- (1 + cor(v1, v2, use="pairwise.complete.obs", method="pearson"))/2;
              }
            else if (method == "manhattan") {
              common_i <- !is.na(v1) & !is.na(v2);
              sum(abs(v1[common_i] - v2[common_i]))/sum(common_i);
              }
            else {
              stop("Distance method error");
              };
          });
      },
      .parallel=TRUE);

    l <- (n/2)*(n-1);
    if (nrow(dists) != l)
      stop("Invalid operation");

    names(dists) <- c("i", "j", "dist");

    structure(
      as.numeric(dists$dist),
      Size = n, Labels = names(x@samples),
      Diag = FALSE, Upper = FALSE,
      method = "user", class = "dist");
  });

setGeneric(
  "plot_tree",
  function(x, types, colors, distance.method="manhattan", method="ward.D2", ...)
    standardGeneric("plot_tree"));

setMethod(
  "plot_tree",
  signature = c("rrbsDataList"),
  function(x, types, colors, distance.method="manhattan", method="ward.D2", ...)
  {
    if (length(x) != length(types))
      stop("Arguments 'x' and 'types' must be the same length.");

    cat("distance.method =", distance.method, ", method =", method, "\n");

    distances <- dist_pairwise(x, method = distance.method);
    hc <- hclust(distances, method=method);
    my.plot.hclust(
      hc, types, colors,
      main=sprintf(
        "CpG methylation clustering\n%s",
        description(x)),
      ...
    );
  });

setMethod(
  "plot_tree",
  signature = c("rrbsMatrix"),
  function(x, types, colors, distance.method="manhattan", method="ward.D2", ...)
  {
    if (ncol(x) != length(types))
      stop("Arguments 'x' and 'types' must be the same length.");

    distances <- dist_pairwise(x, method = distance.method);
    hc <- hclust(distances, method=method);
    my.plot.hclust(
      hc, types, colors,
      main=sprintf(
        "CpG methylation clustering\n%s",
        description(x)), ...);
  });

#' Get features adjusted p-value (FDR adjusted) to distinct groups
#' @param x Data
#' @param sample_group Vector with group name for each sample
setGeneric(
  "get_features",
  function(x, sample_group, ...)
    standardGeneric("get_features"));

setMethod(
  "get_features",
  signature = c("rrbsMatrix", "character"),
  function(x, sample_group, meth_trshld = 0.2)
  {
    if (length(sample_group) != ncol(x))
      stop("Argument 'sample_group' length must agree data sample count");

    group_names <- with(list(v=names(table(sample_group))), { v[v != "" & !is.na(v)]; });

    groups <- {
      r <- lapply(group_names, function(group_name) { which(sample_group == group_name);});
      names(r) <- group_names;
      r; }

    m <- melt(x, varnames=c("coord", "Sample"), value.name="M");

    g_meth <- aaply(
      x@.Data, 1,
      function(x_row){
        sapply(groups, function(group_samples) { sum(x_row[group_samples] > meth_trshld, na.rm=TRUE); });
      });

    g_unmeth <- aaply(
      x@.Data, 1,
      function(x_row){
        sapply(groups, function(group_samples) { sum(x_row[group_samples] <= meth_trshld, na.rm=TRUE); });
      });

    g_val <- aaply(
      x@.Data, 1,
      function(x_row){
        sapply(groups, function(group_samples) { mean(x_row[group_samples], na.rm=TRUE); });
      });

    system.time(
    res <- ldply(
      1:nrow(x),
      function(i_row, meth, unmeth, val, rows){
        tbl <- rbind(
          meth[i_row, ],
          unmeth[i_row, ]
        );

        p.value=fisher.test(tbl, workspace = 2e9)$p.value;
        c(val[i_row,], p.value=p.value);
      },
      meth=g_meth[match(rownames(g_meth), rownames(x)),],
      unmeth=g_unmeth[match(rownames(g_unmeth), rownames(x)),],
      val=g_val[match(rownames(g_val), rownames(x)),],
      rows=rownames(x),
      .parallel=TRUE));
    rownames(res) <- rownames(x);
    res$adj.p.value <- p.adjust(res$p.value);
    res;
  });

# --- cpgs ---
split_coords <- function(v){
  vals <- strsplit(gsub(":|-|/", ",", v), ",")

  GRanges(
    seqnames = sapply(vals, "[", 1),
    ranges = IRanges(
      start = as.integer(sapply(vals, "[", 2)),
      end = as.integer(sapply(vals, "[", 3))),
    strand = sapply(vals, "[", 4)
  );
};

load_coverage_epi <- function(filename, ranges){
  # функциЯ загружает данные для ежика покрытия иЗ файла .bam.bed
  reads <- import.bed(filename, asRangedData=FALSE);

  cov <- ranges;
  cov$coverage <- countOverlaps(cov, reads);
  cov <- as.data.frame(cov[cov$coverage > 0]);
  cov$width <- cov$end - cov$start;
  cov;
}

plot.lib <- function(lib_data_list, max_len=300, bin_width=10,
                     xlab="XmaI fragment length (between sites)", ...)
{
  for (i in 1:length(lib_data_list)) {
    #i <- (page_i-1)*page_size + page_smpl_i;

    if (i <= length(lib_data_list)){
      with(
        list(
          d = lib_data_list[[i]],
          sample_name = names(lib_data_list)[[i]]
        ), {
          d <- d[d$coverage > 0 & d$width < max_len,]

          g <- llply(
            1:(max_len %/% bin_width),
            function(i) {
              d[(i-1)*bin_width <= d$width & d$width < i*bin_width, ]$coverage;
            });

          y_max = quantile(d$coverage, probs=c(0.99))
          #hist(d$coverage, breaks=1000)

          plot(
            c(), xlim=c(0, max_len/bin_width), ylim=c(0,y_max), xaxt="n", yaxt="n",
            xlab=xlab,
            ylab="Coverage (read count)",
            main=sample_name,
            ...);

          boxplot(
            g, add=T, outline=FALSE, las = 2,
            names=sapply(
              1:(max_len %/% bin_width),
              function(i){ if ((i %% 1)==0) { (i)*bin_width;} else { NA; };}
            ),
            ...);
        });
    }
  }
}
