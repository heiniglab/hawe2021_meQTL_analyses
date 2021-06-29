#' Gets trans-associated CpGs for a sentinel SNP
#'
#' Get all trans cpgs for an individual sentinel SNP using a "cosmo"
#' object.
#'
#' @param sentinel id of the sentinel to analyze
#' @param trans.meQTL the pruned table of trans-associations
#' @param cosmo the cosmo object containing the individual associations
#' @param cpgs the probe ranges of the cpgs to check
#' @param cosmo.idxs flag whether to return the indices in the cosmo object rather than the ones
#' for the cpg list probe ranges.
#'
get.trans.cpgs <- function(sentinel, trans.meQTL, cosmo, cpgs=NULL, cosmo.idxs=F) {

  pairs = which(trans.meQTL[,"sentinel.snp"] == sentinel)

  trans.snp.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.snp"], sep=""),
                             ranges=IRanges(trans.meQTL[,"interval.start.snp"],
                                            trans.meQTL[,"interval.end.snp"]))
  trans.cpg.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.cpg"], sep=""),
                             ranges=IRanges(trans.meQTL[,"interval.start.cpg"],
                                            trans.meQTL[,"interval.end.cpg"]))

  pair.snps = GRanges(seqnames=paste("chr", cosmo[,"snp.chr"], sep=""),
                      ranges=IRanges(cosmo[,"snp.pos"], width=1))
  pair.cpgs = GRanges(seqnames=paste("chr", cosmo[,"cpg.chr"], sep=""),
                      ranges=IRanges(cosmo[,"cpg.pos"], width=2))

  pairs = pair.snps %over% trans.snp.ranges[pairs] &
    pair.cpgs %over% trans.cpg.ranges[pairs]

  if(is.null(cpgs) | cosmo.idxs) {
    return(pairs)
  } else {
    is.meQTL = cpgs %over% pair.cpgs[pairs]
    return(is.meQTL)
  }
}


#' Gets trans-associated SNPs for a sentinel CpG
#'
#' Get all trans SNPs for an individual sentinel CpG using a "cosmo"
#' object.
#'
#' @param sentinel id of the sentinel CpG to analyze
#' @param trans.meQTL the pruned table of trans-associations
#' @param cosmo the cosmo object containing the individual associations
#'
get.trans.snps <- function(sentinel, trans.meQTL, cosmo) {

  pairs = which(trans.meQTL[,"sentinel.cpg"] == sentinel)

  trans.snp.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.snp"], sep=""),
                             ranges=IRanges(trans.meQTL[,"interval.start.snp"],
                                            trans.meQTL[,"interval.end.snp"]))
  trans.cpg.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.cpg"], sep=""),
                             ranges=IRanges(trans.meQTL[,"interval.start.cpg"],
                                            trans.meQTL[,"interval.end.cpg"]))

  pair.snps = GRanges(seqnames=paste("chr", cosmo[,"snp.chr"], sep=""),
                      ranges=IRanges(cosmo[,"snp.pos"], width=1))
  pair.cpgs = GRanges(seqnames=paste("chr", cosmo[,"cpg.chr"], sep=""),
                      ranges=IRanges(cosmo[,"cpg.pos"], width=2))

  pairs = pair.snps %over% trans.snp.ranges[pairs] &
    pair.cpgs %over% trans.cpg.ranges[pairs]

  return(pairs)
}

graph.union <- function(g1, g2) {
  add.nodes.to.g1 = setdiff(nodes(g2), nodes(g1))
  g1 = addNode(add.nodes.to.g1, g1)
  add.nodes.to.g2 = setdiff(nodes(g1), nodes(g2))
  g2 = addNode(add.nodes.to.g2, g2)

  if (F) {
    ## this should be replaced with the edgeMatrix function
    e = edgeMatrix(g2)
    from = nodes(g2)[e[1,]]
    to = nodes(g2)[e[2,]]
    g = addEdge(from, to, g)
  }

  return(union(g1, g2))
}


#' Get CpGs that match the mean and standard deviation of the methylation level
#'
#' Given a CpG name or index and a table of mean and standard deviations return
#' all CpGs that are similar
#'
#' @param cpg index or name of a CpG
#' @param msd data.frame with columns 'meanBeta' and 'sdBeta' for each CpG
#' @param eps.m tolerance (absolute) on the mean value
#' @param eps.sd tolerance (absolute) on the sd value
#'
#' @return a set of indices or rownames of matching CpGs
#' @export
get.matched.cpgs <- function(cpg, msd, eps.m=0.05, eps.sd=0.05) {
  selected = which(abs(msd[,"meanBeta"] - msd[cpg, "meanBeta"]) < eps.m &
    abs(msd[,"sdBeta"] - msd[cpg, "sdBeta"]) < eps.sd)
  return(setdiff(selected, match(cpg, rownames(msd))))
}

#' Get SNPs that match the mean and standard deviation of the allele frequency
#'
#' Given a SNP name or index and a table of mean and standard deviations return
#' all SNPs that are similar
#'
#' @param snp index or rsids
#' @param allele.freq data.frame with columns 'rsid' and 'AF' for each SNP
#' @param eps tolerance (absolute) on the allele frequency
#' @param af_column Optional. Custom column name for allele frequency
#'
#' @return a set of indices or rownames of matching SNPs
#' @export
get.matched.snps <- function(snp, allele.freq, eps=0.05, af_column = "AF") {
  handle.characters = is.character(snp)
  if (handle.characters) {
    snp.id = snp
    snp = match(snp.id, allele.freq[,"rsid"])
  }
  selected = which(abs(allele.freq[,af_column] - allele.freq[snp, af_column]) < eps)
  idx = setdiff(selected, snp)
  if (handle.characters) {
    idx = allele.freq[idx,"rsid"]
  }
  return(idx)
}

#' More generic function to get matched entries
#'
#' Given a rowname or index and a table of values return
#' all rows that are similar
#'
#' @param query index or name of the query row
#' @param data data.frame containing the info to match for
#' @param cols the columns to match (if NULL all columns)
#' @param eps tolerance (absolute) on each colum (length 1 or ncol(data))
#'
#' @return a set of indices or rownames of matching rows
#' @export
get.matched <- function(query, data, cols=NULL, eps=0.05) {
  if (is.character(query)) {
    query = match(query, rownames(data))
  }
  if (is.null(cols)) {
    cols = 1:ncol(data)
  }
  stopifnot(length(eps) == length(cols) || length(eps) == 1)
  if (length(eps) == 1) {
    eps = rep(eps, length(cols))
  }
  selected = 1:nrow(data)
  for (i in 1:length(cols)) {
    col = cols[i]
    selected = intersect(which(abs(data[,col] - data[query, col]) < eps[i]), selected)
  }
  return(setdiff(selected, query))
}



#' Resampling based P-values for TFBS enrichment
#'
#' Use resampling on matched background sets to compute empirical P-values
#' for tfbs enrichment. P-values are obtained by looking at the absolute value
#' of the difference of the two proportions.
#'
#' @param is.meQTL boolean vector of the length of CpGs indicating an meQTL
#' @param is.bound boolean vector of the length of CpGs indicating a TFBS
#' @param matched list for all CpGs each containing indices of matched CpGs
#' @param alternative can be "greater" or "less"
#' @param n.resample number of resamplings
#' @param batches number of resamplings before checking
#' @param over percentage of resampled enrichments more extreme than observed
#'             to stop the resampling ahead of time
#'
#' @return empirical P-value
#' @export
empirical.enrichment <- function(is.meQTL, is.bound, matched, alternative="greater", n.resample=10000, batches=100, over=0.2) {
  require(parallel)
  obs = table(is.meQTL, is.bound)

  if (alternative == "greater") {
    cmp <- function(x, y) {
      as.numeric(x >= y)
    }
  } else {
    cmp <- function(x, y) {
      as.numeric(x <= y)
    }
  }

  bg = 0
  for (i in 1:n.resample) {
    set.seed(i) ## for the parallel environment to be reproducible
    bg.set = sapply(matched, sample, 1)
    tab = table(is.bound[bg.set])
    bg = bg + cmp(tab[2], obs[2,2])

    ## check for early breaking
    if (i %% batches == 0) {
      if (bg > over * i) {
        break
      }
    }
  }
  p = (bg + 1) / (i + 1)
  return(p)
}


## this function is from destiny and uses the arpack function
eig.decomp <- function(M, n.eigs, sym) {
  n = nrow(M)
  f <- function(x, extra=NULL) as.matrix(M %*% x)
  wh <- if (sym) 'LA' else 'LM'
  ## constraints: n >= ncv > nev
  ar <- arpack(f, sym = sym, options = list(
               which = wh, n = n, ncv = min(n, 4*n.eigs), nev = n.eigs + 1))
  if (!sym) {
    ar$vectors <- Re(ar$vectors)
    ar$values  <- Re(ar$values)
  }

  ## check for negative values and flip signs
  neg = which(ar$values < 0)
  for (n in neg) {
    ar$values[n] = -ar$values[n]
    ar$vectors[,n] = -ar$vectors[,n]
  }

  return(ar)
}

graph2sparseMatrix <- function(g) {
  em = edgeMatrix(g)
  Asparse = sparseMatrix(em[1,], em[2,], x=1, dims=rep(numNodes(g), 2), dimnames=list(nodes(g), nodes(g)))

  ## make symmetric
  Asparse = Asparse + t(Asparse)
  Asparse[Asparse > 1] = 1

  return(Asparse)
}


## approximate the infinite sum
## in the paper and notes this is matrix is called M
## n.eigs gives the number of eigenvectors used for the approximation
## if this is smaller than 2, the pseudo inverse will be computed
## from and to are the nodes for which we need the transition probabilties
## this avoids to compute the full propagation matrix which is to large
## when from and to are NULL the full matrix will be computed
## sum can be "none", "from", "to" or "both" and will sum up the propagation
## matrix over one or the other or both lists
## "none" returns a "from" by "to" matrix
## "from" returns a vector of length nrow(Asparse)
## "to" returns a vector of length nrow(Asparse)
## "both" returns a nrow(Asparse) by 2 matrix with the from and to vectors
propagation <- function(Asparse, n.eigs=20, from=NULL, to=NULL, sum="none") {
  require(igraph)
  require(Matrix)
  ## transition matrix
  ## transition = Asparse / rowSums(Asparse)

  ## symmetric transition matrix
  D.root = Diagonal(nrow(Asparse), 1 / sqrt(Matrix::rowSums(Asparse)))
  Psym = D.root %*% Asparse %*% D.root

  if (is.null(from)) {
    from = 1:nrow(Asparse)
  }
  if (is.null(to)) {
    to = 1:nrow(Asparse)
  }
  if (is.logical(from)) {
    from = which(from)
  }
  if (is.logical(to)) {
    to = which(to)
  }
  dnames = list(from, to)
  ## if we have character we need to match
  if (is.character(from)) {
    from = match(from, rownames(Asparse))
  }
  if (is.character(to)) {
    to = match(to, colnames(Asparse))
  }

  ## setup different ways of summarizing the propagation matrix
  if (sum == "from") {
    transform = rowSums
    propagation = double(0, length(to))
    names(propagation) = dnames[[2]]
  } else if (sum == "to") {
    transform = colSums
    propagation = double(0, length(from))
    names(propagation) = dnames[[1]]
  } else if (sum == "both") {
    propagation = matrix(0, nrow=nrow(Asparse), ncol=2, dimnames=list(rownames(Asparse), c("from", "to")))
  } else if (sum == "none") {
    transform = function(x) {return(x)}
    propagation = matrix(0, nrow=length(from), ncol=length(to), dimnames=dnames)
  }

  if (n.eigs > 0) {
    ## approximate using the eigen vectors

    if (n.eigs < nrow(Psym)) {
      ## just use a subset of eigen vectors
      eig.M <- eig.decomp(Psym, n.eigs, TRUE)

    } else {
      ## use all eigen vectors (only for small matrices)
      eig.M <- eigen(Psym)
    }

    ## for both computations we discard the first eigenvector as it
    ## represents the stationary distribution
    if (sum %in% c("none", "from", "to")) {
      for (i in 2:n.eigs) {
        increment = (eig.M$values[i] / (1 - eig.M$values[i])) *
          eig.M$vectors[from,i] %*% t(eig.M$vectors[to,i])
        propagation = propagation + transform(increment)
      }

    } else {
      ## when sum is "both" we iterate over the "from" and "to" nodes and add
      ## the from and to vectors incrementally to save memory
      for (i in 2:n.eigs) {
        weight = (eig.M$values[i] / (1 - eig.M$values[i]))
        ## increment the "from" sum vector
        for (ff in from) {
          propagation[,1] = propagation[,1] + weight * (eig.M$vectors[ff,i] * eig.M$vectors[,i])
        }
        ## increment the "to" sum vector
        for (tt in to) {
          propagation[,2] = propagation[,2] + weight * (eig.M$vectors[tt,i] * eig.M$vectors[,i])
        }
      }
    }
  } else {
    ## compute using the pseudo inverse

    ## the first eigenvector corresponds to stationary distribution defined
    ## by the degrees of the nodes
    ## Attention: this is actually the first eigenvector of the asymmetric
    ## Psym matrix
    d = rowSums(Asparse)
    v = sum(d)
    phi0 = d / v

    ## there is an equivalence of the eigenvectors of the symmetric and
    ## asymetric matrix:
    ## EV(sym) = D^-0.5 EV(asym)
    phi0 =  phi0 / sqrt(d)

    ## in the dtp package there is a length normalization step that matches
    ## the vectors exactly (normalizing to unit length vectors)
    phi0 = phi0 / sqrt(sum(phi0^2))

    n <- nrow(Psym)
    inv <- solve(Diagonal(n) - Psym + phi0 %*% t(phi0))
    propagation = inv - Diagonal(n)
    propagation = propagation[from, to]
    if (sum %in% c("none", "from", "to")) {
      propagation = transform(propagation)
    } else {
      stop("sum = 'both' not implemented yet for pseudo inverse")
    }
  }
  return(propagation)
}





## alternatively we can compute expected hitting times for the random walks
## n.eigs gives the number of eigenvectors used for the approximation
## if this is smaller than 2, the pseudo inverse will be computed
## from and to are the nodes for which we need the transition probabilties
## this avoids to compute the full propagation matrix which is to large
## when from and to are NULL the full matrix will be computed
## sum can be "none", "from", "to" or "both" and will sum up the propagation
## matrix over one or the other or both lists
## "none" returns a "from" by "to" matrix
## "from" returns a vector of length nrow(Asparse)
## "to" returns a vector of length nrow(Asparse)
## "both" returns a nrow(Asparse) by 2 matrix with the from and to vectors
## the i-th component of the "from sum" represents the average hitting time
## going from node i to any of the "to" nodes
## the i-th component of the "to sum" represents the average hitting time
## going from any "from" nodes to node i
hitting.time <- function(Asparse, n.eigs=20, from=NULL, to=NULL, sum="none") {

  ## we need the number of edges and the number of nodes and the degree of nodes
  numEdges = sum(Asparse) / 2 + sum(diag(Asparse))
  numNodes = nrow(Asparse)
  deg = rowSums(Asparse)

  ## symmetric transition matrix
  D.root = Diagonal(nrow(Asparse), 1 / sqrt(rowSums(Asparse)))
  Psym = D.root %*% Asparse %*% D.root

  if (is.null(from)) {
    from = 1:numNodes
  }
  if (is.null(to)) {
    to = 1:numNodes
  }
  if (is.logical(from)) {
    from = which(from)
  }
  if (is.logical(to)) {
    to = which(to)
  }
  dnames = list(from, to)
  ## if we have character we need to match
  if (is.character(from)) {
    from = match(from, rownames(Asparse))
  }
  if (is.character(to)) {
    to = match(to, colnames(Asparse))
  }

  ## setup different ways of summarizing the hitting.time matrix
  if (sum == "from") {
    transform = rowMeans
    hitting.time = double(0, length(to))
    names(hitting.time) = dnames[[2]]
  } else if (sum == "to") {
    transform = colMeans
    hitting.time = double(0, length(from))
    names(hitting.time) = dnames[[1]]
  } else if (sum == "both") {
    hitting.time = matrix(0, nrow=numNodes, ncol=2, dimnames=list(rownames(Asparse), c("from", "to")))
  } else if (sum == "none") {
    transform = function(x) {return(x)}
    hitting.time = matrix(0, nrow=length(from), ncol=length(to), dimnames=dnames)
  }

  ## approximate using the eigen vectors

  if (n.eigs < nrow(Psym)) {
    ## just use a subset of eigen vectors
    eig.M <- eig.decomp(Psym, n.eigs, TRUE)

  } else {
    ## use all eigen vectors (only for small matrices)
    eig.M <- eigen(Psym)
  }

  ## Hitting time computed according to Theorem 3.1 from
  ## Lov??sz, L. (1993). Random walks on graphs. Combinatorics.

  compute.hitting.time <- function(from.nodes, to.nodes) {
    sapply(to.nodes, function(tt)
           sapply(from.nodes, function (ff) {
             ht = 2 * numEdges * sum(sapply(2:n.eigs, function(i) {
               v = eig.M$vectors[,i]
               return(1 / (1 - eig.M$values[i]) *
                      (v[tt]^2 / deg[tt] - v[ff] *
                       v[tt] / sqrt(deg[ff] * deg[tt])))
             }))
           }))
  }

  if (sum %in% c("none", "from", "to")) {
    for (i in 2:n.eigs) {
      increment = compute.hitting.time(from, to)
      hitting.time = hitting.time + transform(increment)
    }

  } else {
    ## increment the "from" sum vector
    for (tt in to) {
      hitting.time[,1] = hitting.time[,1] + compute.hitting.time(1:numNodes, tt)
    }
    hitting.time[,1] = hitting.time[,1] / length(to)

    ## increment the "to" sum vector
    for (ff in from) {
      hitting.time[,2] = hitting.time[,2] + compute.hitting.time(ff, 1:numNodes)
    }
    hitting.time[,2] = hitting.time[,2] / length(from)
  }

  return(hitting.time)
}


## find the shortest path with minimal node weight
min.node.weight.path <- function(g, weights, from, to) {
  ## the problem can be transformed into a directed graph problem where all
  ## incoming edges are assigned the node weight

  require(graph)
  require(RBGL)

  h = graphNEL(nodes(g), edgemode="directed")
  em = edgeMatrix(g)
  h = addEdge(nodes(g)[em[1,]], nodes(g)[em[2,]], h, weights[em[2,]])
  h = addEdge(nodes(g)[em[2,]], nodes(g)[em[1,]], h, weights[em[1,]])

  return(sp.between(h, from, to))
}


## plotting networks
plot.propagation.network <- function(g, prop, cpgs, best.trans, trans.genes, sentinel, prefix, additional.genes=NULL, additional.bordercol=NULL) {
  require(graph)
  require(Rgraphviz)
  require(RColorBrewer)

  if (!file.exists(dirname(prefix))) {
    dir.create(dirname(prefix), recursive=T)
  }

  ## now we get the node weights by adding up the from and to weights
  node.weight = rowSums(prop)

  ## finally we would like to find the shortest path with maximal weight
  ## we have an algorithm that finds minimum node weight paths so we need
  ## to turn the weighting around
  ## in addition weights need to be non-negative
  node.weight = max(node.weight) - node.weight + 1

  sp = min.node.weight.path(g, node.weight, from=cpgs, to=best.trans)

  plot.nodes = unique(c(setdiff(unlist(lapply(sp, "[", "path_detail")), NA), sentinel, trans.genes, additional.genes))

  plot.edges = NULL
  for (i in sp) {
    if (length(i$path_detail) > 0) {
      for (j in 1:(length(i$path_detail) - 1)) {
        plot.edges = rbind(plot.edges,
          data.frame(from=i$path_detail[j], to=i$path_detail[j + 1],
                     stringsAsFactors=F))
      }
    }
  }

  ## plot.graph = graphNEL(nodes=unique(c(plot.nodes, trans.genes)))
  ## plot.graph = addEdge(sentinel, trans.genes, plot.graph)
  ## plot.graph = addEdge(plot.edges[,1], plot.edges[,2], plot.graph)
  plot.graph = g
  plot.graph = graph::addNode(sentinel, plot.graph)
  plot.graph = graph::addEdge(sentinel, trans.genes, plot.graph)
  plot.graph = graph::subGraph(plot.nodes, plot.graph)

  ## instead of selecting just one shortest path, we could also use a modified
  ## Dijkstra algorithm that finds all shortest paths

  attrs <- list(node=list(fixedsize=TRUE, fontsize=10, style="filled", fontname="helvetica"), graph=list(overlap="true", root=sentinel, outputorder="edgesfirst"))

  shape = rep("ellipse", numNodes(plot.graph))
  names(shape) = nodes(plot.graph)
  shape[grep("cg", nodes(plot.graph))] = "box"
  shape[sentinel] = "box"

  width = rep(0.8, numNodes(plot.graph))
  names(width) = nodes(plot.graph)
  width[grep("cg", nodes(plot.graph))] = 0.2

  height = rep(0.2, numNodes(plot.graph))
  names(height) = nodes(plot.graph)
  height[grep("cg", nodes(plot.graph))] = 0.2

  label = nodes(plot.graph)
  names(label) = nodes(plot.graph)
  label[grep("cg", nodes(plot.graph))] = ""

  ## fill color represents the random walk score
  pal = brewer.pal(9, "Blues")
  score = prop[match(nodes(plot.graph), rownames(prop)), "from"]
  names(score) = nodes(plot.graph)
  breaks = seq(min(score, na.rm=T), max(score, na.rm=T), length.out=length(pal)+1)
  col = pal[findInterval(score, breaks, rightmost.closed=T)]
  names(col) = nodes(plot.graph)
  col[sentinel] = "#ff0000"
  col[grep("cg", nodes(plot.graph))] = "#00ff00"

  penwidth = rep(1, numNodes(plot.graph))
  names(penwidth) = nodes(plot.graph)
  penwidth[best.trans] = 2

  bordercol = rep("black", numNodes(plot.graph))
  names(bordercol) = nodes(plot.graph)
  bordercol[best.trans] = "red"

  if (!(is.null(additional.genes) || is.null(additional.bordercol))) {
    bordercol[additional.genes] = additional.bordercol
    penwidth[additional.genes] = 2
  }

  nAttrs = list(shape=shape, label=label, width=width, height=height, fillcolor=col, penwidth=penwidth, color=bordercol)

  ecol = rep("black", numEdges(plot.graph))
  names(ecol) = edgeNames(plot.graph)
  ecol[grep(sentinel, names(ecol))] = "red"

  eAttrs = list(color=ecol)


  pdf.file = paste(prefix, ".pdf", sep="")
  dot.file = paste(prefix, ".dot", sep="")
  rdata.file = paste(prefix, ".RData", sep="")

  ## direct plotting plots edges over nodes
  pdf(file=pdf.file)
  plot(plot.graph, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  dev.off()

  ## also use the graphviz command line for plotting
  toDot(plot.graph, dot.file, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  system(paste("twopi -Tpdf -O", dot.file))

  ## plot the color legend
  pdf(file=paste(prefix, "-color-legend.pdf", sep=""))
  plot(c(0, 1), c(0, 9))
  rect(rep(0, 9), 0:8, rep(1, 9), 1:9, col=pal)
  dev.off()

  ## also save the data for later
  save(list=c("nAttrs", "eAttrs", "attrs", "plot.graph"), file=rdata.file)

  ## browser()
}


#' Augment grahs
#'
#' Add TF to CpG edges and trans gene to SNP edges to graphs. This function
#' assumes that you have an object tfbs.ann in the environment. It is not
#' passed as argument (bad style) becuase it is too big (passing by value).
#'
#' @param graphs list of graph objects to add nodes and edges to
#' @param sentinel character id of the sentinel SNP
#' @param trans.genes character vector of genes in the trans locus
#' @param trans.cpgs character vector of CpGs with trans meQTLs
#' @param tfbs.ann logical indicator matrix with CpGs in the rows and
#'        transcription factors as columns
#'
#' @return list of graph objects with the nodes and edges added
#'
#' @imports reshape
#' @imports graph
#' @export
add.to.graphs <- function(graphs, sentinel, trans.genes, trans.cpgs, tfbs.ann) {
  require(reshape)
  require(graph)

  tf.edges = melt(tfbs.ann[trans.cpgs,])
  colnames(tf.edges) = c("cpg", "condition", "adjacent")
  tf.edges = tf.edges[tf.edges[,"adjacent"],]
  for (i in 1:2) {
    tf.edges[,i] = as.character(tf.edges[,i])
  }
  tf = as.character(sapply(strsplit(tf.edges[,"condition"], ".", fixed=T), "[", 1))
  tf.edges = data.frame(tf.edges, tf, stringsAsFactors=F)
  tf.edges = tf.edges[!duplicated(paste(tf.edges[,"tf"], tf.edges[,"cpg"])),]

  ## put together the graphs
  out = list()

  for (graph.idx in 1:length(graphs)) {
    locus.graph = graphs[[graph.idx]]

    ## filter for TFs that are in the graph already
    use.tf.edges = tf.edges[tf.edges[,"tf"] %in% nodes(locus.graph),]

    new.nodes = unique(c(use.tf.edges[,"tf"], trans.cpgs, sentinel, trans.genes))
    new.nodes = setdiff(new.nodes, nodes(locus.graph))

    locus.graph = graph::addNode(new.nodes, locus.graph)

    ## also add some meta data (the type of nodes)
    nodeDataDefaults(locus.graph, "tf") = FALSE
    nodeData(locus.graph, unique(use.tf.edges[,"tf"]), "tf") = TRUE

    nodeDataDefaults(locus.graph, "cpg") = FALSE
    nodeData(locus.graph, trans.cpgs, "cpg") = TRUE

    nodeDataDefaults(locus.graph, "snp") = FALSE
    nodeData(locus.graph, sentinel, "snp") = TRUE

    nodeDataDefaults(locus.graph, "trans.gene") = FALSE
    nodeData(locus.graph, trans.genes, "trans.gene") = TRUE

    ## add edges for the tfbs
    locus.graph = addEdge(use.tf.edges[,"tf"], use.tf.edges[,"cpg"], locus.graph)

    ## add edges for the connection of the locus and its genes
    locus.graph = addEdge(rep(sentinel, length(trans.genes)), trans.genes, locus.graph)

    out[[graph.idx]] = locus.graph
  }

  return(out)
}




#' Match random loci
#'
#' Match random background loci with similar size and number of genes
#'
#' @param width locus widths to be matched
#' @param ngenes number of genes to be matched (same length)
#' @param genes GRanges object with the gene coordinates for the
#'              whole genome (or subset of chromosomes).
#' @param eps tolerance for width matching (relative to width)
#'
#' @return GRanges of matching loci
#'
#' @export
match.loci <- function(width, ngenes, genes, eps, chroms=NULL) {

  if (is.null(chroms)) {
    chroms = seqlevels(genes)
  }

  ## build candidate loci with ngenes
  candidates = GRanges()
  ## do this by chromosome
  for (chr in chroms) {
    chr.genes = genes[which(seqnames(genes) == chr)]
    strand(chr.genes) = "*"
    chr.genes = sort(chr.genes)


    # @TODO: verify this 'try/catch'
    if (ngenes > length(chr.genes)) {
      next()
    }



#>>>>>>> 8954a6f55683565e1ea98f1666e5f2f31a162fc9
    start = start(chr.genes)[-(length(chr.genes) - (ngenes):0)]
    end = end(chr.genes[(ngenes):length(chr.genes)])
    chr.candidates = GRanges(seqnames=chr, ranges=IRanges(start, end))

    candidates = c(candidates, chr.candidates)
  }

  ## make sure that the number of genes really matches
  ng = countOverlaps(candidates, genes)
  candidates = candidates[ng == ngenes]

  ## select by width
  if (eps < Inf) {
    w = width(candidates)
    candidates = candidates[abs(w - width) < eps * width]
  }
  return(candidates)
}


#' Candidate prioritization by random walk analysis
#'
#' @param locus.graph graphNEL object as generated in network_preprocess.R
#' @param largest.cc graphNEL object with the largest connected component
#'        of the protein - protein interaction network
#' @param string.nodes character vector with the names of the nodes in the PPI
#'        network
#' @param tf.nodes character vector with the names of the transcription factors
#' @param n.eigs number of eigen vectors to use for the approximation
#' @param mode character "propagation" or "hitting.time"
#' @param return.g logical should the graph used for analysis be returned?
#'
#' @return list with sentinel, cpgs, tfs, trans.genes, prop, best.trans.gene, best.tf, g
prioritze.by.randomwalk <- function(locus.graph, largest.cc, string.nodes, tf.nodes, n.eigs=500, mode="propagation", return.g=FALSE) {
  require(graph)
  require(igraph)

  n = nodes(locus.graph)

  ## get the sentinel
  sentinel = n[unlist(nodeData(locus.graph, n, "snp"))]

  ## get all the trans locus genes
  trans.genes = n[unlist(nodeData(locus.graph, n, "trans.gene"))]
  trans.genes = intersect(trans.genes, graph::nodes(largest.cc))
  if (length(trans.genes) == 0) {
    warning("no trans genes found in locus graph")
    return(NULL)
  }

  ## also get the cpgs
  cpgs = n[unlist(nodeData(locus.graph, n, "cpg"))]
  ## cpgs = cpgs[grep("cg", cpgs)]
  cpgs = intersect(cpgs, graph::nodes(largest.cc))
  if (length(cpgs) == 0) {
    warning("no cpgs found in locus graph")
    return(NULL)
  }

  ## get the sub graph for the locus (removing snps and all other cpgs)
  ## KAP1 is not connected in the STRING network
  nodes.set = c(string.nodes, setdiff(tf.nodes, "KAP1"), trans.genes, cpgs)
  g = subGraph(intersect(nodes(largest.cc), nodes.set), largest.cc)

  ## check again that we have one connected component
  ig = graph_from_graphnel(g)
  cl = clusters(ig)
  keep = nodes(g)[cl$membership == which.max(cl$csize)]
  if (!all(nodes(g) %in% keep)) {
    missing = setdiff(nodes(g), keep)
    cat("Some nodes are not in the connected component! Please check input!\n")
    cat(missing, sep="\n")
    cat("\n")
    missing = cbind(absolute=sapply(list(trans.genes=trans.genes,
                      string=string.nodes, tfs=tf.nodes), function(x)
                      length(intersect(x, missing))),
      relative=sapply(list(trans.genes=trans.genes, string=string.nodes,
        tfs=tf.nodes), function(x)
        length(intersect(x, missing)) / length(x)))
    print(missing)
    g = subGraph(keep, largest.cc)
    cpgs = intersect(cpgs, nodes(g))
    trans.genes = intersect(trans.genes, nodes(g))
    tf.nodes = intersect(tf.nodes, nodes(g))
  }

  rg = NULL
  if (return.g) {
    rg = g
  }

  ## also get the tfs
  tfs = unique(unlist(adj(g, cpgs)))

  ## now we would like to get the propagation values from the cpgs to the
  ## genes in the trans locus and the tfs
  if (mode == "propagation") {
    prop = propagation(graph2sparseMatrix(g), n.eigs=n.eigs, from=cpgs, to=c(trans.genes, tfs), sum="both")

    best.trans = trans.genes[which.max(prop[trans.genes,"from"])]

    return(list(sentinel=sentinel, cpgs=cpgs, tfs=tfs, trans.genes=trans.genes, prop=prop, best.trans.gene=best.trans, best.tf=tfs[which.max(prop[tfs,"from"])], g=rg))
  }

  if (opt$mode == "hitting.time") {
    ht = hitting.time(graph2sparseMatrix(g), n.eigs=500, from=cpgs, to=c(trans.genes, tfs))
    return(list(sentinel=sentinel, cpgs=cpgs, tfs=tfs, trans.genes=trans.genes, ht=ht, best.trans.gene=trans.genes[which.max(colSums(ht[,trans.genes,drop=F]))], best.tf=tfs[which.max(colSums(ht[,tfs,drop=F]))], g=rg))
  }

}


#' Get nodes by type
#'
#' Convenience function to retrieve nodes of a certain type
#' @param g graphNEL locus graph object
#' @param type character indicating the node type
#'
#' @return character vector with the nodes of given type
get.nodes.by.type <- function(g, type) {
  require(graph)
  n = nodes(g)
  return(n[unlist(nodeData(g, n, type))])
}


#' Prioritization summary table
#'
#' Extract summary information from the random walk analysis
#'
#' @param prioritze list of random walk analysis results generated with the
#'                  network_batch.R script
#' @param full boolean indicating whether all genes in the trans locus should
#'             be reported (TRUE) or only the best (FALSE)
#' @param locus.graphs a list of graphNEL objects with more details on the
#'                     graphs used for analysis (for full summary)
#' @param subset list of subsets of trans genes for each locus this is useful
#'               to recycle permutation runs on maximal intervals
#'
#' @return data.frame with the best candidate gene per locus
#'
#' @export
prioritization.table <- function(prioritize, full=FALSE, locus.graphs=NULL, subset=NULL) {
  require(graph)
  tab = NULL
  for (sentinel in names(prioritize)) {
    if (length(prioritize[[sentinel]]$best.trans.gene) == 0) {
      cat("problem with SNP:", sentinel, "\n")
      next
    }

    if (full) {
      if (is.null(locus.graphs)) {
        stop("need to specify the locus.graphs list when full == TRUE")
      }

      ## get the full list of trans genes (also the ones not in string)
      all.tgenes = get.nodes.by.type(locus.graphs[[sentinel]], "trans.gene")
    } else {
      all.tgenes = with(prioritize[[sentinel]], trans.genes)
    }
    ## check if we are considering only a subset
    if (!is.null(subset)) {
      sid = strsplit(sentinel, ".", fixed=T)[[1]][1]
      ## if (sid == "rs7783715") {
      ##   browser()
      ## }
      all.tgenes = intersect(all.tgenes, subset[[sid]])
    }

    ## put together the details for the tested trans genes
    tested.tgenes = intersect(with(prioritize[[sentinel]], trans.genes),
      all.tgenes)
    if (length(tested.tgenes) == 0) {
      cat("No tested trans genes for SNP:", sentinel, "\n")
      next
    }
    best.tgene = with(prioritize[[sentinel]],
      tested.tgenes[which.max(prop[tested.tgenes,"from"])])

    p = with(prioritize[[sentinel]], prop[tested.tgenes,"from"])
    rescaled = p / sum(p)
    entropy = sum(-rescaled * log(rescaled))
    ## also compute the maximum possible
    u = rep(1/length(p), length(p))
    max.entropy = sum(-u * log(u))

    trans.gene.details = with(prioritize[[sentinel]],
      data.frame(trans.gene=tested.tgenes, prop=prop[tested.tgenes,"from"],
                 p, in.string=TRUE))

    ## add the genes that were not in string
    not.in.string = setdiff(all.tgenes, tested.tgenes)
    if (length(not.in.string) > 0) {
      trans.gene.details = rbind(trans.gene.details,
        data.frame(trans.gene=not.in.string, prop=NA, p=NA, in.string=FALSE))
    }

    smry = with(prioritize[[sentinel]],
      data.frame(sentinel, ncpg=length(cpgs), ntrans=length(tested.tgenes),
                 best.trans=best.tgene, best.tf=best.tf,
                 trans.entropy=entropy, max.entropy=max.entropy,
                 rel.entropy=entropy/max.entropy,
                 best.trans.p=max(rescaled), best.trans.prop=max(p),
                 trans.gene.details))

    ## just get the best trans gene
    if (!full) {
      smry = smry[smry[,"trans.gene"] == best.tgene,]
    }

    tab = rbind(tab, smry)
  }
  return(tab)
}



prioritize.quick <- function(locus.graphs, largest.cc.string, n.eigs=500, return.g=FALSE) {

  ## get all the trans genes
  trans.genes = unique(unlist(lapply(locus.graphs, get.nodes.by.type, "trans.gene")))

  trans.genes.in.string = intersect(trans.genes, nodes(largest.cc.string))

  ## get all the transcription factors connected to any CpG
  tfs = unique(unlist(lapply(locus.graphs, get.nodes.by.type, "tf")))
  tfs.in.string = intersect(tfs, nodes(largest.cc.string))

  ## compute the propagation from transcription factors to all trans genes
  prop = propagation(graph2sparseMatrix(largest.cc.string), n.eigs=n.eigs, from=tfs.in.string, to=trans.genes.in.string, sum="none")

  ## for each locus graph
  rg = NULL
  prio = list()
  for (sentinel in names(locus.graphs)) {

    lg = locus.graphs[[sentinel]]

    ## compute the trans gene scores by weighting each tf by the transsitions
    ## from the connected CpGs
    cpgs = get.nodes.by.type(lg, "cpg")
    cpg.weight = graph::degree(lg, cpgs)

    ## we split the weight on all cpgs with degree > 0
    ## the weight is further split to each TF that is connected to the CpG
    cpg.weight = 1 / (sum(cpg.weight > 0) * cpg.weight)
    cpg.weight[cpg.weight == Inf] = 0
    if (all(cpg.weight == 0)) {
      prio[[sentinel]] = NULL
      next
    }

    ## sum up all the incoming weights for each tf
    tf.weight = sapply(adj(lg, tfs.in.string), function(n) {
      adj.cpgs = intersect(n, cpgs)
      return(sum(cpg.weight[adj.cpgs]))
    })

    ## compute the weighted sums on the trans.genes
    tg = get.nodes.by.type(lg, "trans.gene")
    tg = intersect(tg, nodes(largest.cc.string))
    if (length(tg) == 0) {
      prio[[sentinel]] = NULL
      next
    }
    weighted = colSums(tf.weight * prop[,tg,drop=F])

    best.trans = tg[which.max(weighted)]

    ## put together the results
    if (return.g) {
      rg = lg
    }

    ## make a prop that looks like the one returned by the other method
    prop2col = cbind(to=rep(NA, length(tg)), from=weighted)

    prio[[sentinel]] = list(sentinel=sentinel, cpgs=cpgs, tfs=tfs.in.string, trans.genes=tg, full.prop=prop, prop=prop2col, best.trans.gene=best.trans, best.tf=tfs.in.string[which.max(tf.weight)], g=rg)
  }
  return(prio)
}



#' Gets ranges in one object close by a set of other ranges
#'
#' Gets the first ranges in subject, which are up-/down-stream and overlapping
#' a range in the query
#'
#' @param query ranges for which to get nearby genes
#' @param subject ranges in which to look for nearby genes
#'
#' @return Either the idx of the hits in the subject if idxs=T, or the identified ranges (idxs=F)
#'
get.nearby.ranges <- function(query, subject) {

  nearby <- function(q,s){
    # get preceding, following and ovberlapping instances of any range in query within subject ranges
    pre <- precede(q, s, select="all", ignore.strand=T)
    fol <- follow(q, s, select="all", ignore.strand=T)
    ove <- findOverlaps(q, s, select="all", ignore.strand=T)
    # combine hits
    h <- unique(c(subjectHits(pre), subjectHits(fol), subjectHits(ove)))

    return(list(hits=h, ranges=s[h]))
  }
  res <- lapply(query, function(q) {
    n <- nearby(q, subject)
    n$ranges$distance <- rep(-1, times=length(n$ranges))
    for(i in 1:length(n$ranges)) {
      d <- distance(q,n$ranges[i])
      n$ranges[i]$distance <- d
    }
    return(n$ranges)
  })
  return(res)
}

#' Get human gene annotation with symbols.
#'
#' This will get the gene annotation annotated with respective gene SYMBOL from UCSC for the hg19 build
#' TODO: allow other builds as well, currently this method is not very flexible...
#'
#' @return GRanges object, where names(obj) are entrez-ids, and obj$SYMBOLS contains respective gene symbols
#'
get.gene.annotation <- function(drop.nas=TRUE) {
  library(Homo.sapiens)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  ucsc2symbol = AnnotationDbi::select(Homo.sapiens, keys=keys(Homo.sapiens, keytype="GENEID"),
                       keytype="GENEID", columns="SYMBOL")
  ga = genes(txdb)
  ga$SYMBOL <- ucsc2symbol[match(values(ga)[,"gene_id"],
                                              ucsc2symbol[,"GENEID"]),"SYMBOL"]
  if(drop.nas){
    ga <- ga[!is.na(ga$SYMBOL)]
  }
  return(ga)
}


#' Generate random graphs with matched CpGs instead of the meQTL CpGs
#'
#' This function builds a graph from the largest connected component in string.
#' It uses the transcription factor binding site data from the global tfbs.ann
#' data.frame to connect random CpGs (generated with the get.matched.cpgs
#' fuction instead of the actual trans meQTL CpGs.
#'
#' @param sentinel character indicating which locus (identifed by the sentinel
#'                 SNP) to process
#' @param joined.graph graphNEL object that has the full network info (see
#'                     network construction same also for the following)
#' @param locus.graphs list of graphNEL objects for each locus (should be named
#'                     with the sentinel SNP names
#' @param largest.cc.string graphNEL object with the largest connected
#'                          component of the STRING network
#' @param msd data.frame with the mean and standard deviation of methylation
#'            values for each CpG site
#' @param n.resample either a single integer giving the number of resamplings
#'                   or an integer vector of indices to use for the resampling
#'                   (this is useful if you want to add further resamplings
#'                   later, as the index is also used in naming)
#' @param tf.nodes character vector of node names of transcription factors
#' @param tfbs.ann logical indicator matrix with CpGs in the rows and
#'        transcription factors as columns
#' @param seed optional random seed
#'
#' @return an environment containing the radomized joined.graph and locus.graphs
#'
#' @author Matthias Heinig
#'
#' @export
random.graph.with.matched.cpgs <- function(sentinel, joined.graph, locus.graphs, largest.cc.string, msd, n.resample, tf.nodes, tfbs.ann, seed=NULL) {
  require(graph)
  require(igraph)

  ## we save the randomized locus graphs in an extra environment (not to
  ## interfere with the actual locus.graphs)
  bg.env = new.env()
  assign("locus.graphs", list(), envir = bg.env)
  jg = largest.cc.string

  ## get the trans genes in the graph
  trans.genes = intersect(adj(joined.graph, sentinel)[[1]], graph::nodes(largest.cc.string))

  ## get the cpgs in the graph
  g = locus.graphs[[sentinel]]
  cpgs = nodes(g)
  cpgs = cpgs[unlist(nodeData(g, cpgs, "cpg"))]
  if (length(cpgs) < 5) {
    cat("not enough CpGs\n")
    return(NULL)
  }
  lg = graphNEL(nodes=tf.nodes)

  ## get a matching background cpgs for each cpg with meQTL
  matched = lapply(cpgs, get.matched.cpgs, msd)
  if (is.null(seed)) {
    set.seed(match(sentinel, names(locus.graphs)))
  }

  ## generate graphs for the background
  if (length(n.resample) == 1) {
    n.resample = 1:n.resample
  }
  for (i in n.resample) {
    graph.ok = FALSE
    while (!graph.ok) {
      name = paste(sentinel, i, sep=".")
      bg.set = rownames(msd)[sapply(matched, sample, 1)]

      ## add to the string graph to produce a locus graph and a joined one
      graphs = add.to.graphs(graphs=list(lg, jg), name, trans.genes, bg.set, tfbs.ann)
      bg.g = graphs[[1]]
      jg = graphs[[2]]

      ## check the graph if at least one cpg is connected to a TF
      if (length(unlist(adj(jg, bg.set))) > 0) {
        graph.ok = TRUE
      }
    }
    with(bg.env, {locus.graphs[[name]] <- bg.g})

    ## call the garbage collector to run smoothly on the big compute nodes
    gc()
  }

  ## remove singletons
  jg = subGraph(nodes(jg)[graph::degree(jg) > 0], jg)

  assign("joined.graph", jg, envir = bg.env)

  ## just to make sure we have one connected component (this should be the case
  ## anyway, but just make sure)
  ig = graph_from_graphnel(jg)
  cl = clusters(ig)
  keep = nodes(jg)[cl$membership == which.max(cl$csize)]
  assign("largest.cc", subGraph(keep, jg), envir = bg.env)

  return(bg.env)

}

#'
#' Define a convenience function to get linear model p-values
#'
#' @param modelobject The linear model for which to get the p-value
#'
#' @return The p-value to the given linear model
#'
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


#' Extension to the waitForJobs function of the BatchJobs package which shows
#' some strange behaviour when waiting for jobs (database locked)
#' so we need to make it extra failsafe.
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param reg BatchJobs registry
#' @param waittime time to wait before updating job status
#' @param nretry number of time to retry getting the job status before throwing
#'        an error
#'
#' @imports BatchJobs
myWaitForJobs <- function(reg, waittime=3, nretry=100) {
  require(BatchJobs)
  success = FALSE
  while (nretry > 0 && !success) {
    status = tryCatch({
      while (TRUE) {
        status = showStatus(reg)
        if (status$done + status$expired == status$n) {
          cat("done\n")
          return(list(success=TRUE, nretry=nretry))
        }
        Sys.sleep(waittime)
      }
      return(list(success=FALSE, nretry=nretry))
    }, error=function(e) {
      cat("Error while waiting for jobs:\n")
      print(e)
      cat("\nnumber of retries left: ", nretry - 1, "\n")
      Sys.sleep(waittime + runif(1, 0, 3))
      return(list(success=FALSE, nretry=nretry - 1))
    })
    success = status$success
    nretry = status$nretry
    cat("success after the tryCatch block:", success, "\n")
    cat("nretry after the tryCatch block:", nretry, "\n")
  }


  if (!success) {
    err.msg = paste("Error during batch processing in registry")
    save(envir=sys.frame(), list=ls(envir=sys.frame()), file=file.path(dir, "error_image.RData"))
    stop(err.msg)
  }
}


#' helper function to run on the grid and clean up afterwards
#'
#' @param rfun function to run on the grid the first argument will be iterated
#' @param idx these are the values for the first argument to be iterated through
#' @param more.args list of named arguments for rfun should be list() if none
#' @param name a name for the batch process
#' @param dir directory where to store everything for job handling
#' @param clean.up boolean indicating whether to delete intermediate files
#' @param resources list with resources requested for the submitJobs function
#' @param n.chunks integer indicating the number of chunks. This is useful when
#'                 many jobs are submitted to avoid heavy load on the scheduler
#'                 as well as on the job registry. NULL: no chunks.
#'
#' @return list of results of rfun for every element of idx
#'
#' @imports BatchJobs
#' @export
run.batchjobs <- function(rfun, idx, more.args, name, dir, clean.up=TRUE, resources=list(), n.chunks=NULL, reuse.registry=FALSE) {
  require(BatchJobs)
  ## wrap another run function around the original one to source R/lib.R
  ## also call the garbage collector to avoid memory problems when processing
  ## chunks
  rrfun <- function(...) {
    source("R/lib.R")
    res = rfun(...)
    gc()
    return(res)
  }

  if (reuse.registry) {
    load(paste0(dir, "/registry.RData"))
  } else {
    reg = makeRegistry(name, file.dir=dir)
  }
  batchMap(reg, rrfun, idx, more.args=more.args)

  ## jobs can be packed together in chunks if there are too many
  chunked = getJobIds(reg)
  if (!is.null(n.chunks)) {
    chunked = chunk(chunked, n.chunks=n.chunks, shuffle = TRUE)
  }

  ## job delay to prevent concurrent access to the database by too many jobs
  submitJobs(reg, ids=chunked, resources=resources, job.delay=function(n, i) runif(1, 0, 0.1 * n))
  Sys.sleep(20)

  ## also a custom wait function that is more error tolerant with many jobs
  myWaitForJobs(reg, waittime=30, nretry=100)
  res = reduceResultsList(reg)
  if (clean.up) {
    removeRegistry(reg, ask="no")
  }
  return(res)
}

#' -----------------------------------------------------------------------------
#' helper function to run on the grid and clean up afterwards
#'
#' @param rfun function to run on the grid the first argument will be iterated
#' @param idx these are the values for the first argument to be iterated through
#' @param more.args list of named arguments for rfun should be list() if none
#' @param name a name for the batch process
#' @param dir directory where to store everything for job handling
#' @param clean.up boolean indicating whether to delete intermediate files
#' @param resources list with resources requested for the submitJobs function
#' @param n.chunks integer indicating the number of chunks. This is useful when
#'                 many jobs are submitted to avoid heavy load on the scheduler
#'                 as well as on the job registry. NULL: no chunks.
#' @param source List of files which should be sourced within each job
#'
#' @return list of results of rfun for every element of idx
#'
#' @imports batchtools
#' @export
#' -----------------------------------------------------------------------------
run.batchtools <- function(rfun, idx, more.args, name, dir,
                           clean.up=TRUE, resources=list(),
                           n.chunks=NULL, reuse.registry=FALSE,
                           source=character(0L)) {
  require(batchtools)

  reg = makeRegistry(file.dir=dir, source=source)

  ids <- batchMap(rfun, idx, more.args=more.args)

  ## jobs can be packed together in chunks if there are too many
  if (!is.null(n.chunks)) {
    ids[, chunk := chunk(job.id, n.chunks=n.chunks)]
  }

  submitJobs(ids=ids, resources=resources)

  all_done <- waitForJobs()
  if(!all_done)  {
    resub <- function() {
      print("Not all jobs successful! Tying once more.")

      # jobs which are not yet submitted or not terminated successfully
      ids <- ujoin(findNotDone(), findNotStarted(), all.y=T)

      # resubmit
      submitJobs(ids, resources = resources)
      all_done <- waitForJobs()
      return(all_done)
    }
    all_done <- resub()
    if(!all_done) {
      # once more
      all_done <- resub()
      if(!all_done) warning("Not all jobs successful!")
    }
  }

  res = reduceResultsList()
  if (clean.up & all_done) {
#    removeRegistry(wait=0)
  }
  return(res)
}

#' Evaluate permutations on the grid
#'
#' @param idx run index (iterated by batchMap)
#' @param runs data.frame with two columns: 'sentinel' and 'permutation'
#' @param prefix character indicating the path prefix for the files
#' @param redo logical indicating whether to rerun even if a results file exists
#' @param subset list for subsetting the trans genes passed to
#'               prioritization.table (see the documentation there)
#'
#' @return list with the best candidate for each permutation run
#' @export
evaluate.perm <- function(idx, runs, prefix, redo, subset=NULL) {

  ## get the parameters of the run
  sentinel = runs[idx, "sentinel"]
  permutation = runs[idx, "permutation"]
  name = paste(sentinel, permutation, sep=".")

  outname = paste0(prefix, "_background_cpgs/", sentinel, "/", name, ".RData")

  if (!file.exists(outname) || redo) {
    gfile = paste0(prefix, "_background_cpgs/", sentinel, ".RData")
    load(gfile)

    locus.graph = locus.graphs[[name]]

    prio = prioritze.by.randomwalk(locus.graph, largest.cc=largest.cc, string.nodes=string.nodes, tf.nodes=tf.nodes, n.eigs=500, mode="propagation", return.g=FALSE)

    ## save the full result
    dir.create(dirname(outname), recursive=T)
    save(prio, file=outname)
  } else {
    load(outname)
  }

  ## return only the summary
  prio = list(prio)
  names(prio) = name
  return(prioritization.table(prio, subset=subset))
}

evaluate.quick.perm <- function(idx, runs, prefix, redo, largest.cc.string, subset=NULL) {

  ## get the parameters of the run
  sentinel = runs[idx, "sentinel"]

  outname = paste0(prefix, "_background_cpgs/", sentinel, "/quick_perm.RData")

  if (!file.exists(outname) || redo) {
    gfile = paste0(prefix, "_background_cpgs/", sentinel, ".RData")
    load(gfile)

    prio = prioritize.quick(locus.graphs, largest.cc.string=largest.cc.string,  n.eigs=500, return.g=FALSE)

    ## save the full result
    dir.create(dirname(outname), recursive=T)
    save(prio, file=outname)
  } else {
    load(outname)
  }

  ## return only the summary
  return(prioritization.table(prio, subset=subset))
}

prioritization.empirical.pvalues <- function(tab, full.tab, permuted.prioritization, runs, mode="propagation") {
  tab = data.frame(tab, empirical.p=NA, n.perm=NA)
  full.tab = data.frame(full.tab, empirical.p=NA, n.perm=NA)

  if (mode == "propagation") {
    select.fun <- function(x) return(1)
  } else {
    select.fun <- function(x) return(1:nrow(x))
  }

  sentinels = unique(runs[,"sentinel"])
  for (sentinel in sentinels) {
    ## in case some jobs are missing we use the run id (character)
    run.idx = as.character(which(runs[,"sentinel"] == sentinel))
    missing = setdiff(run.idx, names(permuted.prioritization))
    if (length(missing) > 0) {
      cat("Missing runs for ", sentinel ,":", length(missing), "\n")
      run.idx = setdiff(run.idx, missing)
    }
    bg = lapply(permuted.prioritization[run.idx], function(x) if (!is.null(x)) return(x[select.fun(x),"best.trans.prop"]) else return(NULL))
    if (all(sapply(bg, length) == 0)) {
      cat("all permutations null for sentinel", sentinel, "\n")
      next
    }

    bg.ecdf = ecdf(unlist(bg))

    n.perm = length(unlist(bg))

    idx = tab[,"sentinel"] == sentinel
    tab[idx, "empirical.p"] = 1 - bg.ecdf(tab[idx, "best.trans.prop"])
    tab[idx, "n.perm"] = n.perm

    idx = full.tab[,"sentinel"] == sentinel
    full.tab[idx, "empirical.p"] = 1 - bg.ecdf(full.tab[idx, "prop"])
    full.tab[idx, "n.perm"] = n.perm
  }
  return(list(tab=tab, full.tab=full.tab))
}



tfbs.heatmap <- function(res, qcol, threshold, prune, max.or=Inf, min.or=0) {
  ## make a nice heatmap
  library(ggplot2)
  library(reshape)
  library(scales)

  ## select all experiments with at least one significant enrichment
  max.or = min(c(max(res$test.estimate[res$test.estimate < Inf]), max.or))
  min.or = max(c(min(res$test.estimate[res$test.estimate > 0]), min.or))
  selected = res[res[,qcol] < threshold,]
  selected$test.estimate[selected$test.estimate > max.or] = max.or
  selected$test.estimate[selected$test.estimate < min.or] = min.or

  ## order loci by similar tf profiles
  selected.mat = cast(res, sentinel ~ tf, value=qcol)
  rownames(selected.mat) = selected.mat[,1]
  selected.mat = selected.mat[,-1]
  selected.mat = selected.mat < threshold

  ord = hclust(dist(as.matrix(selected.mat)))$order
  selected$sentinel = factor(as.character(selected$sentinel), levels=rownames(selected.mat)[ord])
  ## reverse the tf levels to get top down alphabetic order
  ## selected$tf = factor(as.character(selected$tf), rev(unique(sort(as.character(selected$tf)))))

  ## order tfs by similar locus profiles
  ord = hclust(dist(t(as.matrix(selected.mat))))$order
  selected$tf = factor(as.character(selected$tf), levels=colnames(selected.mat)[ord])

  ## also try a faceting by cell line / condition
  condition = sapply(strsplit(as.character(selected$tf), ".", fixed=T), "[", 2)
  selected = cbind(selected, condition)

  if (!prune) {

    return(qplot(sentinel, tf, geom="tile", fill=test.estimate, data=selected, ylab="Transcription factor", xlab="Locus") + scale_fill_gradient2(low="blue", mid="white", high="green", trans="log", limits=c(min.or, max.or), guide=guide_colourbar(title="Odds ratio")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8)) + theme(axis.text.y = element_text(size=8))) # + facet_grid(. ~ condition, scales="free_x", space="free_x")

  } else {
    ## also make a pruned list for a clearer picture
    keep = names(which(table(selected$tf) > 0))
    selected$tf = factor(selected$tf, levels=levels(selected$tf)[levels(selected$tf) %in% keep])
    keep = names(which(table(selected$sentinel) > 0))
    selected$sentinel = factor(selected$sentinel, levels=levels(selected$sentinel)[levels(selected$sentinel) %in% keep])

    return(qplot(sentinel, tf, geom="tile", fill=test.estimate, data=selected, ylab="Transcription factor", xlab="Locus") + scale_fill_gradient2(low="blue", mid="white", high="green", trans="log", limits=c(min.or, max.or), guide=guide_colourbar(title="Odds ratio")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8)) + theme(axis.text.y = element_text(size=8))) # + facet_grid(. ~ condition, scales="free_x", space="free_x")
  }
}

#' Gets a single snp range from the gtex snp database
#'
#' @author Johann Hawe
#'
get.snp.range <- function(snp){
  library(data.table)
  snps <- fread(paste0("results/current/gtex-snp-locations.txt"))
  snps <- snps[RS_ID_dbSNP142_CHG37p13 %in% snp,]
  if(nrow(snps)<1){
    warning("SNP not found on gtex db.\n")
    return(NULL)
  } else {
    r <- with(snps,
              GRanges(paste0("chr", Chr),
                      IRanges(as.integer(Pos),width=1)))
    names(r) <- snps$RS_ID_dbSNP142_CHG37p13
    return(r)
  }
}

#' Gets the chip-seq context for a specific set of cpgs
#'
#' @param cpgs ID of cpgs for which to get the context
#'
#' @author Johann Hawe
#'
get.chipseq.context <- function(cpgs){
  load("results/current/cpgs_with_chipseq_context_100.RData")
  tfbs.ann <- tfbs.ann[rownames(tfbs.ann) %in% cpgs,,drop=F]

  return(tfbs.ann[cpgs,,drop=F])
}

#' Creates a graphNEL object from a given bdgraph result for a defined cutoff
#'
#' @paran result The bdgraph result
#' @param cutoff Cutoff to be used for posterior probability of edges
#'
#' @value graph-nel object created from the bdgraph result
#'
#' @author Johann Hawe
#'
graphNEL.from.result <- function(result, cutoff){
  library(BDgraph)
  library(graph)
  library(igraph)
  g.adj <- BDgraph::select(result, cut = cutoff)
  graph.bd <- as_graphnel(graph.adjacency(g.adj, mode="undirected", diag=F))
  return(graph.bd)
}

#' Plot a GGM result graph
#'
#' Plots the a built graph (estimated from the sentinel data) using the
#' twopi-visualization. Can be used to retrieve only the plot graph/node/edge
#' attributes when using pdf.out=NULL, dot.out=NULL and plot.on.device=F
#'
#' @param graph: Graph to be plotted (graphNEL)
#'
#' @param id The Id of the sentinel snp
#' @param ranges The list of sentinel associated ranges (we need snp.genes
#' and cpg.genes, enriched.tfs etc.)
#' @param dot.out File to which to write the graph in dot format
#'
#' @return List of plot graph attributes and the underlying graph structure
#'
#' @author Johann Hawe
#'
plot.ggm <- function(g, id, dot.out=NULL){
  library(graph)
  library(Rgraphviz)
  library(GenomicRanges)

  # remove any unconnected nodes (primarily for bdgraph result, since
  # such nodes are already removed for genenet)
  if(any(graph::degree(g) == 0)){
    g <- removeNode(names(which(graph::degree(g) == 0)), g)
  }

  # add sentinel to network if its is not in there yet or it has been removed
  if(!(id %in% nodes(g))) {
    g <- graph::addNode(c(id), g)
  }

  n <- nodes(g)

  # get trans and cpg gene symbols
  snp.genes <- n[unlist(nodeData(g,n,"snp.gene"))]
  cpg.genes <- n[unlist(nodeData(g,n,"cpg.gene"))]
  tfs <- n[unlist(nodeData(g,n,"tf"))]

  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=14,
                          style="filled", fontname="helvetica"),
                graph=list(overlap="false", root=id, outputorder="edgesfirst"))

  shape = rep("ellipse", numNodes(g))
  names(shape) = n
  shape[grep("^cg", n)] = "box"
  shape[grep("^rs", n)] = "box"

  width = rep(0.8, numNodes(g))
  names(width) = n
  width[grep("cg", n)] = 0.4

  height = rep(0.3, numNodes(g))
  names(height) = n
  height[grep("cg", n)] = 0.4

  label = n
  names(label) = n
  label[grep("cg", n)] = ""

  col = rep("#ffffff", numNodes(g))
  names(col) = n
  col[grep("^rs", n)] = "#fab4ad";
  col[grep("^cg", n)] = "#e4d7bc";
  if(!is.null(tfs)){
    col[tfs] = "green"
  }

  penwidth = rep(1, numNodes(g))
  names(penwidth) = n
  penwidth[snp.genes] = 3
  penwidth[cpg.genes] = 3
  if(!is.null(tfs)){
    penwidth[tfs] = 3
  }

  bordercol = rep("black", numNodes(g));
  names(bordercol) = n;
  bordercol[cpg.genes] = "#e4d7bc";
  bordercol[id] = "#fab4ad";

  nAttrs = list(shape=shape, label=label, width=width,
                height=height, penwidth=penwidth, fillcolor=col,
                color=bordercol)

  # default color for edges: black
  ecol = rep("black", numEdges(g))
  names(ecol) = edgeNames(g)
  for(i in snp.genes) {
    # color any edge from a SNP to one of our snp genes red
    ecol[grepl("^rs|~rs", names(ecol)) & grepl(i, names(ecol))] = "#b3cde2"
  }

  # set also color for cpgs
  for(cg in cpg.genes){
    # color any edge from a cpg to one of its cpg genes blue (proximity edges)
    ecol[grepl("^cg|~cg", names(ecol)) & grepl(cg, names(ecol))] = "#b3cde2"
  }

  # check edgeData and add to colors
  for(edge in names(ecol)){
    n1 <- strsplit(edge,"~")[[1]][1]
    n2 <- strsplit(edge,"~")[[1]][2]

    if(unlist(graph::edgeData(g,n1,n2, "isPPI"))){
      ecol[edge] <- "#decae3"
    }
    if(unlist(graph::edgeData(g,n1,n2, "isChipSeq"))){
      ecol[edge] <- "#ccebc5"
    }
  }

  dir = rep("none", numEdges(g))
  names(dir) = edgeNames(g)

  eAttrs = list(color=ecol, dir=dir)

  if(numEdges(g)>500){
    warning("Skipping plotting on device due to large amount of edges")
  } else{
    plot(g, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }

  if(!is.null(dot.out)){
    # output the dot-information
    #    toDot(graph, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs, subGList=sgs)
    toDot(g, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }

  # return the list of created plotattributes and the possibly modified graph
  # object
  return(list(graph=g, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs))
}

#' Get  STRING interactions (only experimental and db supported)
#'
#' Loads the STRING db saved on disc to a graph R object.
#' Uses a caching mechanism for faster loading of the graph information
#'
#' @return nothing
#'
load.string.db <- function(reload=F) {
  library(data.table)
  library(igraph)
  library(graph)

  cat("Loading string db.", "\n")

  fcache <- "results/current/string.v9.expr.RData"
  if(reload || !file.exists(fcache)){
  # load db anew
  string.all <- fread(paste0("data/current/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt"),
                        data.table=F, header=T, stringsAsFactors=F)
  string.inter <- string.all[string.all$experimental>=1 | string.all$database>=1,]
  rm(string.all)

  string.nodes <- unique(c(string.inter[,1], string.inter[,2]))
  string.db <- graphNEL(nodes=string.nodes)
  string.db <- addEdge(string.inter[,1],
                       string.inter[,2],
                       string.db)

  expr = read.csv("data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",
                        sep="\t",
                        skip=2,
                        stringsAsFactors=F)

  expressed = expr[expr[,"Whole.Blood"] > 0.1, "Name"]
  expressed = sapply(strsplit(expressed, "." ,fixed=T) ,"[", 1)
  library(Homo.sapiens)
  expressed.symbols = select(Homo.sapiens, keys=expressed, keytype="ENSEMBL", columns="SYMBOL")
  expressed.symbols = unique(expressed.symbols[,"SYMBOL"])

  string.nodes = intersect(nodes(string.db), expressed.symbols)
  string.db = subGraph(string.nodes, string.db)

  # get largest connected component
  ig = graph_from_graphnel(string.db)
  cl = clusters(ig)
  keep = nodes(string.db)[cl$membership == which.max(cl$csize)]
  string.db = subGraph(keep, string.db)

  save(file=fcache, string.db)
  } else {
    load(fcache)
  }
  # set to global environment
  assign("STRING.DB", string.db, envir=.GlobalEnv)

  cat("Done.\n")
}

#' Gets residuals based on linear models from a data matrix
#'
#' Calculates the residual matrix for a given matrix considering available covariates.
#' Uses linear model for methylation data and linear mixed
#' model with plate as random effect for the expression data.
#'
#' @param data the matrix for which to calculate the covariates
#' @param data.type the type of data in the matrix: either "meth" or "expr". Depending on the type
#' different formulas are used for the linear model.
#' @param col.names The col.names over which to iterate in the dataframe to calculate
#' the residuals on (e.g. probe.ides, gene.names,..)
#' @param use.cohort Flag whether to include a cohort variable in the model (expression only)
#'
#' @return A matrix  where in the colums are the measured entities (e.g.probes)
#' and in the rows are the samples, containing the calculated residuals
#'
get.residuals <- function(data, data.type,
                          col.names=NULL,
                          use.cohort=F,
                          threads=1) {
  cols <- col.names

  if(data.type=="meth") {
    if(is.null(cols)) {
      cols <- colnames(data)[grepl("^cg", colnames(data))]
    }
    res <- mclapply(cols, function(n) {
      fm <- as.formula(paste0(n,"~",
                             paste0("1+CD4T+CD8T+NK+Bcell+Mono+",
                                       paste0("PC",
                                              paste(1:20, "cp", sep="_"),
                                              collapse="+"))))
      return(lm(fm, data=data, na.action=na.exclude))
    }, mc.cores=threads)
  } else if(data.type == "expr") {
    if(is.null(cols)){
      cols <- colnames(data)[grepl("^ILMN_", colnames(data))]
    }
    res <- mclapply(cols, function(n) {
      if(use.cohort){
        fm <- as.formula(paste0(n, "~",
                           "1+age+sex+RIN+batch1+batch2+cohort"))
      } else {
        fm <- as.formula(paste0(n, "~",
                           "1+age+sex+RIN+batch1+batch2"))
      }
      return(lm(fm,data=data, na.action=na.exclude))
    }, mc.cores=threads)
  } else {
    stop("Data type not supported for residual calculation.")
  }
  # build the full residual matrix from model results
  residual.mat <- matrix(data=unlist(lapply(res, resid)), nrow=nrow(data))
  colnames(residual.mat) <- cols
  rownames(residual.mat) <- rownames(data)

  return(residual.mat)
}

#' Method to quickly filter an edge matrix for only those edges, which are within
#' a specified graph
#'
#' @author Johann Hawe
#'
#' @date 2017/03/13
#'
filter.edge.matrix <- function(g, em){
  e <- graph::edges(g)
  out <- matrix(ncol=2,nrow=0)
  for(i in 1:nrow(em)){
    e1 <- em[i,1]
    e2 <- em[i,2]
    if(e1 %in% names(e)){
      if(e2 %in% e[[e1]]){
        out <- rbind(out,c(e1,e2))
      }
    }
  }
  return(out)
}

#' Quantile normalization
#'
#' @param x ngenes x nsamples matrix to be normalized
#' @return quantile normalized matrix
#' @export
normalize.quantile <- function(x) {
  x = as.matrix(x)
  o = apply(x, 2, order)
  xsort = x
  for (col in 1:ncol(x)) {
    xsort[,col] = x[o[,col], col]
  }
  means = apply(xsort, 1, mean)
  normalized = matrix(NA, nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
  for (col in 1:ncol(x)) {
    normalized[o[,col], col] = means
  }
  return(normalized)
}


#' Annotate positions with their epigenetic states
#'
#' @param cpg.ranges GRanges object with the positions to annotate
#' @param ids character vector with the roadmap epigenome ids to use
#' @param dir the directory where chromHMM files are stored
#'
#' @return character matrix with length(cpg.ranges) rows and length(ids) columns
#'         with the state annotation of the range in each epigenome
#'
#' @export
chromHMM.annotation <- function(cpg.ranges,
                                ids,
                                dir="data/current/roadmap/chromHMM/15state/",
                                suffix="_15_coreMarks_mnemonics.bed.bgz") {
  require(Rsamtools)

  annotation = sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(cpg.ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=cpg.ranges[avail])

    avail.ann = unlist(lapply(avail.ann, function(x) {
      strsplit(x, "\t")[[1]][4]}))
    ann = rep(NA, length(cpg.ranges))
    ann[avail] = avail.ann
    return(ann)
  })

  colnames(annotation) = ids
  rownames(annotation) = names(cpg.ranges)

  return(annotation)
}

# ------------------------------------------------------------------------------
# Get chromHMM annotations. Tuned to work with a large number of GenomicRanges.
# ------------------------------------------------------------------------------
get_chromHMM_annotation <- function(ranges,
                                    ids,
                                    dir="data/current/roadmap/chromHMM/15state/",
                                    suffix="_15_coreMarks_mnemonics.bed.bgz") {
  require(Rsamtools)
  require(data.table)

  annotation = sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(ranges) %in% seqnamesTabix(file))
    states <- fread(cmd=paste0("zcat ", file), header=F)
    states <- with(states, GRanges(V1, IRanges(V2, V3), state=V4))

    st <- states[findOverlaps(ranges[avail], states, select = "first")]$state

    ann = rep(NA, length(ranges))
    ann[avail] = st

    return(ann)
  })

  colnames(annotation) = ids
  rownames(annotation) = names(ranges)

  return(annotation)
}

#' For a set of symbols, gets the emsembl-ids of the corresponding genes from
#' the AnnotationHub.
#' Note: Currently only implemented for human genes
#'
#' @return A data.frame containng the SYMBOL in one column and the ensemble
#' id in another column
#'
#' @author Johann Hawe
#'
#' @date 2017/03/27
#'
ensembl.from.symbol <- function(symbols, na.drop=T, org="hs"){
  result <- NULL
  if("hs" %in% org){
    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    result <- select(hs,
                     keytype="SYMBOL",
                     keys=c(symbols),
                     columns=c("SYMBOL", "ENSEMBL"))
  } else {
    warning("Organism not supported yet: ", org)
  }
  if(na.drop){
    result <- result[!is.na(result$ENSEMBL),,drop=F]
  }

  return(result)
}

#' For a set of ensemble ids, gets the symbols of the corresponding genes from
#' the AnnotationHub.
#' Note: Currently only implemented for human genes
#'
#' @return A data.frame containng the SYMBOL in one column and the ENSEMBLE
#' id in another column
#'
#' @author Johann Hawe
#'
#' @date 2017/03/28
#'
symbol.from.ensembl <- function(ensemblIDs, na.drop=T, org="hs"){
  result <- NULL
  if("hs" %in% org){
    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    result <- select(hs,
                     keytype="ENSEMBL",
                     keys=c(ensemblIDs),
                     columns=c("SYMBOL", "ENSEMBL"))
  } else {
    warning("Organism not supported yet: ", org)
  }
  if(na.drop){
    result <- result[!is.na(result$ENSEMBL),,drop=F]
  }
  # only get the first match for each id
  return(result[!duplicated(result$ENSEMBL),])
}

#' Plot a simple heatmap
#'
#' Uses 'pheatmap' package to plot a very simple heatmap.
#'
#' @param x The data matrix for which to plot the heatmap
#' @param cannot vector of column annotations (must have same names as colnames
#' of x)
#' @param cluster Flag whether to cluster the data in the heatmap or not. default:F
#' @param cors Flag whether correlations are contained in the input matrix x (to
#' appropriately choose the color palette)
#'
#' @return Grob object returned by pheatmap in the 4th position
#'
simple.heatmap <- function(x, cors=F, cannot=NA, cluster=F) {
  library(pheatmap)
  library(RColorBrewer)

  if(cors){
    breaksList = seq(-1, 1, by = 0.1)
  } else {
    breaksList = seq(min(x), max(x), by = 0.1)
  }

  cols <- rev(brewer.pal(n = 7, name = "RdYlBu"))
  cols <- colorRampPalette(cols)(length(breaksList))

  return(pheatmap(x, annotation_col=cannot,
           cluster_rows=cluster,
           cluster_cols=cluster,
           cex=0.7, color = cols, breaks = breaksList)[[4]])
}

# -----------------------------------------------------------------------------
#' Scans genotype files for SNPs within the provided genomic ranges
#'
#' @param ranges GRanges object containing the ranges which to scan for SNPs
#' @param dir Base directory in which genotype information is stored
#' @param genotype.file Path to and including file which contains the genotypes,
#' relative to the base directory
#' @param id.file Path to and including file which contains the individual ids,
#' relative to the base directory
#'
#' @param filter Vector of individual ids to filter for
#' 
#' @author Johann Hawe
#'
# -----------------------------------------------------------------------------
scan.snps <- function(ranges,
                      dir="/storage/groups/epigenereg01/projects/PV_K14115g_Heinig/results/20160204/genoF4",
                      genotype.file="dosage_combined/MAF001/full_sorted.bgz",
                      id.file="individuals.txt",
                      filter = NULL) {
  require(Rsamtools)
  require(data.table)
  
  # load the genotype individuals
  individuals <- as.character(fread(file.path(dir, id.file))$V1)
  
  res <- scanTabix(file.path(dir, genotype.file), 
                   param=ranges)
  
  # scanTabix will also return 'empty' ranges for which no data was found
  res_not_empty <- res[sapply(res, function(x) !length(x)<1)]
  
  print("Scan done, converting Tabix results to data frame.")
  
  if(length(res_not_empty) < 2 && !length(res_not_empty[[1]]) > 1) {
    collapsed <- paste0(unlist(res_not_empty), "\n")
  } else {
    collapsed <- paste0(unlist(res_not_empty), collapse="\n")
  }
  
  # define columns for filtering
  all_cols <- c("chr", "name", "pos", "orig", "alt", individuals)
  if(!is.null(filter)) {
    filter_cols <- c("chr", "name", "pos", "orig", "alt", filter)
  } else {
    filter_cols <- all_cols
  }
  
  select <- match(filter_cols,
                  all_cols)
  
  # select argument also determines order of columns
  genos <- unique(fread(text=collapsed, sep="\t", header=FALSE, 
                        col.names=filter_cols,
                        stringsAsFactors=F, data.table=FALSE, 
                        select=select))
  
  print("Read data of dimension:")
  print(dim(genos))
  
  # set rownames
  rownames(genos) <- genos$name

  return(list(snpInfo=genos[,c(2,1,3,4,5)], snps=genos[,6:ncol(genos)]))
}

#'
#' Collects methylation, expression and genotype data for a either KORA or LOLIPOP
#' cohort.
#'
#' @param cohort The cohort data to load. Either "KORA" or "LOLIPOP".
#' @param snp.ids SNP ids for which to get the genotype data (for LOLIPOP)
#' @param snp.ranges SNP ranges for which to get the genotype data (for KORA)
#' @param meth.probes Optional list of methylation array probe ids to be retrieved
#' @param expr.probes Optional list of expression array probe ids to be retrieved
#' @param cache.global Flag whether to cache the expression and methylation data
#' in the global environment and reuse them when calling the method again.
#'
#' @return A data matrix containing methylation, expression and genotype data
#' for a specific sentinel
#'
#' @autho Johann Hawe
#'
#' @date 02/06/2017
#'
#' @export
#'
collect_data <- function(cohort=c("kora", "lolipop"), snp.ranges=NULL,
                         snp.ids=NULL, meth.probes=NULL, expr.probes=NULL,
                         cache.global=F) {

  #TODO here we use some old/hardcoded paths -> move this to SNAKEMAKE pipeline
  # based paths
  if("kora" %in% cohort) {
    if(is.null(snp.ranges)){
      stop("SNP ranges need to be supplied when loading KORA data.")
    }
    cat("Loading KORA data.\n")

    # load id mapping to intersect expr/meth/geno data
    ID.MAP <- read.table("~/work/analysis/PV_K14115g_Heinig/data/current/individuals_covariates.csv",
                         header=T, sep=";", stringsAsFactors = F)
    # convert ids to character (instead of int) for easier subsetting of data matrices
    ID.MAP$axio_s4f4 <- as.character(ID.MAP$axio_s4f4)
    ID.MAP$expr_s4f4ogtt <- as.character(ID.MAP$expr_s4f4ogtt)
    ID.MAP$meth_f4 <- as.character(ID.MAP$meth_f4)
    # drop individuals with NAs (e.g. due to missing BMI)
    ID.MAP <- na.omit(ID.MAP)
    # load dummy snp to get the individuals which have geno data available
    dummy <- GRanges("chr4", ranges=IRanges(156902056,width=1))
    dummy <- t(scan.snps(dummy)$snps)

    if(!exists("f4.norm") | !exists("beta")) {
      cat("Preparing KORA raw data.\n")
      load("data/current/kora/expression/kora_f4_normalized.Rdata")
      load("data/current/kora/methylation/KF4_beta_qn_bmiq.RData")

      if(cache.global){
        assign("f4.norm", f4.norm, .GlobalEnv)
        assign("beta", beta, .GlobalEnv)
      }
    }

    # gets as 687 individuals, having all data available (some ids of the id
    # map are not contained within the data frame...)
    toUse <- which(ID.MAP$expr_s4f4ogtt %in% colnames(f4.norm) &
                     ID.MAP$meth_f4 %in% colnames(beta) &
                     ID.MAP$axio_s4f4 %in% rownames(dummy))

    ID.MAP <- ID.MAP[toUse,]
    ID.MAP$utbmi <- NULL
    ID.MAP$ul_wbc <- NULL
    ID.MAP$utalteru <- NULL

    cat("Using ", nrow(ID.MAP), " samples.\n")

    # sort our input data s.t. each row corresponds to the same individual
    # using the created ID mapping table
    if(!is.na(meth.probes)){
      meth <- t(beta[,ID.MAP$meth_f4])
      if(!is.null(meth.probes)){
        meth <- meth[,meth.probes,drop=F]
      }
      # get rid of zero-variance probes
      meth <- meth[,apply(meth,2,var, na.rm=T)!=0, drop=F]
    }

    if(!is.na(expr.probes)){
      # use only those individuals for which we have all data available
      expr <- t(f4.norm[,ID.MAP$expr_s4f4ogtt])
      if(!is.null(expr.probes)){
        expr <- expr[,expr.probes,drop=F]
      }
      # remove zeor-variance probes...
      expr <- expr[,apply(expr,2,var)!=0, drop=F]
    }
    if(!is.na(snp.ranges)){
      geno <- get.genotypes(snp.ranges, snp.ids, cohort)
      geno <- geno[ID.MAP$axio_s4f4,,drop=F]
    }
    # get the methylation PCA results
    load(paste0("data/current/kora/methylation/control_probe_pcs_n1727.RData"))
    pcs <- pcs[ID.MAP$meth_f4,]
    colnames(pcs) <- paste(colnames(pcs), "cp", sep="_")
    meth.pcs <- pcs
    rm(pcs)

    # load technical covariates for expression data
    load(paste0("data/current/kora/expression/technical_covariables_kora_f4.Rdata"))
    rownames(covars.f4) <- covars.f4[,1]
    covars.f4 <- covars.f4[ID.MAP$expr_s4f4ogtt,c(2:6)]
    covars.f4$sex <- as.factor(covars.f4$sex)

    # load houseman blood count data for methylation
    houseman <- read.table("data/current/kora/methylation/Houseman/KF4_QN_estimated_cell_distribution_meanimpute473_lessThanOneTRUE.csv",
                           sep=";", header=T,row.names=1)
    houseman <- houseman[ID.MAP$meth_f4,]

    # create initial data frame with covariates
    covars.f4 <- cbind(as.data.frame(covars.f4), houseman, meth.pcs)

    # replace some colnames to easily match with the lolipop data
    cc <- colnames(covars.f4)
    cc[grepl("storage.time",cc)] <- "batch2"
    cc[grepl("plate",cc)] <- "batch1"
    colnames(covars.f4) <- cc
    covars.f4[,"batch1"] <- factor(covars.f4[,"batch1"])

    data <- cbind.data.frame(covars.f4, stringsAsFactors=F)
    if(!is.na(meth.probes)) {
      data <- cbind.data.frame(data, meth, stringsAsFactors=F)
    }
    if(!is.na(expr.probes)){
      data <- cbind.data.frame(data, expr, stringsAsFactors=F)
    }
    if(!is.na(snp.ranges)) {
      data <- cbind.data.frame(data, geno, geno.ids=rownames(geno),
                               stringsAsFactors=F)
    }

    return(data)

  } else if ("lolipop" %in% cohort) {
    if(is.null(snp.ids)){
      stop("SNP ids need to be supplied when loading LOLIPOP data.")
    }

    cat("Loading LOLIPOP data.\n")
    load("data/current/meQTLs/ggmdata_201017.RData")

    if(!is.na(snp.ids)){
      # currently we only use 1 of the 3 available sets of genotypes.
      # at some point we want to check on the individual ones for our analysis
      geno <- get.genotypes(snp.ranges, snp.ids, cohort)
    } else {
      # we need some genotype data to order our expr/meth data
      # this is however discarded in the end
      geno <- t(dosage.oX_ia)
    }

    if(!is.na(meth.probes)) {
      meth <- t(beta)
      # sort using the genotype ordering
      meth <- meth[rownames(geno),]
      if(!is.null(meth.probes)){
        meth <- meth[,meth.probes,drop=F]
      }
    }

    if(!is.na(expr.probes)) {
      expr <- t(exma)
      # sort using the genotype ordering
      expr <- expr[rownames(geno),]
      if(!is.null(expr.probes)){
         expr.probes <- unique(expr.probes)
         idx <- which(expr.probes %in% colnames(expr))
         if(length(idx) != length(colnames(expr))){
          warning("Some expression probes where not available in expression data")
         }
         if(length(idx) == 0){
           warning("No expression probes available, skipping.")
           return(NULL)
         }
         expr.probes <- expr.probes[idx]
         expr <- expr[,expr.probes,drop=F]
      }
      expr <- expr[,!grepl("NA", colnames(expr)),drop=F]
    }

    pheno <- phe
    rownames(pheno) <- pheno$Sample.ID
    pheno <- pheno[rownames(geno),,drop=F]
    # change some colnames to correpons to the same names used in the KORA data
    cnames <- colnames(pheno)
    colnames(pheno)[grepl("Sex",cnames)] <- "sex"
    colnames(pheno)[grepl("Age",cnames)] <- "age"
    colnames(pheno)[grepl("RNA_conv_batch",cnames)] <- "batch1"
    colnames(pheno)[grepl("RNA_extr_batch",cnames)] <- "batch2"
    pheno[,"batch1"] <- factor(pheno[,"batch1"])
    pheno[,"batch2"] <- factor(pheno[,"batch2"])

    data <- cbind.data.frame(pheno, stringsAsFactors=F)
    if(!is.na(expr.probes)) {
      data <- cbind.data.frame(data,expr,stringsAsFactors=F)
    }
    if(!is.na(meth.probes)){
      data <- cbind.data.frame(data, meth,stringsAsFactors=F)
    }
    if(!is.na(snp.ids)) {
      data <- cbind.data.frame(data, geno, geno.ids=rownames(geno),
                               stringsAsFactors=F)
    }
    return(data)

  } else {
    stop("Cohort not supported.")
  }
}

#' Method to get the genotypes for a certain set of SNPs in the specified
#' cohort
#'
#' @param snp.ranges A GRanges of SNPs for which to get genotypes. Needed for
#' KORA cohort only
#' @param snp.ids A vector of SNP-ids for which to get genotypes. Needed for
#' LOLIPOP cohort only
#' @param cohort A string specifiying the cohort for which to load the genotypes
#' either "kora" or "lolipop"
#'
#' @return Matrix of genotypes with SNPs in the columns and subjects in the rows
#'
get.genotypes <- function(snp.ranges, snp.ids, cohort){

  if("kora" %in% cohort){
    if(is.null(snp.ranges) | is.na(snp.ranges)) {
      return(NULL)
    }
    geno <- scan.snps(snp.ranges)
    if(nrow(geno$snps) != length(snp.ranges)){
      warning("Some SNPs were NA in genotype data.")
    }
    geno <- t(geno$snps)
    # get rid of non-changing snps
    geno <- geno[,apply(geno,2,var)!=0, drop=F]

  } else {
    if(!exists("dosage.oX_ia")){
      load("data/current/meQTLs/ggmdata_201017.RData")
    }
    geno <- t(dosage.oX_ia)
    snp.ids <- unique(snp.ids)
    idx <- which(snp.ids %in% colnames(geno))
    if(length(idx) != length(colnames(geno))){
      warning("Some SNPs where not available in geno type data")
    }
    snp.ids <- snp.ids[idx]
    # get rid of non-changing snps
    geno <- geno[,snp.ids,drop=F]
    geno <- geno[,apply(geno,2,var) != 0, drop=F]
  }

  return(data.matrix(geno))
}

#'
#' For a list of symbols, gets all annotated array ids from the illuminaHumanV3.db
#'
#' @param symbols List of symbols for which to get ids
#' @param mapping Flag whether to return a mapping ID->SYMBOL (default: F)
#' @param as.list Flag whether to return result as list or vector in case mapping=F. (default:F)
#'
#' @author Johann Hawe
#'
probes.from.symbols <- function(symbols, annot=NULL, mapping=F, as.list=F) {
  require(illuminaHumanv3.db)
  if(is.null(annot)) {
    annot <- as.list(illuminaHumanv3ALIAS2PROBE)
    annot <- annot[!is.na(annot)]
  }

  annot <- annot[which(names(annot) %in% symbols)]
  if(length(annot) > 0){
    if(mapping) {
      probes <- unlist(annot)
      uprobes <- unique(probes)
      if(length(uprobes) != length(probes)){
        dprobes <- unname(probes[duplicated(probes)])
        warning(paste0("Caution: Some probes had more than one gene annotated,
                       using only one symbol for those probes:\n",
                       dprobes))
      }
      map <- matrix(nrow=length(unlist(annot)))
      rownames(map) <- unlist(annot)
      temp <- lapply(symbols, function(s) {
        ids <- annot[[s]]
        if(!is.null(ids)) {
          map[ids,1] <<- s
        }
      })
      map <- map[!is.na(map[,1]),,drop=F]
      colnames(map) <- "symbol"
      dropped <- length(symbols) - nrow(map)
      if(dropped>0) {
        cat("Dropped",dropped,"symbols since they had no probe.id available.\n")
      }
      return(map)
    }
    if(as.list) {
      res <- annot[symbols]
      names(res) <- symbols
      return(res)
    } else {
      return(unlist(annot[symbols]))
    }
  } else {
    return(NULL)
  }
}

#'
#' Method to get gene level estimates from expression probes
#' (=rowMeans on respective probe.ids/gene)
#'
#' @param m The expression matrix with probe ids in columns
#' @param symbols The gene symbols to be used for the summarization. Only probe ids
#' belonging to a gene symbol in this list a used. Rest of probes is dropped.
#' @param threads Number of threasd which can be used
#'
#' @value Matrix of summarized probe values (per gene probe levels)
#'
#' @author Johann Hawe
#'
#' @date 02/06/2017
#'
#' @export
#'
summarize <- function(m, symbols, threads=1){
  library(illuminaHumanv3.db)
  annot <- as.list(illuminaHumanv3ALIAS2PROBE)
  annot <- annot[!is.na(annot)]
  n <- mclapply(symbols, function(sym) {

    gid <- probes.from.symbols(sym, annot=annot)

    # for LST1 symbol (and possibly others) there is a probe which should likely
    # be only appointed to the NFKBIL1 gene (according to UCSC genome browser).
    # Might be an error in annotation. Remove probe for LST1 manually for now
    if(sym=="LST1"){
      gid <- gid[which(!(gid %in% "ILMN_2046344"))]
    }
    if(is.null(gid)){
      return(NULL)
    }
    if(length(gid)>0){
      if(!all(gid %in% colnames(m))){
        gid <- gid[which(gid %in% colnames(m))]
        if(length(gid)<1){
          warning("None of the probeids for ",
                  sym,
                  " where available in expression data.")
          return(NULL)
        } else {
          warning("Some probes where not available in expression data.")
        }
      }

  	  g <- m[,gid,drop=F]
  	  s <- as.matrix(rowMeans(g), ncol=1)
  	  colnames(s) <- sym
      return(s)
    }
    return(NULL)
  }, mc.cores=threads)

  names(n) <- symbols
  n <- n[!sapply(n,is.null)]
  result <- matrix(unlist(n), ncol=length(n), byrow=F)
  colnames(result) <- names(n)
  rownames(result) <- rownames(m)
  return(result)
}

#' Uses the SNIPA interface to get LD snps for a specified snp. Then checks whether
#' any of the LD snps or their aliases are in the GWAS catalog and returns the
#' respective subset of the catalog.
#'
#' @param SNP RS-id of the SNP to check
#'
#' @return The subset of the Gwas catalog (data.frame) associated with the SNP
#'
#' @author Johann Hawe
#'
#' @date 2017/04/11
#'
#'
get.gwas.hits <- function(SNP) {
  source("R/snipe.R")

  library(data.table)

  cat("Preparing GWAS catalog\n")
  gwas.cat <- suppressWarnings(fread("data/current/gwas/gwas_catalog_v1.0.1-associations_e88_r2017-04-03.tsv",
                    sep="\t", header=T, data.table=F))

  ld.snps <- snipa.get.ld.by.snp(SNP)
  cat("Got", nrow(ld.snps), "rows in SNP-LD data frame\n")

  # get all related rs-ids for the SNP
  if(!all(is.na(ld.snps$RSALIAS))){
    aliases <- unlist(strsplit(ld.snps$RSALIAS, ","))
  } else {
    aliases <- ""
  }
  snps <- unique(trimws(c(ld.snps$QRSID, ld.snps$RSID, aliases)))
  # get gwas overlap with rsids
  result <- gwas.cat[which(gwas.cat$SNPS %in% snps), ,drop=F]

  return(result)
}

#' Gets the correlation pvalue matrix between all entities of the supplied
#' data matrix.
#'
#' @param data The data matrix with the individuals in the rows and the measurements
#' in the columns
#'
#' @param qvalue Flag whether to return qvalues instead of pvalues in the matrix.
#' Default: False
#'
#' @return A matrix of p-values or q-values between each possible pair
#' of data entities
#'
#' @author Johann Hawe
#'
get.correlation.pval.matrix <- function(data, qvalue=F){
  library(reshape)
  N <- ncol(data)
  result <- matrix(-1, ncol=N, nrow=N)
  rownames(result) <- colnames(result) <- colnames(data)
  diag(result) <- 1
  for(i in 1:(ncol(data)-1)){
    for(j in (i+1):ncol(data)){
      pval <- cor.test(data[,i], data[,j])$p.value
      result[i,j] <- pval
    }
  }
  if(qvalue){
    # replace pvalues in matrix with the qvalues
    melted <- melt(result)
    # remove self tests
    melted <- melted[which(melted$X1 != melted$X2),]
    melted <- melted[!melted$value==-1,]
    melted$q <- p.adjust(melted$value, "BH")
    for(i in 1:nrow(melted)){
      result[melted[i,"X1"], melted[i,"X2"]] <- melted[i,"q"]
      result[melted[i,"X2"], melted[i,"X1"]] <- melted[i,"q"]
    }
  } else {
    # copy upper tri matrix to lower tri (automatically done when using qvalues)
    for(i in 1:nrow(result)) {
      for(j in 1:(i-1)) {
        result[i,j] <- result[j,i]
      }
    }
  }
  return(result)
}

#' From a symmetric matrix of (correlation) pvalues of measurements,
#' creates a a graphNEL structure for a given cutoff.
#'
#' @param cor.pvals The matrix of correlation pvalues
#' @param cutoff The pvalue cutoff to be used
#'
#' @return A graphNEL object containing significant associations between the
#' entities
#'
#' @author Johann Hawe
#'
create.correlation.network <- function(cor.pvals, cutoff) {
  library(graph)

  if(!isSymmetric(cor.pvals)){
    stop("Matrix of correlation pvalues must be symmetric!")
  }
  edges <- matrix(ncol=3, nrow=0)
  colnames(edges) <- c("n1", "n2", "pval")

  cn <- colnames(cor.pvals)
  N <- ncol(cor.pvals)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      pv <- cor.pvals[i,j]
      if(pv < cutoff){
        edges <- rbind(edges, c(cn[i], cn[j], pv))
      }
    }
  }
  graph <- graphNEL(union(edges[,1], edges[,2]),edgemode = "undirected")
  graph <- addEdge(edges[,1],
                         edges[,2],
                         weights=1-as.numeric(edges[,3]),
                         graph)
  return(graph)
}

#' Gets genomic ranges for a list of cpg-identifiers using the
#' infiniummethylation.hg19 annotation
#'
#' @param cpgs IDs of cpgs for whcih to get ranges
#'
#' @author Johann Hawe
#'
get.cpg.ranges <- function(cpgs){
  library(FDb.InfiniumMethylation.hg19)
  ranges <- features(FDb.InfiniumMethylation.hg19)
  ranges <- ranges[cpgs]
  return(ranges)
}

#' Annotates a regulatory graph with appropriate node and
#' edge attributes
#'
#' @param g The graph to be annotated
#' @param ranges The ranges to be used for the annotation. Contains
#' information on which genes are TFs, what the CpG genes are etc.
#'
#' @return The same graph instance as given, but annotated with specific
#' node and edge attributes
#'
#' @author Johann Hawe
#'
annotate.graph <- function(g, ranges){
  library(graph)
  # work on all graph nodes
  gn <- nodes(g)

  nodeDataDefaults(g,"cpg") <- F
  nodeData(g, gn, "cpg") <- grepl("^cg", gn)

  nodeDataDefaults(g,"snp") <- F
  nodeData(g, gn, "snp") <- grepl("^rs", gn)

  nodeDataDefaults(g,"snp.gene") <- F
  nodeData(g, gn, "snp.gene") <- gn %in% ranges$snp.genes$SYMBOL

  nodeDataDefaults(g,"cpg.gene") <- F
  nodeData(g, gn, "cpg.gene") <- gn %in% ranges$cpg.genes$SYMBOL

  nodeDataDefaults(g, "tf") <- F
  nodeData(g, gn, "tf") <- gn %in% ranges$tfs$SYMBOL

  nodeDataDefaults(g, "sp.gene") <- F
  nodeData(g, gn, "sp.gene") <- gn %in% ranges$spath$SYMBOL

  # edge data information
  edgeDataDefaults(g, "isChipSeq") <- FALSE
  edgeDataDefaults(g, "isPPI") <- FALSE

  # we will also set tf2cpg priors at this point, so get all TFs
  tfs <- ranges$tfs$SYMBOL
  # also get the chipseq context for our cpgs
  context <- get.chipseq.context(names(ranges$cpgs))

  em <- matrix(ncol=2,nrow=0)
  # for all cpgs
  for(c in rownames(context)){
    for(tf in tfs) {
      # might be that the TF was measured in more than one cell line
      if(any(context[c,grepl(tf, colnames(context))])) {
        em <- rbind(em,c(c,tf))
      }
    }
  }
  if(nrow(em) > 0){
    em <- filter.edge.matrix(g,em)
    if(nrow(em)>0){
      graph::edgeData(g,em[,1], em[,2],"isChipSeq") <- T
    }
  }
  # ppi edgedata
  load.string.db()
  # get subset of edges which are in our current graph
  STRING.SUB <- subGraph(intersect(nodes(STRING.DB), gn), STRING.DB)
  edges <- t(edgeMatrix(STRING.SUB))
  edges <- cbind(gn[edges[,1]], gn[edges[,2]])
  edges <- filter.edge.matrix(g,edges)
  if(nrow(edges) > 0){
    edgeData(g,edges[,1], edges[,2],"isPPI") <- T
  }
  return(g)
}

#' Takes a gene expression matrix and adjusts the genes' expression for cis-eQTLs
#' using gtex eqtl results.
#'
#' @param g Gene expression matrix with all the genes in the column. Column names
#' need to be gene symbols.
#' @param eqtls List of eqtls for which to adjust the data
#'
#' @return The gene expression matrix adjusted for cis-eQTLs
#'
#' @author  Johann Hawe
#'
adjust.cis.eqtls <- function(g, eqtls, geno.data, cpg_genes, threads=1){

  # for all cpg-genes, adjust for potential cis-eqtls (i.e. snp-effects)
  toUse <- which(eqtls$gene %in% colnames(g) & eqtls$gene %in% cpg_genes)
  if(length(toUse) == 0){
    return(g)
  }

  # for each gene, check whether we have an associated cis-eqtl
  temp <- mclapply(colnames(g), function(x) {
    # indices of ciseqtls
    idxs <- which((eqtls$gene == x) & names(eqtls) %in% colnames(geno.data))
    if(length(idxs) < 1) {
      return(NULL)
    }
    if(x %in% cpg_genes) {
      snpids <- names(eqtls[idxs])
      X <- cbind.data.frame(g[,x,drop=F], geno.data[,snpids,drop=F])
      model <- lm(as.formula(paste0("`", x, "` ~",
                                    paste0(snpids,collapse="+"))),
                  data =  X,
                  na.action=na.exclude)

      r <- resid(model)
      return(r)
    } else {
      return(g[,x])
    }
  }, mc.cores=threads)

  names(temp) <- colnames(g)

  # set results to our data matrix
  cnt <- 0
  for(i in 1:length(temp)){
    x <- names(temp)[i]
    if(!is.null(temp[[i]])) {
      g[,x] <- temp[[i]]
      cnt <- cnt + 1
    }
  }

  cat("Adjusted expression value of", cnt, "genes.\n")
  return(g)
}

#' Takes a methylation matrix and adjusts the cpgs' expression for cis-meQTLs
#' using currently the BONDER meqtl results.
#'
#' @param c Methylation beta matrix with all the cpgs in the columns.
#'
#' @return The methylation matrix adjusted for cis-meQTLs
#'
#' @author  Johann Hawe
#'
adjust.cis.meqtls <- function(c, meqtls, geno.data, threads=1){
  # for all cpgs, adjust for potential cis-meqtls (i.e. snp-effects)
  toUse <- which(meqtls$cpg %in% colnames(c))
  if(length(toUse) == 0){
    return(c)
  }

  meqtls <- meqtls[toUse]

  # for each gene, check whether we have an associated cis-meqtl
  temp <- mclapply(colnames(c), function(x) {
    # indices of cis-meqtls
    idxs <- which((meqtls$cpg == x) & names(meqtls) %in% colnames(geno.data))
    if(length(idxs) < 1) {
      return(NULL)
    }
    snpids <- names(meqtls[idxs])
    X <- cbind.data.frame(c[,x,drop=F], geno.data[,snpids,drop=F])
    # drop data where all values are the same
    X <- X[,unlist(apply(X,2,function(x){ var(x,na.rm = T)!=0 })), drop=F]
    if(ncol(X)<1 | !x %in% colnames(X)){
      return(NULL)
    }
    model <- lm(as.formula(paste0("`", x, "` ~",
                                  paste0(snpids[snpids %in% colnames(X)],collapse="+"))),
                data =  X,
                na.action=na.exclude)

    r <- resid(model)
    return(r)
  }, mc.cores=threads)

  names(temp) <- colnames(c)

  # set results to our data matrix
  cnt <- 0
  for(i in 1:length(temp)){
    x <- names(temp)[i]
    if(!is.null(temp[[i]])) {
      c[,x] <- temp[[i]]
      cnt <- cnt + 1
    }
  }

  cat("Adjusted methylation beta of", cnt, "cpgs\n")
  return(c)
}

#' Get gene symbols for the given (array) probe ids
#'
#' @param ids List of probe ids for which to get the symbols
#' @param mapping Whether to return a mapping (probe->symbol); default: F
#'
#' @author Johann Hawe
#'
#' @return If mapping=F: All symbols referenced by the given probe ids. Otherwise
#' a mapping of all probeids to the corresponding symbols (a data.frame)
#'
symbols.from.probeids <- function(ids, mapping=F){
  library(illuminaHumanv3.db)
  annot <- illuminaHumanv3SYMBOLREANNOTATED
  avail <- mappedkeys(annot)
  ids.avail <- ids[which(ids %in% avail)]

  symbols <- unlist(as.list(annot[ids.avail]))

  if(mapping) {
    return(data.frame(id=ids.avail,symbol=symbols))
  } else {
    return(symbols)
  }
}

#' Load cis-eqtls identified in KORA
#'
#' Loads a list of cis-eqtls identified previously in the KORA F4 cohort.
#' (Supplement Table S2 of Schramm et al.).
#'
#' @param expr.probes Optional. Restrict the eqtls to be loaded to the this list of expression
#' probes obly
#'
#' @author Johann Hawe
#'
#' @return GRanges objects containing the eqtl information
#'
load.eqtls <- function(feqtls, expr.probes=NULL) {
  library(GenomicRanges)
  # columns: top SNP KORA F4;Chr SNP;minor allele KORA F4;MAF KORA F4;Probe_Id;Gene;Chr Gene; etc...
  eqtls <- read.csv(feqtls, sep=";", header=T, stringsAsFactors=F)
  # load snp ranges
  snps <- get_snpPos_biomart(eqtls$top.SNP.KORA.F4)
  snps <- with(snps, GRanges(paste0("chr", chr),
                             IRanges(start, width=1), id=snp))
  names(snps) <- snps$id

  sids <- eqtls$top.SNP.KORA.F4
  keep <- sids %in% names(snps)
  if(is.null(expr.probes)){
    keep <- keep & (eqtls$Probe_Id %in% expr.probes)
  }

  sids <- sids[keep]
  snps <- snps[sids]
  genes <- eqtls$Gene[keep]
  pids <- eqtls$Probe_Id[keep]

  # try to map some additional symbols
  temp <- lapply(1:length(genes), function(i) {
    if(genes[i] == "") {
      sym <- symbols.from.probeids(pids[i])
      genes[i] <<- sym
    }
  })

  # drop duplicates (since the analysis was probebased and we are only on the
  # per gene level)
  keep <- which(!duplicated(paste0(genes, names(snps))))
  eqtls <- snps[keep]
  eqtls$probe_id <- pids[keep]
  eqtls$gene <- genes[keep]
  return(eqtls)
}

#' Load cis meqtls identified in our study
#'
#' @param cpg.probes Optional. Restrict the list of meqtls on those cpg-probes only.
#'
#' @author Johann Hawe
#'
#' @return GRanges object reflecting all found cis-meQTLs
#'
load.meqtls <- function(fmeqtls, cpg.probes=NULL) {
  library(GenomicRanges)

  load(fmeqtls)

  if(is.null(cpg.probes)){
    cosmosub <- cosmo
  } else {
    cosmosub <- cosmo[(cosmo$cpg %in% cpg.probes),,drop=F]
  }
  rm(cosmo)

  if(nrow(cosmosub) > 0){
    # load snp ranges (columns: chr pos id)
    snps <- with(cosmosub, GRanges(paste0("chr", snp.chr),
                                   IRanges(snp.pos, width=1), id=snp, cpg=cpg))
    names(snps) <- snps$id

    return(unique(snps))
  } else {
    return(NULL)
  }
}

#' Handles the NAs of CpGs potentially present in the given
#' data matrix
#'
#' CpGs are removed from the data matrix, if they show NA in
#' more than 1% of the available data points. After the removal,
#' any obeservations still containing NA values for CpGs are removed
#' by using the \code{complete.cases} method
#'
#' @param data The data which to check for NA cpg values+
#'
#' @return The same data matrix, but where CpGs and/or individuals
#' showing to many NAs are removed
#'
#' @author Johann Hawe
#'
handle.na.cpgs <- function(data){
  # devide data (cpgs/non-cpgs)
  cgdata <- data[,grepl("^cg", colnames(data))]
  rdata <- data[,!grepl("^cg", colnames(data))]

  # remove cpgs which show NA in more than 1% of the samples
  N <- nrow(data)*0.01
  cgdata <- cgdata[,!unlist(apply(cgdata,2,function(x) sum(is.na(x))>N))]
  ndata <- cbind(rdata,cgdata)

  # drop incomplete cases
  cc <- complete.cases(ndata)
  ndata <- ndata[cc,,drop=F]

  return(ndata)
}

#' Creates an annotation object, mapping TFBS to TSS
#'
#' Loads all available TFBS collected from public sources (Encode, Remap) and
#' overlaps those with the provided TSS.
#' Code adapted from file R/annotate-cpgs.R.
#'
#'
#' @author Johann Hawe
#'
annotate_tfbs_to_tss <- function(fbinding_sites_remap,
                                 fbinding_sites_encode,
                                 tss) {
  library(rtracklayer)
  library(data.table)

  # get the TFBS regions from remap
  tfbs = import(fbinding_sites_remap)
  ann = t(matrix(unlist(strsplit(values(tfbs)[,"name"], ".", fixed=T)), nrow=3))
  colnames(ann) = c("geo_id", "TF", "condition")
  values(tfbs) = DataFrame(name=values(tfbs)[,"name"],
                           data.frame(ann, stringsAsFactors=F))

  # we write out a table with all conditions and select the blood related ones
  conditions = t(matrix(unlist(strsplit(unique(values(tfbs)[,"name"]), ".",
                                        fixed=T)), nrow=3))
  colnames(conditions) = c("geo_id", "TF", "condition")
  conditions = conditions[order(conditions[,"condition"]),]
  conditions = conditions[,c(1,3)]
  conditions = conditions[!duplicated(paste(conditions[,1], conditions[,2])),]
  conditions = data.frame(conditions, blood.related=F)
  for (term in c("amlpz12_leukemic", "aplpz74_leukemia",
                 "bcell", "bjab", "bl41",
                 "blood", "lcl", "erythroid", "gm",
                 "hbp", "k562", "kasumi",
                 "lymphoblastoid", "mm1s", "p493",
                 "plasma", "sem", "thp1", "u937")) {
    conditions[grep(term, conditions[,2]),"blood.related"] = TRUE
  }

  # select the appropriate blood related TFBS subset
  selected = tfbs[values(tfbs)[,"condition"] %in%
                    conditions[conditions[,"blood.related"],"condition"]]

  # load the encode tfs separately
  encode = as.data.frame(fread(fbinding_sites_encode, header=F))
  encode = with(encode, GRanges(seqnames=V1, ranges=IRanges(V2 + 1, V3),
                                name=paste("ENCODE", V4, tolower(V6), sep="."),
                                geo_id="ENCODE", TF=V4,
                                condition=tolower(V6)))

  # filter blood related cell lines
  encode.lcl = encode[grep("gm", values(encode)[,"condition"])]
  values(encode.lcl)[,"condition"] = "lcl"
  encode.k562 = encode[grep("k562", values(encode)[,"condition"])]
  values(encode.k562)[,"condition"] = "k562"

  # combine remap and encode TFBS
  selected = c(selected, encode.lcl, encode.k562)

  # create an annotation matrix for the TSS
  chip = paste(values(selected)[,"TF"], values(selected)[,"condition"], sep=".")
  chip_exp = unique(chip)

  tfbs_ann = sapply(chip_exp, function(x) overlapsAny(tss,
                                                      selected[chip == x]))
  rownames(tfbs_ann) = names(tss)

  return(tfbs_ann)
}

#' Method to query biomart for snp information by genomic region.
#' Currently only gets the rsIDs.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
get_snps_biomart <- function(ranges, grch=37, additional.info=F) {

  require(biomaRt)

  # create filter values
  chrs <- gsub("chr", "", as.character(seqnames(ranges)))
  start <- start(ranges)
  end <- end(ranges)

  # set filter names
  filt <- c("chr_name", "start", "end")

  # set attributes to get
  attr <- c("refsnp_id", "chr_name", "chrom_start", "chrom_end")
  if(additional.info) {
    attr <- c(attr, "consequence_type_tv")
  }

  # query biomart
  mart <- useEnsembl(biomart="snp", GRCh=grch, dataset="hsapiens_snp")
  snps <- getBM(attr, filt,
                values=list(chrs, start, end),
                mart=mart,
                uniqueRows=T)

  # set new colnames
  if(additional.info) {
    colnames(snps) <- c("snp", "chr", "start", "end", "variant_consequence")
  } else {
    colnames(snps) <- c("snp", "chr", "start", "end")
  }

  return(snps)
}

#' Method to query biomart for snp position information by their rsIds.
#'
#' @param rsIds Vector of rsIds for which to get the genomic position info
#' @param grch The genome assembly version. Default: 37
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
get_snpPos_biomart <- function(rsIds, grch=37, additional.info=F) {

  require(biomaRt)

  # create filter values
  vals <- unique(as.character(rsIds))

  # set filter names
  filt <- c("snp_filter")

  # set attributes to get
  attr <- c("refsnp_id", "chr_name", "chrom_start", "chrom_end")

  # check whether to get some additional variant information
  if(additional.info) {
    attr <- c(attr,  "consequence_type_tv")
  }
  # query biomart
  mart <- useEnsembl(biomart="snp", GRCh=grch, dataset="hsapiens_snp")
  snps <- getBM(attr, filt,
                values=list(vals),
                mart=mart,
                uniqueRows=T)
  # set new colnames
  if(additional.info) {
    colnames(snps) <- c("snp", "chr", "start", "end", "variant_consequence")
  } else {
    colnames(snps) <- c("snp", "chr", "start", "end")
  }

  return(snps)
}

#' Method to query biomart for snp+flanking region sequences.
#'
#' @param rsIds Vector of rsIds for which to get the genomic position info
#' @param downstr 20. Size of downstream flanking region to get sequence for
#' @param upstr 20. Size of upstream flanking region to get sequence for
#' @param grch The genome assembly version. Default: 37
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
get_snp_flanking_biomart <- function(rsIds, downstr=20, upstr=20, grch=37) {
  require(biomaRt)

  # create filter values
  vals <- list(snp_filter=unique(as.character(rsIds)),
               upstream_flank=downstr,
               downstream_flank=upstr)

  # set filter names
  filt <- c("snp_filter", "upstream_flank", "downstream_flank")

  # set attributes to get
  attr <- c("refsnp_id", "snp")

  # query biomart
  mart <- useEnsembl(biomart="snp",
                     GRCh=grch,
                     dataset="hsapiens_snp")
  snps <- getBM(attr, filt, values=vals, mart=mart, uniqueRows=T, checkFilters=F)

  return(snps)
}

# ------------------------------------------------------------------------------
#' Helper function to get distance matched cpgs for a specific SNP
#'
#' @param snp The SNP for which to get the distance matched CpGs
#' @param bg_cpgs All background cpgs
#' @param dl Lower distance bound
#' @param du Upper distance bound
#' @param matchd_cpgs The cpgs matched according to methylation levels
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------
get_distance_matched_cpgs <- function(snp, bg_cpgs,
                                      dl, du, matched_cpgs, is.trans=F) {
  if(is.null(bg_cpgs$chr)) {
    stop("Expecting extra chromosome information as character for CpGs!")
  }
  
  if(is.trans) {
    # ensure different chromosomes for trans analysis
    bg_cpgs_sub <- bg_cpgs[which(bg_cpgs$chr != as.character(seqnames(snp)))]

    # only allow beta-matched cpgs
    use <- which(names(bg_cpgs_sub) %in% matched_cpgs)
    scpgs <- bg_cpgs_sub[use]
  } else {
    # same chromosome for longrange
    bg_cpgs_sub <- bg_cpgs[which(bg_cpgs$chr == as.character(seqnames(snp)))]
    dists <- distance(snp, bg_cpgs_sub)

    # distance filter and only allow beta-matched cpgs
    use <- which(names(bg_cpgs_sub) %in% matched_cpgs &
                   dl <= dists &
                   dists <= du)
    scpgs <- bg_cpgs_sub[use]
  }
  return(scpgs)
}

# ------------------------------------------------------------------------------
# helper to get leveled contingency table for meQTL enrichment analysis
# ------------------------------------------------------------------------------
get_cont_table <- function(meqtl, bg) {
  nmeqtl <- rep("meQTL", length(meqtl))
  m <- ifelse(meqtl, "overlap", "no_overlap")

  nbackground <- rep("background", length(bg))
  bg <- ifelse(bg, "overlap", "no_overlap")

  df <- cbind.data.frame(type=c(nmeqtl, nbackground),
                         affiliation=c(m, bg),
                         stringsAsFactors=F)
  df$type <- factor(df$type,
                    levels=c("meQTL", "background"),
                    ordered=T)
  df$affiliation <- factor(df$affiliation,
                           levels=c("overlap", "no_overlap"),
                           ordered=T)
  return(table(df$type, df$affiliation))
}

# ------------------------------------------------------------------------------
#' Method to quickly set a vector of default plotting colors
#' to be used in all generated plots for consistency
#'
#' TODO: parameter for number of colors to get. Then probably
#' the palettes have to be switched accordingly
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------
set_defaultcolors <- function() {
  library(wesanderson)
  cols <- wes_palette(n=5, "FantasticFox1")
  return(cols)
}

# ------------------------------------------------------------------------------
#' Small helper for a list() function automatically setting
#' the variable names as the names of the created list
#'
#' compare: https://stackoverflow.com/a/21059868
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
listN <- function(...){
  ln <- list(...)
  names(ln) <- as.character(substitute(list(...)))[-1]
  ln
}

# ------------------------------------------------------------------------------
#' Obtain houseman weighted chromHMM annotations.
#' Houseman info should be matched (e.g. when using EPIC cpgs -> use EPIC sample
#' houseman estimates)
#' 
#' @param ranges Ranges which to annotate. Should contain a meta column 'id'
#' to set rownames of result properly.
# ------------------------------------------------------------------------------
get_weighted_chromhmm_annotation <- function(ranges,
                                             fhouseman = "data/current/kora/methylation/Houseman/KF4_QN_estimated_cell_distribution_meanimpute473_lessThanOneTRUE.csv",
                                             froadmap_samples = "data/current/roadmap/sample_info.txt",
                                             fmnemonics = "data/current/roadmap/chromHMM/15state/mnemonics.txt") {
  # load roadmap information
  roadmap.samples = read.csv(froadmap_samples,
                             sep="\t", stringsAsFactors=F)
  
  mnemonics = read.csv(fmnemonics,
                       sep="\t")
  levels = with(mnemonics, paste(STATE.NO., MNEMONIC, sep="_"))
  labels = mnemonics[,"DESCRIPTION"]
  use = roadmap.samples$ANATOMY == "BLOOD"
  ids = roadmap.samples[use,"Epigenome.ID..EID."]
  cells = roadmap.samples[use,"Standardized.Epigenome.name"]
  names(cells) = ids
  
  # match this to the houseman cell types
  type = rep("other", sum(use))
  type[grep("CD4", roadmap.samples[use,5])] = "CD4T"
  type[grep("CD8", roadmap.samples[use,5])] = "CD8T"
  type[grep("CD15", roadmap.samples[use,5])] = "Gran"
  type[grep("CD56", roadmap.samples[use,5])] = "NK"
  type[grep("CD19", roadmap.samples[use,5])] = "Bcell"
  type[grep("CD14", roadmap.samples[use,5])] = "Mono"
  
  # check the average cell composition in our samples from the houseman estimates
  houseman = read.csv(fhouseman, sep=";")
  pop.mean = colMeans(houseman[,-1])
  
  # load the annotation from chromHMM
  ann <- get_chromHMM_annotation(ranges, ids)
  
  # weight by houseman population means
  weighted <- t(apply(ann, 1, function(ann) {
    ## make a cell x state indicator matrix
    tab = sapply(levels, function(l) as.numeric(ann == l))
    ## average by cell type
    by.type = apply(tab, 2, tapply, type, mean)
    ## weight by population means
    weighted = pop.mean %*% by.type[names(pop.mean),]
    return(weighted)
  }))
 
  colnames(weighted) = labels
  if(!is.null(ranges$name)) {
    rownames(weighted) <- rownames(ann) <- ranges$name
  } else if (!is.null(ranges$id)) {
    rownames(weighted) <- rownames(ann) <- ranges$id  
  } else {
    rownames(weighted) <- rownames(ann) <- names(ranges)
  }
  
  return(weighted)
}
