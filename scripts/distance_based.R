library("mccbn")
library("optparse")
parser <- OptionParser(usage="usage: %prog [options] FILES")
parser <-
    add_option(parser, c('-B', '--num_samples'), type="integer", default=100,
               metavar='INT', help="Number of samples [default %default]",
               dest='num_samples')
parser <-
    add_option(parser, c('-l', '--log'), type="character", default="split.log",
               metavar='LOGFILE', help='Log file [default %default]',
               dest='log')
parser <-
    add_option(parser, c('-o', '--output'), type="character", default="M0.rds",
               metavar='OUTFILE', help='Output file [default %default]',
               dest='outfile')
args <- parse_args(parser, positional_arguments=c(1, Inf))
FILES <- args$args
args <- args$options

log <- file(args$log, open="wt")
sink(log)
sink(log, type="message")


jaccard.distance <- function(p1, p2) {
  aux <- p1 + p2
  distance <- 1 - (sum(aux == 2) / sum(aux > 0))
  return(distance)
}

aggregate <- function(samples, D) {
    df <- data.frame()
    posets <- list()  
    for (s in samples) {
      infile <- paste("asa_D", D, "_s", s - 1, ".rds", sep="")
      infile_table <- paste("D", D, "_s", s - 1, "_params.txt", sep="")
      
      if (file.exists(file.path("output", infile))) {
            ret <- readRDS(file.path("output", infile))
      } else {
        cat("skipping sample", s - 1, "\n")
        df[s, 1] <- s - 1
        posets[[s]] <- matrix(0, p, p)
        next
      }
      infile_table <- file.path("output", infile_table)
      df[s, 1] <- s - 1
      df[s, 2] <- ret$eps
      p <- length(ret$lambda)
      j <- 3 + p - 1
      for (i in 1:p) {
        df[s, 3 + i - 1] <- ret$lambda[i]
      }
      df[s, j + 1] <- ret$llhood
      if (file.exists(infile_table)) {
        params <- read.table(infile_table, sep="\t", skip=1)
        df[s, j + 2] <- which.max(params[, 2])
        df[s, j + 3] <- max(params[, 2])
      } else {
        df[s, j + 2] <- NA
        df[s, j + 3] <- NA
      }
      df[s, j + 4] <- sum(ret$poset)
      posets[[s]] <- ret$poset
    }

    colnames(df) <- c("sample", "epsilon", paste("lambda", seq(1, p), sep=""),
                      "obs_llhood_avg_MCEM", "max_llhood_step", "max_llhood",
                      "num_edges")
    return(list(df=df, posets=posets))
}

cat("Aggregating results\n")
num.samples <- args$num_samples
D1 <- aggregate(1:num.samples, "1")
D2 <- aggregate(1:num.samples, "2")

jaccard.dist <- numeric(length=num.samples)
jaccard.dist.closure <- numeric(length=num.samples)

cat("Computing distances among reconstructed posets\n")
for (l in 1:num.samples) {
    poset1 <- D1$posets[[l]]
    poset1.closure <- trans_closure(poset1)
    poset2 <- D2$posets[[l]]
    poset2.closure <- trans_closure(poset2)

    jaccard.dist[l] <- jaccard.distance(poset1, poset2)
    jaccard.dist.closure[l] <- jaccard.distance(poset1.closure, poset2.closure)
}

posetD1 <- readRDS(FILES[2])
posetD2 <- readRDS(FILES[1])
distML <- jaccard.distance(trans_closure(posetD1$poset), trans_closure(posetD2$poset))
pvalue <- sum(jaccard.dist.closure >= distML) / num.samples

saveRDS(list("jaccard.dist"=jaccard.dist, "jaccard.dist.closure"=jaccard.dist.closure,
             "distML"=distML, "pvalue"=pvalue),
        file=args$outfile)

