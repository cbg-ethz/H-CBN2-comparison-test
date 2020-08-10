# read in command line arguments
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
args <- parse_args(parser, positional_arguments=c(1, 2))
FILES <- args$args
args <- args$options

log <- file(args$log, open="wt")
sink(log)
sink(log, type="message")

# set seed for reproducibility
set.seed(47, kind="default")

splits <- function(B, obs_events, N1, N2) {

  # number of observations / genotypes
  n <- nrow(obs_events)
  # number of mutations
  p <- ncol(obs_events)

  for (i in 0:(B - 1)) {
    
    cat("split #", i, "\n")

    idxs <- sample(n, N1, replace=FALSE)
    s1 <- obs_events[idxs, ]
    s2 <- obs_events[-idxs, ]
    stopifnot(nrow(s2) == N2)

    saveRDS(s1, file=file.path("splits", paste("D1_s", i, ".rds", sep="")))
    saveRDS(s2, file=file.path("splits", paste("D2_s", i, ".rds", sep="")))

  }
}

D1 <- readRDS(FILES[1])
D2 <- readRDS(FILES[2])

observations <- rbind(D1, D2)
saveRDS(observations, file=args$outfile)

splits(B=args$num_samples, obs_events=observations, N1=nrow(D1),
       N2=nrow(D2))
