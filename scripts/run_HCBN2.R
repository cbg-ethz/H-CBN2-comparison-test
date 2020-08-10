library("mccbn")
library("optparse")
parser <- OptionParser(usage="usage: %prog [options] FILES")
parser <-
    add_option(parser, c('-L', '--num_samples'), type="integer", default=1000,
               metavar='INT', help="Number of samples [default %default]",
               dest='num_samples')
parser <-
    add_option(parser, c('-s', '--sampling_mode'), type="character",
               default='forward', metavar='MODE', help="Sampling mode [default %default]",
               dest='sampling_mode')
parser <-
    add_option(parser, c('-i', '--iter_MCEM'), type="integer", default=100,
               metavar='INT', help="Number of iterations MCEM [default %default]",
               dest='iter_MCEM')
parser <-
    add_option(parser, c('-t', '--temperature'), type="integer", default=50,
               metavar='INT', help="Initial temperature [default %default]",
               dest='T0')
parser <-
    add_option(parser, c('-a', '--adaptive_rate'), type="numeric", default=0.1,
               metavar='FLOAT', help="Adaptive rate [default %default]",
               dest='adap_rate')
parser <-
    add_option(parser, c('--iter_ASA'), type="integer", default=17000,
               metavar='INT', help="Number of iterations ASA [default %default]",
               dest='iter_ASA')
parser <-
    add_option(parser, '--asa', action="store_true", default=FALSE, 
               help="Indicate whether to run ASA [default %default]",
               dest='adaptive')
parser <-
    add_option(parser, '--seed', type="integer", default=47,
               metavar='INT', help="Seed [default %default]",
               dest='seed')
parser <-
    add_option(parser, '--thrds', type="integer", default=4,
               metavar='INT', help="Threads [default %default]",
               dest='threads')
parser <-
    add_option(parser, '--out_prefix', type="character", default="D1_s0",
               metavar='CHAR', help="Ouput prefix [default %default]",
               dest='out_prefix')
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


observations <- readRDS(FILES[1])
p <- ncol(observations)
N <- nrow(observations)

# Initial poset
posets <- candidate_posets(observations, obs_weights=rep(1, N),
                           min_compatible_geno_fraction=0.8)
poset0 <- tail(posets, n=1)[[1]]
aux <- sum(apply(observations, 1, is_compatible, poset0))
cat("Number of edges (initial poset): ", sum(poset0), "\n")
cat("Number of compatible observations with poset: ", aux, "(",
    round(aux / N, digits=3),  ")\n")

update.step.size <- ceiling(0.8 * args$iter_MCEM)
ret <-
    adaptive.simulated.annealing(
        poset0, observations, times=NULL, lambda.s=1.0, L=args$num_samples,
        sampling=args$sampling_mode, max.iter=args$iter_MCEM,
        update.step.size=update.step.size, tol=0.001, T0=args$T0,
        adap.rate=args$adap_rate, acceptance.rate=1 / p, step.size=100,
        max.iter.asa=args$iter_ASA, adaptive=args$adaptive,
        outdir=args$out_prefix, thrds=args$threads, verbose=TRUE,
        seed=args$seed)

saveRDS(ret, file=args$outfile)
cat("Done!")
