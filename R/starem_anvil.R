## =====================================================================
## SRM 620 -- STAR-EM (with STAR-DR): ANVIL (Purdue RCAC) version
## Selorm Buaka
##
## Same methods as the verified desktop script (SRM620_STAREM_desktop.R); only the
## orchestration is different: the work is split into independent SLURM jobs.
##   STAR_MODE = prep     panel, descriptive tables, choice of K, true policy values,
##                        learned missingness model  -> prep.rds   (run once)
##   STAR_MODE = real     HeartSteps analysis on the fine delta grid (run once)
##   STAR_MODE = sim      one SLURM array task = one block of replications of one
##                        scenario -> sim_XXX.rds  (STAR_PILOT=1: timing pilot)
##   STAR_MODE = combine  all tables and figures from the finished tasks
## Every replication has its own seed (SEED + 100000 * scenario + replication), so
## tasks can run in any order, be resubmitted after a failure or an Anvil outage,
## and a later batch of replications never repeats an earlier one.
## =====================================================================

## ---------------- settings (identical to the verified desktop run unless noted)
NREP      <- 5       # simulation replications per condition (desktop check only; Anvil: 200-500)
SEED      <- 100     # master seed (set.seed(100) convention); STAR-EM starts/refits offset from it
REAL_A    <- FALSE   # FALSE: re-run the trial's own randomisation (P(send)=0.6 when available) in every replication.
                     # TRUE replays the 37 real send sequences, which freezes ONE realised randomisation
                     # and made even the oracle (complete-data) estimator biased by about -126 steps.
N_TRUTH   <- 1000    # resampled participants per policy for the true values (Anvil: 5000)
NBOOT     <- 499     # bootstrap draws per interval (Anvil: 499-1999)
N_TREES   <- 100     # trees per random forest (Anvil: 200)
N_THREADS <- max(1, parallel::detectCores() - 1)
DELTAS_REAL  <- c(-1, -0.5, 0, 0.5, 1)  # delta0 grid for the real data (Anvil: -1.5 to 1.5 by 0.25)
DELTAA_REAL  <- c(0, -0.5)               # extra tilt on days a suggestion was sent (Anvil: -1, -0.5, 0, 0.5)
M_IMP        <- 20                       # multiple imputations per missing reward (kNN-MI comparator)
K_IMP        <- 20                       # nearest observed donors for each missing reward
EM_BINS      <- 10                       # EM: reward histogram bins (1 zero bin + 9 quantile bins of positive steps)
EM_TOL       <- 1e-9                     # EM: stop when the relative change in log-likelihood is below this
EM_MAXIT     <- 2000                     # EM: maximum iterations
STAREM_BINS  <- 10                       # STAR-EM: reward bins (1 zero bin + quantile bins of positive steps)
STAREM_TREES <- 50                       # STAR-EM: trees in the M-step forest (Anvil: 150)
STAREM_MAXIT <- 20                       # STAR-EM: max EM iterations for the main fit (Anvil: 40)
STAREM_MINIT <- 5                        # STAR-EM: minimum EM iterations before the stopping rule applies
STAREM_TOL   <- 0.01                     # STAR-EM: stop when the expected mean of the missing rewards moves < 1% ...
STAREM_PTOL  <- 0.02                     # ... AND the state/missingness parameters move < 0.02, twice in a row
                                         # (the forest M-step makes the log-likelihood jitter, so it is not the stop signal)
STAREM_BROUNDS <- 30                     # STAR-EM missingness model: boosted trees refitted at every M-step
STAREM_BLR   <- 0.1                      # ... learning rate. 30 trees x 0.1 chosen by 5-fold cross-validation on
                                         #     HeartSteps (held-out log-lik -1512 vs -1549 for logistic + splines;
                                         #     flat between 30 and 40 trees; more trees overfit)
STAREM_BDEPTH <- 3                       # ... depth of each tree (up to 3-way interactions)
STAREM_RUNCAP <- 10                      # spell features are capped at 10 decisions (2 days)
STAREM_AVG   <- 3                        # STAR-EM: imputation model = average of the last 3 iterates (stochastic-EM practice)
STAREM_M     <- 5                        # STAR-EM: proper imputations = bootstrap EM refits (Anvil: 20)
STAREM_WARM  <- 3                        # STAR-EM: EM iterations per bootstrap refit (warm-started)
STAREM_KS    <- 1:3                      # STAR-EM: numbers of hidden day types compared (K = 1: no mixture layer)
STAREM_K     <- 2                        # STAR-EM: default number of hidden states; replaced by the K chosen on the real data
STAREM_STARTS <- 2                       # STAR-EM: random starts per fit (local maxima check; Anvil: 5)
## hidden-activity-state world (simulation scenario 7): a latent low/high activity state
## drives BOTH the steps and the missingness (shared-parameter MNAR); a notification can
## lift a person out of a low-activity spell. No distribution is assumed: low-state steps
## are the smaller of two resampled real step counts, high-state steps the larger.
HID_STAY_DAY <- 0.7                      # P(next day starts in the same state as this day started)
HID_UP0      <- 0.05                     # P(low -> high) at the next decision, no notification
HID_UPA      <- 0.30                     # extra P(low -> high) right after a notification
HID_DOWN     <- 0.05                     # P(high -> low) at the next decision
HID_GAM      <- 2                        # log-odds of missingness are higher by this in the low state
BOUND_Q      <- 0.99                     # upper reward bound for the worst-case bounds = 99th pct of observed steps

## ---------------- Anvil overrides (environment variables set in the slurm scripts)
envn <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }
MODE      <- Sys.getenv("STAR_MODE", "sim")
PILOT     <- Sys.getenv("STAR_PILOT") == "1"
TASK      <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
N_THREADS <- as.integer(envn("SLURM_CPUS_PER_TASK", 1))
REPS      <- envn("STAR_REPS", 1000)        # replications per scenario (journal version)
CHUNK     <- envn("STAR_CHUNK", 25)         # replications per array task, N = 37 scenarios
CHUNK_BIG <- envn("STAR_CHUNK_BIG", 2)      # replications per array task, N = 150 scenario
REPS_BIG  <- envn("STAR_REPS_BIG", 300)     # replications of the N = 150 scenario (a scaling check;
                                            # it is about 10x as costly per replication -- set from the pilot)
N_BIG     <- envn("STAR_NBIG", 150)         # participants in the large hidden-state scenario
N_TRUTH   <- envn("STAR_NTRUTH", 5000)      # resampled participants per policy for the true values
STAREM_STARTS <- envn("STAR_STARTS", 5)     # EM starts when K > 1 (local maxima); K = 1 uses one
NBOOT     <- envn("STAR_NBOOT", 499)
DELTAS_REAL <- seq(-1.5, 1.5, by = 0.25)    # finer real-data grid for the paper
DELTAA_REAL <- c(0, -0.5)
NREP <- REPS

## ---------------- output folder, log, packages (installed once by slurm/00_setup.sh)
out_dir <- path.expand(Sys.getenv("STAR_OUT", "~/starem_out"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, sprintf("log_%s%s.txt", MODE, if (MODE == "sim") sprintf("_%03d", TASK) else ""))
say <- function(...) { msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", ...)
                       cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE) }
cat("", file = log_file)
T_START <- Sys.time()
suppressPackageStartupMessages({ library(data.table); library(ranger); library(FNN); library(ggplot2) })
setDTthreads(N_THREADS)
say("mode ", MODE, if (MODE == "sim") paste0(" task ", TASK) else "", "; ", N_THREADS, " cores; output: ", out_dir)
figure <- function(name, p, w = 9, h = 4.2) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  say("saved figure ", name, ".png")
}
data_file <- path.expand(Sys.getenv("STAR_DATA", file.path(out_dir, "suggestions.csv")))
if (!file.exists(data_file)) stop("HeartSteps data not found at ", data_file, " -- run slurm/00_setup.sh on a login node first")
## =====================================================================
## 1. Decision-level panel
##    episode = participant-day (up to 5 decisions); A = suggestion sent;
##    R = Jawbone steps in the 30 min after the decision; O = R observed.
## =====================================================================
prepare_panel <- function(file) {
  raw <- fread(file)
  lg  <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1")
  d <- raw[!is.na(decision.index.nogap),
           .(id = as.integer(user.index), k = as.integer(decision.index.nogap),
             avail = lg(avail), send = lg(send), loc = as.character(dec.location.category),
             temp = suppressWarnings(as.numeric(dec.temperature)),
             precip = suppressWarnings(as.numeric(dec.precipitation.chance)),
             R = suppressWarnings(as.numeric(jbsteps30)))]
  d[, day := k %/% 5L + 1L]; d[, t := k %% 5L + 1L]
  d[, A := as.integer(avail & send)]
  d[, loc_home := as.integer(!is.na(loc) & loc == "home")]
  d[, loc_work := as.integer(!is.na(loc) & loc == "work")]
  d[temp < -100, temp := NA]; d[precip < 0, precip := NA]
  if (max(d$precip, na.rm = TRUE) > 1) d[, precip := precip / 100]
  d[is.na(temp), temp := median(d$temp, na.rm = TRUE)]
  d[is.na(precip), precip := median(d$precip, na.rm = TRUE)]
  d[, O := as.integer(!is.na(R))]
  qc <- c(mean = mean(log1p(d$R[d$O == 1])), sd = sd(log1p(d$R[d$O == 1])))
  setorder(d, id, day, t)
  d[, ep := .GRP, by = .(id, day)]
  d[, sent_today := cumsum(A) - A, by = ep]
  d[, lagobs := shift(O, fill = 0L), by = ep]
  d[, lagq := shift(fifelse(O == 1L, (log1p(R) - qc["mean"]) / qc["sd"], 0), fill = 0), by = ep]
  d[, prevobs := shift(O, fill = 1L), by = id]     # previous decision observed (across days)
  d[, prevmiss := 1L - prevobs]
  d[, prevrun := { m <- 1L - O; r <- integer(.N); cur <- 0L        # length of the missing run so far
                   for (j in seq_len(.N)) { r[j] <- min(cur, 10L); cur <- if (m[j] == 1L) cur + 1L else 0L }; r },
    by = id]
  ## each participant's real missing rate: leave-one-out for fitting the missingness
  ## model (no leakage), full rate carried by the participant into the simulation
  d[, `:=`(pm_full = mean(1 - O), pm_loo = (sum(1 - O) - (1 - O)) / pmax(.N - 1, 1)), by = id]
  panel <- d[, .(id, ep, day, t, avail, A, R, O, loc_home, loc_work, temp, precip,
                 sent_today, lagq, lagobs, prevobs, prevmiss, prevrun, pm_full, pm_loo)]
  attr(panel, "qc") <- qc
  panel
}

## =====================================================================
## 2. STAR-DR and comparison estimators
## =====================================================================
## 02_functions.R
## STAR-DR: Spell-aware, Tilted, Augmented-Reward Doubly Robust off-policy evaluation
## plus comparison estimators and distribution-free inference.
## No distribution is assumed for rewards anywhere in this file.


FEATS  <- c("t", "loc_home", "loc_work", "temp", "precip", "sent_today", "lagq", "lagobs", "prevobs", "day")
NTREE  <- N_TREES
NTHR   <- N_THREADS
GAMMA  <- 1                                   # undiscounted value within a day
P_SEND <- 0.6                                 # HeartSteps randomisation probability

q_fun <- function(R, qc) (log1p(R) - qc[1]) / qc[2]

## ---------------------------------------------------------------- policies
## Each returns P(A = 1 | history). Policies use only information that is always
## observed (availability, context, actions already taken), never the rewards.
POLICIES <- list(
  always   = function(d) as.numeric(d$avail),
  never    = function(d) rep(0, nrow(d)),
  tailored = function(d) as.numeric(d$avail & d$precip < 0.5 & d$t %in% 2:4),
  capped   = function(d) as.numeric(d$avail & d$sent_today < 2)
)
pib1 <- function(d) ifelse(d$avail, P_SEND, 0)         # known behaviour policy
CONTRASTS <- list(c("always", "never"), c("tailored", "always"), c("capped", "always"))

## ---------------------------------------------------------------- random forests
rf_fit <- function(X, y, classif = FALSE) {
  if (classif) {
    y <- factor(y, levels = c(0, 1))
    if (length(unique(y)) < 2) return(list(const = mean(as.numeric(as.character(y)))))
    ranger(x = X, y = y, num.trees = NTREE, min.node.size = 10, probability = TRUE,
           num.threads = NTHR)
  } else {
    if (length(unique(y)) < 2) return(list(const = y[1]))
    ranger(x = X, y = y, num.trees = NTREE, min.node.size = 10, num.threads = NTHR)
  }
}
rf_pred <- function(f, X, classif = FALSE) {
  if (!is.null(f$const)) return(rep(f$const, nrow(X)))
  p <- predict(f, data = X, num.threads = NTHR)$predictions
  if (classif) p[, "1"] else p
}
Xmat <- function(d, a = NULL) {
  X <- as.data.frame(d[, FEATS, with = FALSE])
  X$A <- if (is.null(a)) d$A else rep(a, nrow(d))
  X
}

## participant-level folds (cross-fitting never splits a person)
make_folds <- function(id, K = 5) {
  u <- unique(id); f <- sample(rep_len(1:K, length(u)))
  f[match(id, u)]
}

## ---------------------------------------------------------------- Step 1-3: nuisances
## p_obs = P(O = 1 | H, A); mu1 = E[R | H, A, O = 1];
## with delta(A) = delta0 + deltaA * A (row j of DS):
## e_j = E[exp(delta(A) q(R)) | H, A, O = 1]; nu_j = E[R exp(delta(A) q(R)) | H, A, O = 1]
## (A is a feature of every model, so conditioning on A is exact)
nuisances <- function(d, DS, qc, fold) {
  n <- nrow(d); K <- max(fold); J <- nrow(DS)
  p <- mu1 <- numeric(n)
  e  <- matrix(1, n, J); nu <- matrix(0, n, J)
  qR <- ifelse(d$O == 1, q_fun(d$R, qc), 0)
  for (k in 1:K) {
    tr <- fold != k; te <- fold == k
    X_tr <- Xmat(d[tr]); X_te <- Xmat(d[te])
    p[te] <- rf_pred(rf_fit(X_tr, d$O[tr], TRUE), X_te, TRUE)
    ob <- tr & d$O == 1
    X_ob <- Xmat(d[ob])
    mu1[te] <- rf_pred(rf_fit(X_ob, d$R[ob]), X_te)
    for (j in seq_len(J)) {
      if (DS$d0[j] == 0 && DS$dA[j] == 0) { nu[te, j] <- mu1[te]; next }
      w <- exp((DS$d0[j] + DS$dA[j] * d$A[ob]) * qR[ob])
      e[te, j]  <- pmax(rf_pred(rf_fit(X_ob, w), X_te), 1e-6)
      nu[te, j] <- rf_pred(rf_fit(X_ob, d$R[ob] * w), X_te)
    }
  }
  list(p = pmin(pmax(p, 0.02), 1), mu1 = mu1, e = e, nu = nu, qR = qR)
}

## Step 4: tilted augmented pseudo-reward (defined for every decision)
pseudo_reward <- function(d, nz, j, d0, dA = 0) {
  delta <- d0 + dA * d$A                              # decision-specific tilt
  p <- nz$p; mu0 <- nz$nu[, j] / nz$e[, j]
  m <- p * nz$mu1 + (1 - p) * mu0                     # E[R | H, A] under the tilt
  pid <- 1 / (1 + ((1 - p) / p) * exp(delta * nz$qR) / nz$e[, j])   # P(O = 1 | H, A, R)
  pid <- pmax(pid, 0.05)
  ifelse(d$O == 1, d$R / pid + (1 - 1 / pid) * m, m)
}

## ---------------------------------------------------------------- Step 5: FQE + DR
## importance ratios with the KNOWN behaviour policy
ratios <- function(d, pe1) {
  pb1 <- pib1(d)
  pe <- ifelse(d$A == 1, pe1, 1 - pe1); pb <- ifelse(d$A == 1, pb1, 1 - pb1)
  r  <- pe / pb
  cum  <- ave(r, d$ep, FUN = cumprod)                              # rho_{1:t}
  prev <- ave(r, d$ep, FUN = function(x) c(1, head(cumprod(x), -1)))  # rho_{1:t-1}
  list(cum = cum, prev = prev)
}

## cross-fitted fitted-Q evaluation on pseudo-rewards, backwards in t
fqe <- function(d, Rt, pe1, fold) {
  n <- nrow(d); K <- max(fold); Tm <- max(d$t)
  Qobs <- V <- Vnext <- numeric(n)
  key <- paste(d$ep, d$t)
  for (tt in Tm:1) {
    idx <- which(d$t == tt); if (!length(idx)) next
    y <- Rt[idx] + GAMMA * Vnext[idx]
    for (k in 1:K) {
      tr <- idx[fold[idx] != k]; te <- idx[fold[idx] == k]
      if (!length(te)) next
      if (length(tr) < 20) { Qobs[te] <- mean(y); V[te] <- mean(y); next }
      f  <- rf_fit(Xmat(d[tr]), y[match(tr, idx)])
      q1 <- rf_pred(f, Xmat(d[te], 1)); q0 <- rf_pred(f, Xmat(d[te], 0))
      Qobs[te] <- ifelse(d$A[te] == 1, q1, q0)
      V[te]    <- pe1[te] * q1 + (1 - pe1[te]) * q0
    }
    if (tt > 1) {
      prv <- match(paste(d$ep[idx], tt - 1), key); ok <- !is.na(prv)
      Vnext[prv[ok]] <- V[idx][ok]
    }
  }
  list(Qobs = Qobs, V = V)
}

## Every estimator is returned in the same form: per-participant matrices N and D
## (participants x steps); the estimate is sum_t colSums(N)[t] / colSums(D)[t].
## This covers per-decision IS, its self-normalised form and DR, and lets one
## bootstrap routine serve all estimators.
point <- function(ND) sum(colSums(ND$N) / pmax(colSums(ND$D), 1e-12))

## mean-type estimator: per-episode sums averaged over episodes
nd_mean <- function(d, term) {
  ids <- sort(unique(d$id)); ii <- match(d$id, ids); Tm <- max(d$t)
  N <- matrix(0, length(ids), Tm)
  agg <- data.table(i = ii, t = d$t, n = term)[, .(n = sum(n)), by = .(i, t)]
  N[cbind(agg$i, agg$t)] <- agg$n
  ne <- d[, .(ne = uniqueN(ep)), by = id][match(ids, id), ne]
  list(N = N, D = matrix(ne, length(ids), Tm))
}
## self-normalised (weighted) per-decision estimator
nd_ratio <- function(d, num, den) {
  ids <- sort(unique(d$id)); ii <- match(d$id, ids); Tm <- max(d$t)
  N <- D <- matrix(0, length(ids), Tm)
  agg <- data.table(i = ii, t = d$t, n = num, dd = den)[, .(n = sum(n), dd = sum(dd)), by = .(i, t)]
  N[cbind(agg$i, agg$t)] <- agg$n; D[cbind(agg$i, agg$t)] <- agg$dd
  list(N = N, D = D)
}

## ---------------------------------------------------------------- estimators
## ---------------------------------------------------------------- compute cost
## time (seconds) and peak memory (Mb, from R's garbage collector) of each component,
## so the paper can report cost next to bias and coverage
.clock <- new.env()
timed <- function(label, expr) {
  invisible(gc(reset = TRUE)); t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  g <- gc(); peak <- sum(g[, ncol(g)])                # last column = max used (Mb)
  .clock$rows[[length(.clock$rows) + 1]] <- data.table(component = label,
    secs = proc.time()[["elapsed"]] - t0, peak_mb = peak)
  invisible(val)
}

## ---------------------------------------------------------------- multiple imputation
## kNN-MI: each missing reward is replaced M times by a REAL observed reward drawn from
## its K nearest observed decisions (same availability and action), so the imputations
## carry the whole conditional distribution, not one number (the flaw of mean or single
## imputation). With a tilt, donors are drawn with probability proportional to
## exp(delta(A) q(R_donor)): the pattern-mixture form f(r | missing) ~ f(r | observed)
## exp(delta q(r)) applied to the local neighbourhood. Pooled with Rubin's rules.
mi_draws <- function(d, qc, d0 = 0, dA = 0, M = M_IMP, k = K_IMP) {
  miss <- which(d$O == 0); out <- matrix(NA_real_, length(miss), M)
  if (!length(miss)) return(out)
  X <- as.matrix(d[, FEATS, with = FALSE]); sc <- apply(X, 2, sd); sc[!is.finite(sc) | sc == 0] <- 1
  X <- sweep(X, 2, sc, "/")
  grp <- paste(d$avail, d$A)
  for (g in unique(grp[miss])) {
    m_i <- miss[grp[miss] == g]; don <- which(d$O == 1 & grp == g)
    if (length(don) < k) don <- which(d$O == 1)
    kk <- min(k, length(don))
    nn <- get.knnx(X[don, , drop = FALSE], X[m_i, , drop = FALSE], k = kk)$nn.index
    Rd <- matrix(d$R[don][nn], nrow = length(m_i))            # donor rewards (rows = missing decisions)
    w  <- exp((d0 + dA * d$A[m_i]) * q_fun(Rd, qc)); w <- w / rowSums(w)
    cw <- w; if (kk > 1) for (j in 2:kk) cw[, j] <- cw[, j - 1] + w[, j]
    rows <- match(m_i, miss)
    for (m in 1:M) {
      pick <- pmin(rowSums(cw < runif(length(m_i))) + 1, kk)
      out[rows, m] <- Rd[cbind(seq_along(m_i), pick)]
    }
  }
  out
}
nd_var <- function(ND) {                                     # estimate + wild-type variance
  th  <- colSums(ND$N) / pmax(colSums(ND$D), 1e-12)
  psi <- sweep(ND$N - sweep(ND$D, 2, th, "*"), 2, pmax(colSums(ND$D), 1e-12), "/")
  c(est = sum(th), var = sum(rowSums(psi)^2))
}
## Distribution-free interval for multiple imputation ("MI-Boot pooled", Schomaker & Heumann
## 2018): on EVERY imputed data set, NBOOT wild cluster bootstrap draws (Webb weights, the
## same draws as STAR-DR) are generated around that data set's estimate; all draws are
## pooled and the 2.5% / 97.5% percentiles give the interval. The spread of the pooled draws
## carries both the within- and the between-imputation variability, with no normal or t
## approximation. ND_m: list (one per imputation) of N/D pieces; Wl: list of weight matrices.
mi_boot <- function(ND_m, Wl, level = 0.95) {
  a <- (1 - level) / 2
  ests <- numeric(length(ND_m))
  draws <- unlist(lapply(seq_along(ND_m), function(m) {
    N <- ND_m[[m]]$N; D <- ND_m[[m]]$D
    th <- colSums(N) / pmax(colSums(D), 1e-12); ests[m] <<- sum(th)
    psi <- sweep(N - sweep(D, 2, th, "*"), 2, pmax(colSums(D), 1e-12), "/")
    sum(th) + as.vector(Wl[[m]] %*% rowSums(psi))
  }))
  data.table(est = mean(ests), se = sd(draws),
             lo_w = quantile(draws, a, names = FALSE), hi_w = quantile(draws, 1 - a, names = FALSE))
}
## Rubin's rules (t interval) -- kept ONLY as a comparison column (lo, hi), not the main interval
rubin <- function(ev) {                                      # ev: M x 2 (est, var)
  M <- nrow(ev); Qb <- mean(ev[, 1]); W <- mean(ev[, 2]); B <- var(ev[, 1])
  Tv <- W + (1 + 1 / M) * B
  df <- if (B > 0) (M - 1) * (1 + W / ((1 + 1 / M) * B))^2 else 1e6
  h  <- qt(0.975, df) * sqrt(Tv)
  data.table(est = Qb, se = sqrt(Tv), lo = Qb - h, hi = Qb + h, lo_w = Qb - h, hi_w = Qb + h)
}
mi_label <- function(d0, dA) if (dA != 0) sprintf("kNN-MI-tilt-A(dA=%g)", dA) else if (d0 == 0) "kNN-MI" else "kNN-MI-tilt"
em_label <- function(d0, dA) if (dA != 0) sprintf("EM-cells-tilt-A(dA=%g)", dA) else if (d0 == 0) "EM-cells" else "EM-cells-tilt"
sem_label <- function(d0, dA, K = NULL) paste0(if (dA != 0) sprintf("STAR-EM-A(dA=%g)", dA) else "STAR-EM",
                                               if (!is.null(K)) sprintf(" K=%d", K))

## ---------------------------------------------------------------- EM algorithm
## Nonparametric SELECTION MODEL for the outcome pair (R, M), fitted by EM.
##  * R is coarsened into bins b = 1..B (bin 1 = zero steps, the rest = quantiles of the
##    observed positive steps). Within each context cell c (availability/action x
##    decision time) the reward law is a free HISTOGRAM pi[c, b]: any shape -- skewed,
##    zero-heavy, multi-peaked. No normal or other named distribution.
##  * Missingness: P(M = 1 | c, b) = expit(alpha_c + delta(A_c) * q(rbar_b)), with alpha_c
##    free per cell and delta(A) = delta0 + deltaA * A FIXED (it is not identified).
## Observed-data log-likelihood:
##   l = sum_obs log{ pi[c, b_i] (1 - P(M=1 | c, b_i)) } + sum_miss log{ sum_b pi[c, b] P(M=1 | c, b) }
## E-step: for a missing reward in cell c, w[c, b] = P(bin b | missing, c)
##         proportional to pi[c, b] * P(M=1 | c, b)
## M-step: pi[c, b] = (n_obs[c, b] + n_miss[c] w[c, b]) / n_c ;
##         alpha_c solves sum_b n_c pi[c, b] expit(alpha_c + delta_c q_b) = n_miss[c]
## The log-likelihood is tracked every iteration, so monotone ascent (the EM guarantee)
## is checked, and the number of iterations is reported. With delta fixed the model is
## just-identified, so the maximised likelihood is FLAT in delta: the data cannot choose
## delta, which is why it is a sensitivity parameter (shown in the real-data table).
em_np <- function(d, qc, d0 = 0, dA = 0) {
  obs <- d$O == 1; Ro <- d$R[obs]
  br  <- unique(quantile(Ro[Ro > 0], seq(0, 1, length.out = EM_BINS), names = FALSE))
  binof <- function(r) ifelse(r == 0, 1L, 1L + findInterval(r, br, rightmost.closed = TRUE, all.inside = TRUE))
  B   <- length(br)                                       # 1 zero bin + (length(br) - 1) positive bins
  bin <- rep(NA_integer_, nrow(d)); bin[obs] <- binof(Ro)
  rbar <- vapply(1:B, function(b) { v <- d$R[obs][bin[obs] == b]; if (length(v)) mean(v) else NA_real_ }, 0)
  rbar[is.na(rbar)] <- mean(Ro)
  qb  <- q_fun(rbar, qc)
  cell <- paste(ifelse(!d$avail, "U", d$A), d$t); ucell <- unique(cell)
  cid <- match(cell, ucell); C <- length(ucell)
  dlc <- as.numeric(tapply(d0 + dA * d$A, cid, function(z) z[1]))     # tilt per cell (A is part of the cell)
  nob <- matrix(as.numeric(table(factor(cid[obs], 1:C), factor(bin[obs], 1:B))), C, B)
  nmis <- tabulate(cid[!obs], C); ncell <- rowSums(nob) + nmis
  pi <- (nob + 0.5) / rowSums(nob + 0.5)                  # start: observed histogram (smoothed)
  alpha <- qlogis(pmin(pmax(nmis / ncell, 0.005), 0.995))
  ll_old <- -Inf; lls <- numeric(0)
  for (it in 1:EM_MAXIT) {
    pm <- plogis(alpha + outer(dlc, qb))                  # C x B: P(M = 1 | cell, bin)
    ll <- sum((nob * log(pi * (1 - pm)))[nob > 0]) + sum((nmis * log(rowSums(pi * pm)))[nmis > 0])
    lls <- c(lls, ll)
    if (abs(ll - ll_old) < EM_TOL * (1 + abs(ll))) break
    ll_old <- ll
    w   <- pi * pm; w <- w / rowSums(w)                   # E-step
    tot <- nob + nmis * w                                 # expected complete-data counts
    pi  <- tot / rowSums(tot)                             # M-step: histogram
    for (c in seq_len(C)) alpha[c] <- if (nmis[c] == 0) -15 else          # M-step: missingness
      uniroot(function(a) sum(tot[c, ] * plogis(a + dlc[c] * qb)) - nmis[c], c(-40, 40), tol = 1e-10)$root
  }
  w <- pi * plogis(alpha + outer(dlc, qb)); w <- w / rowSums(w)
  list(w = w, cid = cid, bin = bin, B = B, iters = it, loglik = ll,
       monotone = all(diff(lls) >= -1e-7 * (1 + abs(lls[-1]))))
}
## ---------------------------------------------------------------- STAR-EM (main method)
## STAR-EM: a HIDDEN MARKOV SHARED-PARAMETER MODEL for the outcome pair (R, M), fitted by
## the EM algorithm in its Baum-Welch form (Dempster, Laird & Rubin, 1977; Baum et al., 1970).
##
## Each participant's decisions t = 1..T carry a hidden state S_t in {1..K} (for example
## "active", "low activity", "disengaged / device not worn"):
##   initial state     P(S_1 = k) = pi0_k
##   transitions       P(S_t = l | S_{t-1} = k, A_{t-1}, new day) = Tr_g[k, l]
##                     g = (suggestion sent at the previous decision?) x (new day?):
##                     a notification can move the person to another state
##   reward emission   P(bin of R_t = b | S_t = k, H_t, A_t): probability FOREST with the state
##                     as a feature; within a bin R is drawn from the EMPIRICAL distribution of
##                     observed rewards (no named distribution for steps)
##   missingness       logit P(M_t = 1 | S_t = k, H_t, A_t, R_t)
##                         = alpha(H_t, A_t) + gamma_k + delta(A_t) q(R_t)
## The state drives BOTH the reward and the missingness: given the state the data are
## missing at random, but marginally they are MISSING NOT AT RANDOM (low-activity spells
## are the ones that go missing). With ~200 repeated decisions per person this latent
## structure is identified under rank conditions (Allman, Matias & Rhodes 2009; Gassiat,
## Cleynen & Robin 2016 for nonparametric emissions), so the data -- not an assumed delta --
## carry the main MNAR signal. delta remains available for RESIDUAL MNAR within a state.
##
## E-step  forward-backward (scaled): state posteriors gamma_t(k), pair posteriors xi_t(k, l),
##         and bin posteriors of each missing reward given the state
## M-step  pi0 and Tr in closed form (expected counts); forest and logistic model by
##         weighted MLE on the expected complete data (a generalised EM, monitored for
##         monotone ascent)
## Practical safeguards: several random starts (local maxima), K chosen by held-out
## log-likelihood, states ordered by mean steps (no label switching in the output),
## imputation by forward-filtering backward-sampling of WHOLE state paths (keeps the time
## dependence), proper multiple imputation by bootstrap EM refits, Rubin's rules.

## ---------------------------------------------------------------- missingness model
## No formula and no shape is written down for missingness. P(missing | history, action,
## state) is learned by GRADIENT-BOOSTED TREES (Newton / LogitBoost steps, CART trees from
## the recommended package rpart, so nothing extra to install). The MNAR tilt delta(A) q(R)
## enters as a fixed OFFSET, so the sensitivity analysis is unchanged. Missingness is treated
## as a SPELL PROCESS (device non-wear comes in runs): the features include how long the
## current missing spell has lasted, how long the current observed run has lasted, and the
## person's own missing history so far -- all observed, all computed before the decision.
## With K > 1 hidden states the short-run spell features are left out: the latent state is
## then the spell process, and it must explain the runs (that is what links missingness to
## activity and carries the MNAR signal).
spell_features <- function(d) {
  M <- 1L - d$O; id <- d$id; cap <- STAREM_RUNCAP
  runm <- runo <- integer(length(M)); past <- numeric(length(M))
  for (rows in split(seq_along(M), id)) {             # rows are in time order within a person
    rm <- ro <- 0L; cm <- 0; nb <- 0
    for (j in rows) {
      runm[j] <- min(rm, cap); runo[j] <- min(ro, cap)
      past[j] <- if (nb > 0) cm / nb else NA_real_
      if (M[j] == 1L) { rm <- rm + 1L; ro <- 0L } else { ro <- ro + 1L; rm <- 0L }
      cm <- cm + M[j]; nb <- nb + 1
    }
  }
  past[is.na(past)] <- mean(M)                        # first decision of a person: overall rate
  data.frame(run_missing = runm, run_observed = runo, past_missing_rate = past)
}
miss_base <- function(d) data.frame(t = d$t, A = d$A, avail = as.integer(d$avail), loc_home = d$loc_home,
  loc_work = d$loc_work, sent_today = d$sent_today, temp = d$temp, precip = d$precip, day = d$day)

## Newton boosting of the logit with a fixed offset and case weights, started from the best
## constant. Leaf updates are clipped for stability on rare events.
boost_logit <- function(X, y, w, off, rounds = STAREM_BROUNDS) {
  ## start from the best constant (weighted intercept given the offset; 1-D Newton), not 50%
  b0 <- qlogis(min(max(sum(w * y) / sum(w), 1e-4), 1 - 1e-4)) - weighted.mean(off, w)
  for (it in 1:25) { p <- plogis(b0 + off); g <- sum(w * (y - p)); h <- sum(w * p * (1 - p))
                     b0 <- b0 + g / max(h, 1e-8); if (abs(g / max(h, 1e-8)) < 1e-8) break }
  trees <- list(); f <- rep(b0, nrow(X))
  ctrl <- rpart::rpart.control(maxdepth = STAREM_BDEPTH, cp = 0, xval = 0, minbucket = 50, minsplit = 100,
                               maxcompete = 0, maxsurrogate = 0, usesurrogate = 0)
  for (r in seq_len(rounds)) {
    p <- plogis(off + f); h <- pmax(p * (1 - p), 1e-6)
    df <- X; df$.z <- pmin(pmax((y - p) / h, -20), 20)
    tr <- rpart::rpart(.z ~ ., data = df, weights = w * h, control = ctrl, method = "anova")
    f <- f + STAREM_BLR * pmin(pmax(predict(tr, X), -4), 4)
    trees[[length(trees) + 1]] <- tr
  }
  list(b0 = b0, trees = trees)
}
boost_pred <- function(b, X) {
  f <- rep(b$b0, nrow(X)); for (tr in b$trees) f <- f + STAREM_BLR * pmin(pmax(predict(tr, X), -4), 4); f }

## missingness features of `rows` with the hidden state set to k
mfeat <- function(S, rows, K, k) {
  X <- S$MB[rows, , drop = FALSE]
  if (K == 1) cbind(X, S$MS[rows, , drop = FALSE])
  else { X$past_missing_rate <- S$MS$past_missing_rate[rows]; X$S_state <- factor(k, levels = seq_len(K)); X }
}
## logit P(missing) WITHOUT the tilt, from the cross-fitted ensembles (new data: average of both)
star_mlogit <- function(par, S, rows, k, fold = NULL) {
  X <- mfeat(S, rows, par$K, k)
  if (is.null(fold)) return((boost_pred(par$mfit[[1]], X) + boost_pred(par$mfit[[2]], X)) / 2)
  out <- numeric(length(rows))
  for (f in 1:2) { i <- fold == f; if (any(i)) out[i] <- boost_pred(par$mfit[[f]], X[i, , drop = FALSE]) }
  out
}

star_setup <- function(d, qc) {
  obs <- d$O == 1; Ro <- d$R[obs]
  br  <- unique(quantile(Ro[Ro > 0], seq(0, 1, length.out = STAREM_BINS), names = FALSE)); B <- length(br)
  binof <- function(r) ifelse(r == 0, 1L, 1L + findInterval(r, br, rightmost.closed = TRUE, all.inside = TRUE))
  bin <- rep(NA_integer_, nrow(d)); bin[obs] <- binof(Ro)
  rbar <- vapply(1:B, function(b) { v <- Ro[bin[obs] == b]; if (length(v)) mean(v) else mean(Ro) }, 0)
  list(obs = obs, bin = bin, X = Xmat(d), MB = miss_base(d), MS = spell_features(d),
       lr = ifelse(obs, log1p(d$R), NA_real_), qobs = ifelse(obs, q_fun(d$R, qc), 0),
       qb = q_fun(rbar, qc), rbar = rbar, B = B, don_b = split(Ro, bin[obs]), A = d$A, ep = d$ep)
}

## sequences for a set of rows (rows may repeat = bootstrap copies; seqk labels each copy)
seq_info <- function(S, rows, seqk) {
  sq  <- split(seq_along(rows), factor(seqk, levels = unique(seqk)))
  grp <- integer(length(rows))                              # transition group INTO position t
  for (p in sq) if (length(p) > 1) {
    newday <- S$ep[rows[p[-1]]] != S$ep[rows[p[-length(p)]]]
    grp[p[-1]] <- 1L + S$A[rows[p[-length(p)]]] + 2L * newday
  }
  list(sq = sq, grp = grp)
}

## emission probabilities P(bin | state, H, A) from the CROSS-FITTED forests: a row in fold f is
## predicted by the forest trained on the other fold (participants never predict themselves),
## which keeps the likelihood honest (in-bag forest predictions memorise the data and would
## inflate the log-likelihood, break monotone ascent and favour too many states).
## New data (fold = NULL) get the average of the two forests.
star_P <- function(par, X, B, fold = NULL) {
  pr <- function(f, Xs) { pp <- predict(f, Xs, num.threads = NTHR)$predictions
    out <- matrix(0, nrow(Xs), B); cn <- intersect(colnames(pp), as.character(1:B))
    out[, as.integer(cn)] <- pp[, cn, drop = FALSE]; out }
  if (is.null(fold)) return((pr(par$forest[[1]], X) + pr(par$forest[[2]], X)) / 2)
  P <- matrix(0, nrow(X), B)
  for (f in 1:2) { i <- fold == f; if (any(i)) P[i, ] <- pr(par$forest[[f]], X[i, , drop = FALSE]) }
  P
}
## emission log-likelihoods of every row under every state, and bin posteriors of missing rewards
star_emis <- function(par, S, rows, d0, dA, fold = NULL) {
  K <- par$K; B <- S$B; obs <- S$obs[rows]; po <- which(obs); pm <- which(!obs)
  X <- S$X[rows, , drop = FALSE]; dl <- d0 + dA * S$A[rows]
  loge <- matrix(0, length(rows), K); wmb <- array(0, c(length(pm), K, B)); pmk <- numeric(K)
  for (k in seq_len(K)) {
    X$S_state <- factor(k, levels = seq_len(K))
    P <- star_P(par, X, B, fold)
    eta <- star_mlogit(par, S, rows, k, fold); pmk[k] <- mean(plogis(eta))
    if (length(po)) loge[po, k] <- log(pmax(P[cbind(po, S$bin[rows[po]])], 1e-12)) +
        log(pmax(1 - plogis(eta[po] + dl[po] * S$qobs[rows[po]]), 1e-12))
    if (length(pm)) {
      pk <- P[pm, , drop = FALSE] * plogis(eta[pm] + outer(dl[pm], S$qb))
      loge[pm, k] <- log(pmax(rowSums(pk), 1e-300)); wmb[, k, ] <- pk / pmax(rowSums(pk), 1e-300)
    }
  }
  list(loge = loge, wmb = wmb, po = po, pm = pm, pmk = pmk)
}

## scaled forward-backward over every sequence: state posteriors, expected transitions, log-lik
star_fb <- function(par, em, si) {
  K <- par$K; n <- nrow(em$loge); G <- length(par$Tr)
  gam <- matrix(0, n, K); xi <- lapply(seq_len(G), function(g) matrix(0, K, K)); ll <- 0
  for (p in si$sq) {
    Tn <- length(p); m <- apply(em$loge[p, , drop = FALSE], 1, max)
    e  <- exp(em$loge[p, , drop = FALSE] - m)
    a  <- matrix(0, Tn, K); cs <- numeric(Tn)
    a[1, ] <- par$pi0 * e[1, ]; cs[1] <- sum(a[1, ]); a[1, ] <- a[1, ] / cs[1]
    if (Tn > 1) for (t in 2:Tn) {
      a[t, ] <- as.vector(a[t - 1, ] %*% par$Tr[[si$grp[p[t]]]]) * e[t, ]
      cs[t] <- sum(a[t, ]); a[t, ] <- a[t, ] / cs[t] }
    ll <- ll + sum(log(cs)) + sum(m)
    b <- matrix(1, Tn, K)
    if (Tn > 1) for (t in (Tn - 1):1) {
      g <- si$grp[p[t + 1]]
      b[t, ] <- as.vector(par$Tr[[g]] %*% (e[t + 1, ] * b[t + 1, ])) / cs[t + 1]
      xi[[g]] <- xi[[g]] + (a[t, ] %o% (e[t + 1, ] * b[t + 1, ])) * par$Tr[[g]] / cs[t + 1]
    }
    gg <- a * b; gam[p, ] <- gg / rowSums(gg)
  }
  list(gam = gam, xi = xi, ll = ll)
}

## M-step: closed-form initial/transition probabilities; weighted forest; weighted logistic
star_mstep <- function(par, S, rows, si, em, fb, d0, dA, seed, fold) {
  K <- par$K; B <- S$B; po <- em$po; pm <- em$pm; nm <- length(pm)
  first <- vapply(si$sq, `[`, 1L, 1L)
  pi0 <- colMeans(fb$gam[first, , drop = FALSE]) + 1e-6; pi0 <- pi0 / sum(pi0)
  Tr  <- if (is.null(fb$xi)) par$Tr else lapply(fb$xi, function(x) { x <- x + 1e-3; x / rowSums(x) })
  X <- S$X[rows, , drop = FALSE]; Xo <- X[po, , drop = FALSE]; Xm <- X[rep(pm, B), , drop = FALSE]
  Xe <- do.call(rbind, c(lapply(seq_len(K), function(k) { x <- Xo; x$S_state <- k; x }),
                         lapply(seq_len(K), function(k) { x <- Xm; x$S_state <- k; x })))
  Xe$S_state <- factor(Xe$S_state, levels = seq_len(K))
  ye <- factor(c(rep(S$bin[rows[po]], K), rep(rep(1:B, each = nm), K)), levels = 1:B)
  wmiss <- lapply(seq_len(K), function(k) as.vector(fb$gam[pm, k] * em$wmb[, k, ]))
  we <- c(as.vector(fb$gam[po, , drop = FALSE]), unlist(wmiss))
  keep <- we > 1e-8
  pos_e <- c(rep(po, K), rep(rep(pm, B), K))                      # position behind each expanded row
  forest <- lapply(1:2, function(f) {                             # forest f is trained WITHOUT fold f
    tr <- keep & fold[pos_e] != f
    ## with K > 1 the hidden state is a candidate at EVERY split (otherwise mtry = 3 of 12
    ## features means the forest often never sees it, and the state cannot shape the rewards)
    ranger(x = Xe[tr, , drop = FALSE], y = ye[tr], case.weights = we[tr] / sum(we[tr]),
           probability = TRUE, num.trees = STAREM_TREES, min.node.size = 50,
           always.split.variables = if (K > 1) "S_state" else NULL,
           num.threads = NTHR, seed = seed + f) })                    # same seed every iteration: smooth ascent
  ## missingness: boosted trees on the expected complete data, tilt as a fixed offset,
  ## cross-fitted by the same sequence folds; fixed size, refitted at every M-step
  dl <- d0 + dA * S$A[rows]
  Me <- do.call(rbind, c(lapply(seq_len(K), function(k) mfeat(S, rows[po], K, k)),
                         lapply(seq_len(K), function(k) mfeat(S, rows[rep(pm, B)], K, k))))
  me  <- c(rep(0, length(po) * K), rep(1, nm * B * K))
  off <- c(rep(dl[po] * S$qobs[rows[po]], K), rep(rep(dl[pm], B) * rep(S$qb, each = nm), K))
  mfit <- lapply(1:2, function(f) {
    tr <- keep & fold[pos_e] != f
    boost_logit(Me[tr, , drop = FALSE], me[tr], we[tr], off[tr]) })
  list(K = K, pi0 = pi0, Tr = Tr, forest = forest, mfit = mfit)
}

## one EM run: random start (sticky transitions, random state posteriors) or warm start
## data-driven starting values (the lectures: starting values matter): k-means on each
## decision's local missing fraction and local activity (11-decision window in its sequence)
star_kmeans_start <- function(S, rows, si, K) {
  n <- length(rows); f1 <- f2 <- numeric(n); mlr <- mean(S$lr, na.rm = TRUE)
  wsum <- function(x, h = 5) { m <- length(x); cs <- c(0, cumsum(x)); i <- seq_len(m)
                               cs[pmin(m, i + h) + 1] - cs[pmax(1, i - h)] }
  for (p in si$sq) {
    o <- S$obs[rows[p]]; lr <- ifelse(o, S$lr[rows[p]], 0)
    no <- wsum(o); f1[p] <- 1 - no / wsum(rep(1, length(p)))
    f2[p] <- ifelse(no > 0, wsum(lr) / pmax(no, 1), mlr)
  }
  cl <- kmeans(scale(cbind(f1, f2)), centers = K, nstart = 5)$cluster
  gam <- matrix(0.15 / max(K - 1, 1), n, K); gam[cbind(seq_len(n), cl)] <- 0.85; gam
}

## one EM run: data-driven start ("kmeans"), random start, or warm start (`init`).
## Stopping rule (practical EM with a forest M-step): the M-step refits forests, so the
## log-likelihood jitters around its maximum instead of rising forever; it is tracked
## (monotone-ascent diagnostic) but is NOT the stop signal. EM stops when what the analysis
## uses has settled: the expected mean of the missing rewards (relative change < STAREM_TOL)
## and the state / missingness parameters (largest change < STAREM_PTOL), twice in a row,
## after at least STAREM_MINIT iterations. The imputation model is the average of the last
## STAREM_AVG iterates' E-step quantities (a standard smoothing step for stochastic EM).
star_imean <- function(S, rows, em, fb) {
  if (!length(em$pm)) return(NA_real_)
  K <- ncol(fb$gam); v <- 0
  for (k in seq_len(K)) v <- v + fb$gam[em$pm, k] * as.vector(em$wmb[, k, ] %*% S$rbar)
  mean(v)
}
star_pvec <- function(par, em) {
  tr <- if (par$K > 1) unlist(lapply(par$Tr, diag)) else NULL
  c(em$pmk, par$pi0, tr)                              # P(missing) by state, state shares, stay probabilities
}
star_em_run <- function(S, rows, seqk, d0, dA, K, maxit, init = NULL, seed = 1, start = "kmeans") {
  set.seed(seed); si <- seq_info(S, rows, seqk); n <- length(rows); B <- S$B
  fold <- integer(n); for (i in seq_along(si$sq)) fold[si$sq[[i]]] <- (i %% 2L) + 1L   # 2 folds of sequences
  if (is.null(init)) {
    stick <- matrix(0.2 / max(K - 1, 1), K, K); diag(stick) <- if (K > 1) 0.8 else 1
    par <- list(K = K, pi0 = rep(1 / K, K), Tr = rep(list(stick), 4), mfit = NULL)
    gam <- if (K > 1 && start == "kmeans") star_kmeans_start(S, rows, si, K) else {
      g <- matrix(runif(n * K), n, K); g / rowSums(g) }
    pm <- which(!S$obs[rows]); f0 <- tabulate(S$bin[rows[S$obs[rows]]], B); f0 <- f0 / sum(f0)
    em <- list(po = which(S$obs[rows]), pm = pm, wmb = array(rep(f0, each = length(pm) * K), c(length(pm), K, B)))
    fb <- list(gam = gam, xi = NULL)
  } else {
    par <- init; em <- star_emis(par, S, rows, d0, dA, fold); fb <- star_fb(par, em, si)
  }
  path <- numeric(0); ipath <- numeric(0); small <- 0; hist <- list(); pv_old <- NULL; converged <- FALSE
  for (it in 1:maxit) {
    par <- star_mstep(par, S, rows, si, em, fb, d0, dA, seed, fold)          # M-step
    em  <- star_emis(par, S, rows, d0, dA, fold); fb <- star_fb(par, em, si)   # E-step
    path  <- c(path, fb$ll); im <- star_imean(S, rows, em, fb); ipath <- c(ipath, im)
    pv <- star_pvec(par, em)
    hist[[length(hist) + 1]] <- list(par = par, em = em, fb = fb); if (length(hist) > STAREM_AVG) hist[[1]] <- NULL
    if (it > 1) {
      d_im <- if (is.na(im)) 0 else abs(im - ipath[it - 1]) / max(abs(ipath[it - 1]), 1)
      d_pv <- if (length(pv) == length(pv_old)) max(abs(pv - pv_old)) else Inf
      small <- if (d_im < STAREM_TOL && d_pv < STAREM_PTOL) small + 1 else 0
    }
    pv_old <- pv
    if (it >= STAREM_MINIT && small >= 2) { converged <- TRUE; break }
  }
  ## average the last iterates' E-step quantities (state posteriors, bin posteriors, log-lik)
  last <- hist[[length(hist)]]; L <- length(hist)
  if (L > 1) {
    last$em$wmb <- Reduce(`+`, lapply(hist, function(h) h$em$wmb)) / L
    last$em$loge <- Reduce(`+`, lapply(hist, function(h) h$em$loge)) / L
    last$fb$gam <- Reduce(`+`, lapply(hist, function(h) h$fb$gam)) / L
    last$fb$ll  <- mean(vapply(hist, function(h) h$fb$ll, 0))
  }
  list(par = last$par, em = last$em, fb = last$fb, si = si, path = path, ipath = ipath,
       iters = it, converged = converged)
}

## forward-filtering backward-sampling: draw WHOLE hidden state paths from their posterior
star_ffbs <- function(par, em, si) {
  K <- par$K; s <- integer(nrow(em$loge))
  for (p in si$sq) {
    Tn <- length(p); e <- exp(em$loge[p, , drop = FALSE] - apply(em$loge[p, , drop = FALSE], 1, max))
    a <- matrix(0, Tn, K); a[1, ] <- par$pi0 * e[1, ]; a[1, ] <- a[1, ] / sum(a[1, ])
    if (Tn > 1) for (t in 2:Tn) { a[t, ] <- as.vector(a[t - 1, ] %*% par$Tr[[si$grp[p[t]]]]) * e[t, ]
                                  a[t, ] <- a[t, ] / sum(a[t, ]) }
    s[p[Tn]] <- sample.int(K, 1, prob = a[Tn, ])
    if (Tn > 1) for (t in (Tn - 1):1) {
      w <- a[t, ] * par$Tr[[si$grp[p[t + 1]]]][, s[p[t + 1]]]; s[p[t]] <- sample.int(K, 1, prob = w) }
  }
  s
}

## choose K by held-out log-likelihood (fit on half the participants, score the other half)
star_choose_K <- function(d, qc, Ks = STAREM_KS, d0 = 0, dA = 0) {
  S <- star_setup(d, qc); ids <- unique(d$id)
  set.seed(SEED); half <- sample(rep_len(1:2, length(ids)))
  rbindlist(lapply(Ks, function(K) {
    t0 <- proc.time()[["elapsed"]]
    cv <- sum(vapply(1:2, function(h) {
      tr <- which(d$id %in% ids[half != h]); te <- which(d$id %in% ids[half == h])
      run <- star_em_run(S, tr, d$id[tr], d0, dA, K, STAREM_MAXIT, seed = SEED + K)
      star_fb(run$par, star_emis(run$par, S, te, d0, dA), seq_info(S, te, d$id[te]))$ll }, 0))
    full <- star_em_run(S, seq_len(nrow(d)), d$id, d0, dA, K, STAREM_MAXIT, seed = SEED + K)
    data.table(K = K, heldout_loglik = round(cv, 2), insample_loglik = round(full$fb$ll, 2),
               iterations = full$iters, converged = full$converged, seconds = round(proc.time()[["elapsed"]] - t0, 1))
  }))
}

## interpretable description of the fitted states (real data)
star_state_table <- function(d, qc, K = STAREM_K, d0 = 0, dA = 0) {
  S <- star_setup(d, qc); n <- nrow(d)
  run <- star_em_run(S, seq_len(n), d$id, d0, dA, K, STAREM_MAXIT, seed = SEED)
  X <- S$X
  tab <- rbindlist(lapply(seq_len(K), function(k) {
    X$S_state <- factor(k, levels = seq_len(K))
    P <- star_P(run$par, X, S$B)
    data.table(state = k, share_of_decisions = mean(run$fb$gam[, k]),
               mean_steps_if_observed = as.numeric(mean(P %*% S$rbar)),
               P_zero_steps = mean(P[, 1]),
               P_missing = mean(plogis(star_mlogit(run$par, S, seq_len(n), k))),
               stay_same_day_no_notification = run$par$Tr[[1]][k, k],
               stay_same_day_after_notification = run$par$Tr[[2]][k, k])
  }))
  tab <- tab[order(mean_steps_if_observed)][, state := .I]
  tab[, (2:ncol(tab)) := lapply(.SD, round, 3), .SDcols = 2:ncol(tab)]
  tab
}

## full STAR-EM for one tilt: best of several random starts, then proper multiple imputation
## (bootstrap participants, warm-started EM, FFBS state paths, bins given state, REAL rewards)
star_em <- function(d, qc, d0 = 0, dA = 0, K = STAREM_K) {
  S <- star_setup(d, qc); n <- nrow(d); miss <- which(!S$obs)
  ## several starts only matter when there are hidden labels to permute (K > 1); with K = 1
  ## every start reaches the same fit, so one start is used (saves computing time, same answer)
  starts <- lapply(seq_len(if (K == 1) 1 else STAREM_STARTS), function(s)
    star_em_run(S, seq_len(n), d$id, d0, dA, K, STAREM_MAXIT, seed = SEED + 100 * s,
                start = if (s == 1) "kmeans" else "random"))
  lls  <- vapply(starts, function(z) z$fb$ll, 0)
  best <- starts[[which.max(lls)]]
  rows_by <- split(seq_len(n), d$id); ids <- names(rows_by); si0 <- seq_info(S, seq_len(n), d$id)
  imps <- matrix(NA_real_, length(miss), STAREM_M); biters <- integer(STAREM_M)
  for (m in seq_len(STAREM_M)) {
    sel  <- sample(ids, length(ids), replace = TRUE)
    rows <- unlist(rows_by[sel], use.names = FALSE)
    seqk <- rep(seq_along(sel), lengths(rows_by[sel]))                      # each copy is its own sequence
    fb_m <- star_em_run(S, rows, seqk, d0, dA, K, STAREM_WARM, init = best$par, seed = SEED + 1000 + m)
    biters[m] <- fb_m$iters
    em0 <- star_emis(fb_m$par, S, seq_len(n), d0, dA)                       # original data, refitted model
    st  <- star_ffbs(fb_m$par, em0, si0)                                    # one posterior draw of all states
    pb  <- em0$wmb[cbind(rep(seq_along(miss), S$B), rep(st[miss], S$B), rep(seq_len(S$B), each = length(miss)))]
    pb  <- matrix(pb, length(miss), S$B)                                    # P(bin | state, missing)
    cw  <- pb; if (S$B > 1) for (b in 2:S$B) cw[, b] <- cw[, b - 1] + pb[, b]
    bb  <- pmin(rowSums(cw < runif(length(miss)) * cw[, S$B]) + 1L, S$B)
    imps[, m] <- vapply(bb, function(b) { v <- S$don_b[[as.character(b)]]; v[sample.int(length(v), 1)] }, 0)
  }
  inc <- diff(best$path)
  list(imp = imps, iters = best$iters, loglik = best$fb$ll, converged = best$converged, K = K,
       pct_increasing = if (length(inc)) 100 * mean(inc >= -1e-6 * (1 + abs(best$path[-1]))) else 100,
       start_spread = round(diff(range(lls)), 3), weights = paste(round(best$par$pi0, 3), collapse = "/"),
       boot_iters = mean(biters), mean_imputed = mean(imps))
}

## ---------------------------------------------------------------- EM versus direct MLE
## The lectures' point: classical (direct) maximisation of the observed-data likelihood
## is hard when a layer is missing. Same model (the EM-cells selection model), two
## optimisers: EM vs direct numerical MLE (optim, BFGS, softmax/logit parametrisation).
## Reports convergence, iterations, time and the likelihood each one reaches.
em_vs_mle <- function(d, qc, d0 = 0) {
  e0 <- proc.time()[["elapsed"]]; em <- em_np(d, qc, d0, 0); t_em <- proc.time()[["elapsed"]] - e0
  obs <- d$O == 1; Ro <- d$R[obs]
  br  <- unique(quantile(Ro[Ro > 0], seq(0, 1, length.out = EM_BINS), names = FALSE)); B <- length(br)
  bin <- rep(NA_integer_, nrow(d))
  bin[obs] <- ifelse(Ro == 0, 1L, 1L + findInterval(Ro, br, rightmost.closed = TRUE, all.inside = TRUE))
  rbar <- vapply(1:B, function(b) { v <- Ro[bin[obs] == b]; if (length(v)) mean(v) else mean(Ro) }, 0)
  qb <- q_fun(rbar, qc)
  cid <- match(paste(ifelse(!d$avail, "U", d$A), d$t), unique(paste(ifelse(!d$avail, "U", d$A), d$t))); C <- max(cid)
  nob <- matrix(as.numeric(table(factor(cid[obs], 1:C), factor(bin[obs], 1:B))), C, B)
  nmis <- tabulate(cid[!obs], C)
  negll <- function(th) {
    a <- th[1:C]; L <- cbind(0, matrix(th[-(1:C)], C, B - 1))
    pi <- exp(L - apply(L, 1, max)); pi <- pi / rowSums(pi)
    pm <- plogis(a + outer(rep(d0, C), qb))
    -(sum((nob * log(pmax(pi * (1 - pm), 1e-300)))[nob > 0]) + sum((nmis * log(pmax(rowSums(pi * pm), 1e-300)))[nmis > 0]))
  }
  m0 <- proc.time()[["elapsed"]]
  op <- tryCatch(optim(rep(0, C + C * (B - 1)), negll, method = "BFGS", control = list(maxit = 1000)),
                 error = function(e) list(value = NA, convergence = 99, counts = c(NA, NA)))
  t_ml <- proc.time()[["elapsed"]] - m0
  data.table(delta = d0, parameters = C + C * (B - 1),
             method = c("EM algorithm", "Direct MLE (optim BFGS)"),
             loglik = round(c(em$loglik, -op$value), 3),
             converged = c(TRUE, op$convergence == 0),
             iterations = c(em$iters, op$counts[1]), seconds = round(c(t_em, t_ml), 2))
}

## multiple imputation from the EM posterior: draw the bin of each missing reward from
## w[c, ], then a REAL observed step count from that bin (same cell if possible)
em_draws <- function(d, em, M = M_IMP) {
  miss <- which(d$O == 0); out <- matrix(NA_real_, length(miss), M)
  if (!length(miss)) return(out)
  obs <- which(d$O == 1)
  don_cb <- split(d$R[obs], paste(em$cid[obs], em$bin[obs])); don_b <- split(d$R[obs], em$bin[obs])
  cw <- em$w; if (em$B > 1) for (b in 2:em$B) cw[, b] <- cw[, b - 1] + em$w[, b]
  cm <- em$cid[miss]
  for (m in 1:M) {
    bb <- pmin(rowSums(cw[cm, , drop = FALSE] < runif(length(miss))) + 1L, em$B)
    out[, m] <- mapply(function(c, b) {
      v <- don_cb[[paste(c, b)]]; if (is.null(v)) v <- don_b[[as.character(b)]]
      v[sample.int(length(v), 1)] }, cm, bb)
  }
  out
}

## d must contain the analyst-observed data; returns a named list of ND objects
## (names "estimator|delta|policy").
## DS = data.table(d0, dA): the assumed tilts. Labels: MAR-DR (0, 0), STAR-DR (d0, 0),
## STAR-DR-A(dA=x) (d0, x). Distribution-free bounds are added for every policy.
ds_label <- function(d0, dA) if (dA != 0) sprintf("STAR-DR-A(dA=%g)", dA) else if (d0 == 0) "MAR-DR" else "STAR-DR"
estimate_all <- function(d, qc, DS = data.table(d0 = 0, dA = 0), policies = names(POLICIES), K = 5,
                         R_true = NULL, bounds = TRUE, sem_K = STAREM_K) {
  stopifnot(!is.unsorted(d$ep), all(d[, !is.unsorted(t), by = ep]$V1))
  .clock$rows <- list()
  fold <- make_folds(d$id, K)
  nz   <- timed("Nuisance forests (shared by IPW, RF-Impute, DR)", nuisances(d, DS, qc, fold))
  out  <- list(); bl <- list(); bu <- list(); mi_nd <- list()
  r_hi <- quantile(d$R[d$O == 1], BOUND_Q, names = FALSE)      # rewards assumed in [0, r_hi]
  imp  <- list()
  timed(sprintf("kNN-MI: drawing %d imputations", M_IMP),
        for (j in seq_len(nrow(DS))) imp[[j]] <- mi_draws(d, qc, DS$d0[j], DS$dA[j]))
  emfit <- list(); imp_em <- list()
  for (j in seq_len(nrow(DS))) timed("EM: fitting the selection model (per tilt)",
        emfit[[j]] <- em_np(d, qc, DS$d0[j], DS$dA[j]))
  timed(sprintf("EM-MI: drawing %d imputations", M_IMP),
        for (j in seq_len(nrow(DS))) imp_em[[j]] <- em_draws(d, emfit[[j]]))
  semfit <- list()
  for (j in seq_len(nrow(DS))) semfit[[j]] <- lapply(sem_K, function(K)
        timed(sprintf("STAR-EM (K = %d): EM fit + %d bootstrap EM refits (per tilt)", K, STAREM_M),
              star_em(d, qc, DS$d0[j], DS$dA[j], K)))
  klab <- function(K) if (length(sem_K) > 1) K else NULL           # label K only when several are fitted
  miss <- which(d$O == 0)
  for (pn in policies) {
    pe1 <- POLICIES[[pn]](d); rr <- ratios(d, pe1)
    Robs <- ifelse(d$O == 1, d$R, 0)
    ## 1. complete-case weighted PDIS (ignores missingness)
    timed("CC-WPDIS", out[[paste("CC-WPDIS", 0, pn, sep = "|")]] <-
      nd_ratio(d, rr$cum * d$O * Robs, rr$cum * d$O))
    ## 2. MAR inverse-probability weighting
    timed("MAR-IPW", out[[paste("MAR-IPW", 0, pn, sep = "|")]] <-
      nd_mean(d, rr$cum * d$O * Robs / nz$p))
    ## 3. single random-forest imputation + PDIS (one number per missing value)
    timed("RF-Impute (single)", out[[paste("RF-Impute", 0, pn, sep = "|")]] <-
      nd_mean(d, rr$cum * ifelse(d$O == 1, d$R, nz$mu1)))
    ## 3b. kNN multiple imputation + PDIS, one ND per imputation (pooled below)
    timed("kNN-MI: analysing M datasets", for (j in seq_len(nrow(DS))) {
      lab <- mi_label(DS$d0[j], DS$dA[j])
      key <- paste(lab, DS$d0[j], sep = "@")
      mi_nd[[key]][[pn]] <- lapply(seq_len(M_IMP), function(m) {
        Rc <- Robs; Rc[miss] <- imp[[j]][, m]; nd_mean(d, rr$cum * Rc) })
      attr(mi_nd[[key]], "d0") <- DS$d0[j]
    })
    ## 3d. STAR-EM (main method): proper multiple imputation from the semiparametric EM + PDIS
    timed("STAR-EM: analysing M datasets", for (j in seq_len(nrow(DS))) for (sf in semfit[[j]]) {
      lab <- sem_label(DS$d0[j], DS$dA[j], klab(sf$K))
      key <- paste(lab, DS$d0[j], sep = "@")
      mi_nd[[key]][[pn]] <- lapply(seq_len(STAREM_M), function(m) {
        Rc <- Robs; Rc[miss] <- sf$imp[, m]; nd_mean(d, rr$cum * Rc) })
      attr(mi_nd[[key]], "d0") <- DS$d0[j]
    })
    ## 3c. EM (cells) multiple imputation + PDIS
    timed("EM-MI: analysing M datasets", for (j in seq_len(nrow(DS))) {
      lab <- em_label(DS$d0[j], DS$dA[j])
      key <- paste(lab, DS$d0[j], sep = "@")
      mi_nd[[key]][[pn]] <- lapply(seq_len(M_IMP), function(m) {
        Rc <- Robs; Rc[miss] <- imp_em[[j]][, m]; nd_mean(d, rr$cum * Rc) })
      attr(mi_nd[[key]], "d0") <- DS$d0[j]
    })
    ## 4. STAR-DR for each assumed tilt (0, 0 is the MAR doubly robust estimator)
    for (j in seq_len(nrow(DS))) timed("DR: fitted-Q + pseudo-reward (per tilt)", {
      Rt <- pseudo_reward(d, nz, j, DS$d0[j], DS$dA[j])
      fq <- fqe(d, Rt, pe1, fold)
      term <- (GAMMA^(d$t - 1)) * (rr$cum * (Rt - fq$Qobs) + rr$prev * fq$V)
      out[[paste(ds_label(DS$d0[j], DS$dA[j]), DS$d0[j], pn, sep = "|")]] <- nd_mean(d, term)
    })
    ## 5. worst-case bounds: every missing reward set to the smallest (0) or the largest
    ##    plausible value (r_hi). The ratios are non-negative, so these bracket the
    ##    policy value under ANY joint law of (R, O) -- no tilt, no model.
    if (bounds) {
      bl[[pn]] <- nd_mean(d, rr$cum * ifelse(d$O == 1, d$R, 0))
      bu[[pn]] <- nd_mean(d, rr$cum * ifelse(d$O == 1, d$R, r_hi))
      out[[paste("Bound-lower", 0, pn, sep = "|")]] <- bl[[pn]]
      out[[paste("Bound-upper", 0, pn, sep = "|")]] <- bu[[pn]]
    }
    ## 6. oracle (simulation only): PDIS with the true, complete rewards
    if (!is.null(R_true))
      out[[paste("Oracle-PDIS", 0, pn, sep = "|")]] <- nd_mean(d, rr$cum * R_true)
  }
  ## bounds for contrasts a - b: [lower(a) - upper(b), upper(a) - lower(b)]
  if (bounds) for (cc in CONTRASTS) if (all(cc %in% names(bl))) {
    lab <- paste0(cc[1], "-", cc[2])
    out[[paste("Bound-lower", 0, lab, sep = "|")]] <- list(N = bl[[cc[1]]]$N - bu[[cc[2]]]$N, D = bl[[cc[1]]]$D)
    out[[paste("Bound-upper", 0, lab, sep = "|")]] <- list(N = bu[[cc[1]]]$N - bl[[cc[2]]]$N, D = bl[[cc[1]]]$D)
  }
  ## pool the multiple imputations, for policies and contrasts.
  ## Main interval (lo_w, hi_w): pooled wild-bootstrap percentiles (distribution-free).
  ## Comparison only (lo, hi): Rubin's rules with a t reference distribution.
  n_id <- uniqueN(d$id)
  Wl <- lapply(seq_len(max(M_IMP, STAREM_M)), function(m) matrix(webb(NBOOT * n_id), NBOOT, n_id))
  pool_one <- function(ND_m) {
    b <- mi_boot(ND_m, Wl); r <- rubin(t(vapply(ND_m, nd_var, numeric(2))))
    data.table(est = b$est, se = b$se, lo = r$lo, hi = r$hi, lo_w = b$lo_w, hi_w = b$hi_w)
  }
  mi_rows <- list()
  for (key in names(mi_nd)) {
    lab <- sub("@.*", "", key); d0 <- attr(mi_nd[[key]], "d0"); L <- mi_nd[[key]]
    for (pn in intersect(policies, names(L)))
      mi_rows[[length(mi_rows) + 1]] <- cbind(data.table(estimator = lab, delta = d0, policy = pn), pool_one(L[[pn]]))
    for (cc in CONTRASTS) if (all(cc %in% names(L)))       # paired: same weights for both policies
      mi_rows[[length(mi_rows) + 1]] <- cbind(data.table(estimator = lab, delta = d0, policy = paste0(cc[1], "-", cc[2])),
        pool_one(lapply(seq_along(L[[cc[1]]]), function(m)
          list(N = L[[cc[1]]][[m]]$N - L[[cc[2]]][[m]]$N, D = L[[cc[1]]][[m]]$D))))
  }
  attr(out, "mi") <- rbindlist(mi_rows)
  attr(out, "em") <- rbindlist(lapply(seq_len(nrow(DS)), function(j) rbindlist(lapply(semfit[[j]], function(sf) data.table(
    d0 = DS$d0[j], dA = DS$dA[j], iterations = emfit[[j]]$iters,
    max_loglik = round(emfit[[j]]$loglik, 4), monotone_ascent = emfit[[j]]$monotone,
    starem_K = sf$K, starem_converged = sf$converged,
    starem_iterations = sf$iters, starem_loglik = round(sf$loglik, 3),
    starem_pct_increasing = round(sf$pct_increasing, 1),
    starem_boot_iterations = sf$boot_iters, starem_mean_imputed = round(sf$mean_imputed, 1),
    starem_start_spread = sf$start_spread, starem_weights = sf$weights)))))
  attr(out, "timing") <- rbindlist(.clock$rows)[, .(secs = sum(secs), peak_mb = max(peak_mb)), by = component]
  out
}

## ---------------------------------------------------------------- Step 6: inference
## Two distribution-free intervals, neither refits the nuisance models:
##  (a) participant (cluster) bootstrap   (b) wild cluster bootstrap, Webb weights
webb <- function(n) sample(c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5)), n, TRUE)

infer <- function(ND_list, B = 999, level = 0.95) {
  nm  <- names(ND_list); n_id <- nrow(ND_list[[1]]$N)
  C   <- t(rmultinom(B, n_id, rep(1 / n_id, n_id)))        # B x participants
  W   <- matrix(webb(B * n_id), B, n_id)
  a   <- (1 - level) / 2
  draws_boot <- draws_wild <- list(); pts <- c()
  for (k in nm) {
    N <- ND_list[[k]]$N; D <- ND_list[[k]]$D
    th <- colSums(N) / pmax(colSums(D), 1e-12); pts[k] <- sum(th)
    draws_boot[[k]] <- rowSums((C %*% N) / pmax(C %*% D, 1e-12))
    psi <- sweep(N - sweep(D, 2, th, "*"), 2, pmax(colSums(D), 1e-12), "/")  # influence
    draws_wild[[k]] <- pts[k] + as.vector(W %*% rowSums(psi))
  }
  row <- function(lab, est, db, dw) data.table(
    target = lab, est = est, se = sd(db),
    lo = 2 * est - quantile(db, 1 - a, names = FALSE), hi = 2 * est - quantile(db, a, names = FALSE),
    lo_w = quantile(dw, a, names = FALSE), hi_w = quantile(dw, 1 - a, names = FALSE))
  res <- rbindlist(lapply(nm, function(k) row(k, pts[k], draws_boot[[k]], draws_wild[[k]])))
  ## policy contrasts, using the same bootstrap draws (paired)
  for (k in nm) {
    s <- strsplit(k, "|", fixed = TRUE)[[1]]
    if (startsWith(s[1], "Bound")) next              # bound contrasts are built directly
    for (cc in CONTRASTS) if (s[3] == cc[1]) {
      k2 <- paste(s[1], s[2], cc[2], sep = "|")
      if (!k2 %in% nm) next
      res <- rbind(res, row(paste(s[1], s[2], paste0(cc[1], "-", cc[2]), sep = "|"),
                            pts[k] - pts[k2], draws_boot[[k]] - draws_boot[[k2]],
                            draws_wild[[k]] - draws_wild[[k2]]))
    }
  }
  res[, c("estimator", "delta", "policy") := tstrsplit(target, "|", fixed = TRUE)]
  res[, target := NULL]
  res[, delta := as.numeric(delta)]
  setcolorder(res, c("estimator", "delta", "policy"))
  res[]
}

## =====================================================================
## 3. Distribution-free simulation generator
## =====================================================================
## 03_generator.R
## Distribution-free plasmode ("digital twin") generator built from HeartSteps.
##  * Participants are resampled with their REAL context histories
##    (availability, location, weather, time of day, study day).
##  * Rewards are REAL step counts: each simulated reward is drawn from the k
##    real decisions most similar in context, history and action
##    (nearest-neighbour conditional resampling). No reward distribution is written down.
##  * Missingness is imposed with known mechanisms, so the truth is known.


KNN_FEATS <- c("t", "loc_home", "loc_work", "temp", "precip", "sent_today", "lagq", "day")
K_NN <- 20

make_pool <- function(panel) {
  pool <- panel[O == 1]
  X  <- as.matrix(pool[, KNN_FEATS, with = FALSE])
  sc <- apply(X, 2, sd); sc[!is.finite(sc) | sc == 0] <- 1
  X  <- sweep(X, 2, sc, "/")
  skel <- panel[, .(id, day, t, avail, loc_home, loc_work, temp, precip, pm_full,
                    A_real = A, M_real = 1L - O)]
  setorder(skel, id, day, t)
  ## three separate donor pools: unavailable (already active / driving), available
  ## but not sent, and sent -- unavailable decisions have far higher step counts
  g1 <- pool$A == 1; g0 <- pool$A == 0 & pool$avail; gu <- !pool$avail
  list(X1 = X[g1, , drop = FALSE], R1 = pool$R[g1],
       X0 = X[g0, , drop = FALSE], R0 = pool$R[g0],
       XU = X[gu, , drop = FALSE], RU = pool$R[gu],
       sc = sc, skel = skel, rows = split(seq_len(nrow(skel)), skel$id),
       qc = attr(panel, "qc"))
}

draw_reward <- function(Xq, a, avail, pool) {
  R <- numeric(nrow(Xq))
  grp <- ifelse(!avail, "U", ifelse(a == 1, "1", "0"))
  for (gg in c("U", "0", "1")) {
    s <- which(grp == gg); if (!length(s)) next
    Xd <- pool[[paste0("X", gg)]]; Rd <- pool[[paste0("R", gg)]]
    nn <- get.knnx(Xd, Xq[s, , drop = FALSE], k = K_NN)$nn.index
    R[s] <- Rd[nn[cbind(seq_along(s), sample.int(K_NN, length(s), TRUE))]]
  }
  R
}

## Generate a cohort of N resampled participants.
## policy = NULL -> HeartSteps behaviour policy (P(send) = 0.6 when available);
## otherwise a function from POLICIES (used to compute the true value).
gen_cohort <- function(N, pool, policy = NULL, hidden = FALSE) {
  ids <- sample(names(pool$rows), N, replace = TRUE)
  rr  <- pool$rows[ids]
  d   <- pool$skel[unlist(rr)]
  d[, sid := rep(seq_len(N), lengths(rr))]
  d[, id := sid]
  setorder(d, sid, day, t)
  d[, ep := .GRP, by = .(sid, day)]
  d[, `:=`(A = 0L, R = 0, sent_today = 0, lagq = 0, S_hidden = 0L)]
  key <- paste(d$ep, d$t)
  if (hidden) {                                       # state at the start of each day: Markov over days
    eps <- unique(d[, .(sid, ep)]); eps[, dr := seq_len(.N), by = sid]
    S0 <- integer(nrow(eps)); f <- which(eps$dr == 1)
    S0[f] <- 1L + rbinom(length(f), 1, 0.5)            # 1 = low activity, 2 = high activity
    for (r in seq_len(max(eps$dr))[-1]) { i <- which(eps$dr == r)
      S0[i] <- ifelse(runif(length(i)) < HID_STAY_DAY, S0[i - 1], 3L - S0[i - 1]) }
    S0 <- S0[match(d$ep, eps$ep)]
  }
  for (tt in 1:5) {
    idx <- which(d$t == tt); if (!length(idx)) next
    if (tt > 1) {                                     # carry within-day history forward
      prv <- match(paste(d$ep[idx], tt - 1), key); ok <- !is.na(prv)
      set(d, idx[ok], "sent_today", d$sent_today[prv[ok]] + d$A[prv[ok]])
      set(d, idx[ok], "lagq", q_fun(d$R[prv[ok]], pool$qc))
    }
    if (is.null(policy) && REAL_A) {
      set(d, idx, "A", d$A_real[idx])                 # replay the trial's real randomisation
    } else {
      p1 <- if (is.null(policy)) ifelse(d$avail[idx], P_SEND, 0) else policy(d[idx])
      set(d, idx, "A", as.integer(rbinom(length(idx), 1, p1)))
    }
    Xq <- sweep(as.matrix(d[idx, KNN_FEATS, with = FALSE]), 2, pool$sc, "/")
    if (!hidden) { set(d, idx, "R", draw_reward(Xq, d$A[idx], d$avail[idx], pool)); next }
    ## hidden activity state of this decision: day-start state at t = 1, then a Markov step
    ## from the previous decision; a notification at the previous decision can lift a low state
    s <- S0[idx]
    if (tt > 1) {
      prv <- match(paste(d$ep[idx], tt - 1), key); ok <- which(!is.na(prv)); sp <- d$S_hidden[prv[ok]]
      up <- HID_UP0 + HID_UPA * d$A[prv[ok]]; u <- runif(length(ok))
      s[ok] <- ifelse(sp == 1L, ifelse(u < up, 2L, 1L), ifelse(u < HID_DOWN, 1L, 2L))
    }
    set(d, idx, "S_hidden", s)
    r1 <- draw_reward(Xq, d$A[idx], d$avail[idx], pool); r2 <- draw_reward(Xq, d$A[idx], d$avail[idx], pool)
    set(d, idx, "R", ifelse(s == 1L, pmin(r1, r2), pmax(r1, r2)))
  }
  d[]
}

## ---------------------------------------------------------------- missingness
## The missing indicator is the second half of the outcome pair (R, O). Its law is
## LEARNED from HeartSteps instead of written down:
##   P(M = 1 | context, length of the current missing spell, A, person) = probability forest fitted to the
##   real data (features below), using each participant's leave-one-out missing rate,
## so device-non-wear spells, time-of-day / study-day trends and person-to-person
## heterogeneity are those of the real study. MNAR is added as a tilt on the logit:
##   logit P(M = 1 | ..., R) = c + logit(p_forest) + (delta0 + deltaA * A) * q_gen(R)
## c is calibrated to the target missing rate. q_gen is "linear" (= the q used by
## STAR-DR) or "threshold" (sedentary vs active, to test a misspecified tilt shape).
fit_miss_model <- function(panel) {
  X <- as.data.frame(panel[, .(t, day, loc_home, loc_work, temp, precip, sent_today, A,
                               avail = as.integer(avail), prevrun, pm = pm_loo)])
  ranger(x = X, y = factor(1L - panel$O, levels = 0:1), probability = TRUE, num.trees = 300,
         min.node.size = 25, num.threads = NTHR, seed = 100)
}
q_gen <- function(R, qc, shape = "linear")
  if (shape == "threshold") ifelse(log1p(R) >= qc[1], 1, -1) else q_fun(R, qc)

learned_missing <- function(d, mm, rate, d0 = 0, dA = 0, shape = "linear", qc, extra = 0) {
  pos <- d[, seq_len(.N), by = sid]$V1                   # d is sorted by sid, day, t
  X <- as.data.frame(d[, .(t, day, loc_home, loc_work, temp, precip, sent_today, A,
                           avail = as.integer(avail), prevrun = 0L, pm = pm_full)])
  lg <- function(p) qlogis(pmin(pmax(p, 0.005), 0.995))
  ## logit of P(missing) for every possible length of the current missing spell (0 = previous observed)
  L <- vapply(0:10, function(r) { X$prevrun <- r; lg(predict(mm, X, num.threads = NTHR)$predictions[, "1"]) },
              numeric(nrow(d)))
  tilt <- (d0 + dA * d$A) * q_gen(d$R, qc, shape) + extra   # extra: hidden-state term (scenario 7)
  U <- runif(nrow(d)); P <- max(pos)
  run <- function(c) {
    M <- integer(nrow(d)); rl <- integer(nrow(d)); i <- which(pos == 1)
    M[i] <- as.integer(U[i] < plogis(c + L[cbind(i, 1L)] + tilt[i]))
    for (p in 2:P) { i <- which(pos == p)
      rl[i] <- ifelse(M[i - 1] == 1L, pmin(rl[i - 1] + 1L, 10L), 0L)
      M[i] <- as.integer(U[i] < plogis(c + L[cbind(i, rl[i] + 1L)] + tilt[i])) }
    M
  }
  c0 <- uniroot(function(c) mean(run(c)) - rate, c(-20, 20), tol = 1e-3)$root
  run(c0)
}

## mech: "MCAR" (benchmark) or "learned" (the forest + tilt above; d0 = dA = 0 is MAR)
## "transplant": no random missingness draws at all. Each simulated participant keeps
## the REAL missing/observed pattern of the participant they were resampled from
## (real spells, real rate, real time-of-day structure). For MNAR the whole pattern is
## shifted by s whole days (cyclically), choosing the shift with probability
##   P(s)  proportional to  exp( sum_t M_t(s) * (delta0 + deltaA * A_t) * q_gen(R_t) ),
## i.e. an exponential tilt on the JOINT law of the binary missingness vector and the
## rewards: low-step days are more likely to fall inside real non-wear spells when
## delta < 0. delta = 0 gives a uniformly random shift (patterns unrelated to R).
transplant_missing <- function(d, d0 = 0, dA = 0, shape = "linear", qc) {
  tilt <- (d0 + dA * d$A) * q_gen(d$R, qc, shape)
  M <- integer(nrow(d))
  for (rows in split(seq_len(nrow(d)), d$sid)) {
    L <- length(rows); m0 <- d$M_real[rows]
    shifts <- seq(0, L - 1, by = 5)                    # whole days keep the time-of-day slot
    cand <- vapply(shifts, function(k) m0[((seq_len(L) - 1 + k) %% L) + 1], numeric(L))
    sc <- as.vector(crossprod(cand, tilt[rows]))       # sum_t M_t(s) * tilt_t for each shift
    pr <- exp(sc - max(sc)); pr <- pr / sum(pr)
    M[rows] <- cand[, sample.int(length(shifts), 1, prob = pr)]
  }
  M
}

apply_missing <- function(d, mech, rate, d0 = 0, dA = 0, shape = "linear", qc, mm, extra = 0) {
  M <- switch(mech,
    "MCAR"       = rbinom(nrow(d), 1, rate),
    "transplant" = transplant_missing(d, d0, dA, shape, qc),
    learned_missing(d, mm, rate, d0, dA, shape, qc, extra))
  d[, R_true := R]
  d[, O := 1L - as.integer(M)]
  d[O == 0L, R := NA_real_]
  ## the analyst's history only uses OBSERVED rewards
  d[, lagobs := shift(O, fill = 0L), by = ep]
  d[, lagq := shift(fifelse(O == 1L, q_fun(R, qc), 0), fill = 0), by = ep]
  d[, prevobs := shift(O, fill = 1L), by = sid]
  d[]
}

## true policy values by large-scale Monte Carlo (complete rewards, target policy acting)
true_values <- function(pool, N = 5000, policies = POLICIES, hidden = FALSE) {
  rbindlist(lapply(names(policies), function(pn) {
    d <- gen_cohort(N, pool, policies[[pn]], hidden)
    ep_sum <- d[, .(v = sum(R)), by = .(sid, ep)]
    pid <- ep_sum[, .(v = sum(v), ne = .N), by = sid]
    data.table(policy = pn, truth = mean(ep_sum$v),
               mcse = sd(pid$v / pid$ne) / sqrt(N))
  }))
}

## =====================================================================
## 4. RUN
## =====================================================================
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
prep_file <- file.path(out_dir, "prep.rds")

## scenarios: the 7 desktop scenarios (N = 37, the trial's size) + the hidden-state scenario at N = 150
grid <- data.table(
  scenario   = c("MAR (learned)", "MNAR (learned)", "MNAR, action-dependent", "MNAR, wrong tilt shape",
                 "MAR (real patterns)", "MNAR (real patterns)", "MNAR, hidden activity states",
                 sprintf("MNAR, hidden activity states (N = %d)", N_BIG)),
  mech       = c(rep("learned", 4), rep("transplant", 2), "learned", "learned"), rate = 0.30,
  N          = c(rep(37L, 7), as.integer(N_BIG)),
  d0_true    = c(0, -1, -0.5, -1, 0, -1, 0, 0),
  dA_true    = c(0, 0, -1, 0, 0, 0, 0, 0),
  shape      = c("linear", "linear", "linear", "threshold", "linear", "linear", "linear", "linear"),
  world      = c(rep("base", 6), "hidden", "hidden"))
grid[, cond := .I]
## tilts the analyst assumes: MAR, the true tilt ignoring the action, the true action-dependent
## tilt, and (hidden-state scenarios) a delta = -1 sensitivity tilt
assumed <- function(g) unique(rbind(data.table(d0 = 0, dA = 0),
  if (g$d0_true != 0 || g$dA_true != 0) data.table(d0 = g$d0_true, dA = 0),
  if (g$dA_true != 0) data.table(d0 = g$d0_true, dA = g$dA_true),
  if (g$world == "hidden") data.table(d0 = -1, dA = 0)))
## array tasks: N = 37 scenarios first, the N = 150 scenario last (can be submitted separately)
tasks <- rbindlist(lapply(grid$cond, function(c) {
  big <- grid$N[c] > 37; ch <- if (big) CHUNK_BIG else CHUNK; nr <- if (big) REPS_BIG else REPS
  s <- seq(1, nr, by = ch)
  data.table(cond = c, r1 = s, r2 = pmin(s + ch - 1, nr)) }))
tasks[, task := .I]

## one replication of scenario cnd (identical to the desktop loop body)
run_rep <- function(cnd, r) {
  g <- grid[cond == cnd]
  set.seed(100 + 100000 * cnd + r)
    t0 <- proc.time()[["elapsed"]]
  out <- tryCatch({
      hid <- g$world == "hidden"
      d <- gen_cohort(g$N, pool, hidden = hid)
      d <- apply_missing(d, g$mech, g$rate, g$d0_true, g$dA_true, g$shape, qc, mm,
                         extra = if (hid) HID_GAM * (d$S_hidden == 1L) else 0)
      nd <- estimate_all(d, qc, DS = assumed(g), R_true = d$R_true,
                         sem_K = if (hid) sort(unique(c(1L, 2L, STAREM_K))) else STAREM_K)
      x <- rbind(infer(nd, B = NBOOT), attr(nd, "mi"), use.names = TRUE)
      cost_rows[[length(cost_rows) + 1]] <<- attr(nd, "timing")[, `:=`(cond = cnd, rep = r)]
      em_rows[[length(em_rows) + 1]] <<- attr(nd, "em")[, `:=`(cond = cnd, rep = r)]
      x[, `:=`(cond = cnd, rep = r, miss = 1 - mean(d$O), secs = proc.time()[["elapsed"]] - t0)]
    }, error = function(e) { say("  rep ", r, " failed: ", conditionMessage(e)); NULL })
  out
}

## ---------------------------------------------------------------- prep
if (MODE == "prep") {
panel <- prepare_panel(data_file); qc <- attr(panel, "qc")
desc <- data.table(
  quantity = c("participants", "episodes (days)", "decisions", "available",
               "P(send | available)", "rewards missing", "zeros among observed",
               "median steps (observed)", "mean steps (observed)", "skewness (observed)",
               "P(missing | previous missing)", "P(missing | previous observed)"),
  value = with(panel, c(uniqueN(id), uniqueN(ep), length(id), mean(avail), mean(A[avail]), 1 - mean(O),
    mean(R[O == 1] == 0), median(R[O == 1]), mean(R[O == 1]),
    { y <- R[O == 1]; mean((y - mean(y))^3) / sd(y)^3 },
    mean(O[prevobs == 0] == 0), mean(O[prevobs == 1] == 0))))
stopifnot(nrow(desc) == 12)
desc[, value := round(value, 3)]
print(desc); fwrite(desc, file.path(out_dir, "table1_data_description.csv"))
by_grp <- panel[O == 1, .(n = .N, mean_steps = round(mean(R)), median_steps = median(R)),
                by = .(group = ifelse(!avail, "unavailable", ifelse(A == 1, "sent", "available, not sent")))]
print(by_grp); fwrite(by_grp, file.path(out_dir, "table1b_steps_by_group.csv"))
say("panel built: ", nrow(panel), " decisions")

## Figure 1: the real reward distribution (why no distribution is assumed)
y <- panel$R[panel$O == 1]
p1 <- ggplot(data.frame(y = log10(y + 1)), aes(y)) +
  geom_histogram(bins = 60, fill = "grey60", colour = "white") +
  labs(x = "log10(steps + 1) in the 30 min after a decision", y = "Number of decisions",
       title = sprintf("Observed rewards: %.0f%% zeros, median %.0f, mean %.0f, max %.0f steps",
                       100 * mean(y == 0), median(y), mean(y), max(y))) + theme_bw()
figure("fig1_reward_distribution", p1)

## STAR-EM: choose the number of hidden day types K (held-out log-likelihood, MAR fit)
say("STAR-EM: choosing K (hidden day types) among ", paste(STAREM_KS, collapse = ","), " by held-out log-likelihood ...")
kt <- star_choose_K(panel, qc)
print(kt); fwrite(kt, file.path(out_dir, "table3e_STAR-EM_choose_K.csv"))
STAREM_K <- kt$K[which.max(kt$heldout_loglik)]
say("STAR-EM: chosen K = ", STAREM_K, " (used for the real data and the simulation)")
say("true policy values with ", N_TRUTH, " resampled participants per policy ...")
pool <- make_pool(panel)
make_truth <- function(hidden) {
  truth <- true_values(pool, N = N_TRUTH, hidden = hidden)
  tv <- setNames(truth$truth, truth$policy); tm <- setNames(truth$mcse, truth$policy)
  rbind(truth[, .(policy, truth, mcse)],
    data.table(policy = c("always-never", "tailored-always", "capped-always"),
               truth  = c(tv["always"] - tv["never"], tv["tailored"] - tv["always"], tv["capped"] - tv["always"]),
               mcse   = sqrt(c(tm["always"]^2 + tm["never"]^2, tm["tailored"]^2 + tm["always"]^2,
                               tm["capped"]^2 + tm["always"]^2))))[, world := if (hidden) "hidden" else "base"][]
}
truth_all <- rbind(make_truth(FALSE), make_truth(TRUE))
print(truth_all); fwrite(truth_all, file.path(out_dir, "table4_true_values.csv"))
say("fitting the missingness model on the real data ...")
mm <- fit_miss_model(panel)
## realism check: does the learned generator reproduce the real missingness structure
## (overall rate, spells, person-to-person differences) when run at the real rate?
set.seed(99)
dchk <- gen_cohort(200, pool)
dchk <- apply_missing(dchk, "learned", 1 - mean(panel$O), 0, 0, "linear", qc, mm)
mstat <- function(O, prevobs, id) c(missing_rate = 1 - mean(O),
  P_miss_after_miss = mean(O[prevobs == 0] == 0), P_miss_after_obs = mean(O[prevobs == 1] == 0),
  sd_of_person_rates = sd(tapply(1 - O, id, mean)))
chk <- data.table(statistic = names(mstat(panel$O, panel$prevobs, panel$id)),
                  real = round(mstat(panel$O, panel$prevobs, panel$id), 3),
                  simulated = round(mstat(dchk$O, dchk$prevobs, dchk$sid), 3))
print(chk); fwrite(chk, file.path(out_dir, "table4b_missingness_realism.csv"))
saveRDS(list(panel = panel, qc = qc, pool = pool, truth_all = truth_all, mm = mm, STAREM_K = STAREM_K), prep_file)
fwrite(tasks, file.path(out_dir, "tasks.csv"))
say(sprintf("prep done: K = %d; %d array tasks (tasks 1-%d: N = 37 scenarios; %d-%d: N = %d)",
            STAREM_K, nrow(tasks), max(tasks[cond < 8]$task), min(tasks[cond == 8]$task), nrow(tasks), N_BIG))
}

## ---------------------------------------------------------------- real
if (MODE == "real") {
P <- readRDS(prep_file); panel <- P$panel; qc <- P$qc; STAREM_K <- P$STAREM_K
## real-data analysis (fine delta grid)
DS_REAL <- CJ(d0 = DELTAS_REAL, dA = DELTAA_REAL)
say("real-data analysis: ", nrow(DS_REAL), " (delta0, deltaA) pairs x 4 policies + bounds ...")
t0 <- Sys.time()
st_tab <- star_state_table(panel, qc)
print(st_tab); fwrite(st_tab, file.path(out_dir, "table3f_STAR-EM_hidden_states.csv"))
## the lectures' point, on the real data: EM versus direct numerical MLE of the same model
evm <- rbind(em_vs_mle(panel, qc, 0), em_vs_mle(panel, qc, -1))
print(evm); fwrite(evm, file.path(out_dir, "table3d_EM_vs_direct_MLE.csv"))
nd_real  <- estimate_all(panel, qc, DS = DS_REAL)
res_real <- rbind(infer(nd_real, B = NBOOT), attr(nd_real, "mi"), use.names = TRUE)
em_real <- attr(nd_real, "em")
print(em_real); fwrite(em_real, file.path(out_dir, "table3c_EM_diagnostics_profile_likelihood.csv"))
say("EM: max log-likelihood ranges ", round(diff(range(em_real$max_loglik)), 4),
    " across the delta grid (about 0 = flat: the data cannot identify delta)")
cost_real <- attr(nd_real, "timing")
print(cost_real); fwrite(cost_real, file.path(out_dir, "table7a_compute_cost_realdata.csv"))
fwrite(res_real, file.path(out_dir, "table2_realdata_estimates.csv"))
say("real-data analysis done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
print(res_real[delta == 0, .(estimator, policy, est = round(est, 1), lo_w = round(lo_w, 1), hi_w = round(hi_w, 1))])

## distribution-free bounds (no assumption on how R and O are related)
bnd <- dcast(res_real[startsWith(estimator, "Bound")], policy ~ estimator, value.var = c("est", "lo_w", "hi_w"))
bnd <- bnd[, .(policy, lower = round(`est_Bound-lower`, 1), upper = round(`est_Bound-upper`, 1),
               ci_lower = round(`lo_w_Bound-lower`, 1), ci_upper = round(`hi_w_Bound-upper`, 1))]
print(bnd); fwrite(bnd, file.path(out_dir, "table3b_worst_case_bounds.csv"))

res_real[, dA := fifelse(grepl("dA=", estimator), suppressWarnings(as.numeric(sub(".*dA=([-0-9.]+).*", "\\1", estimator))), 0)]
curve <- res_real[estimator %in% c("MAR-DR", "STAR-DR") & grepl("-", policy)]
tip <- curve[, {
  base <- .SD[delta == 0]; sig0 <- base$lo_w > 0 | base$hi_w < 0
  incl <- .SD[lo_w <= 0 & hi_w >= 0]
  .(est_MAR = round(base$est, 1), CI_MAR = sprintf("[%.1f, %.1f]", base$lo_w, base$hi_w),
    significant_at_MAR = sig0,
    tipping_delta = if (sig0 && nrow(incl)) incl$delta[which.min(abs(incl$delta))] else NA_real_)
}, by = policy]
print(tip); fwrite(tip, file.path(out_dir, "table3_tipping_points.csv"))
curve2 <- res_real[grepl("DR", estimator) & grepl("-", policy)]
curve2[, dA_lab := factor(dA)]
bl <- res_real[startsWith(estimator, "Bound") & grepl("-", policy), .(lw = min(est), up = max(est)), by = policy]
p2 <- ggplot(curve2, aes(delta, est, colour = dA_lab, fill = dA_lab)) +
  geom_rect(data = bl, aes(xmin = -Inf, xmax = Inf, ymin = lw, ymax = up), inherit.aes = FALSE,
            fill = "grey90", colour = NA) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lo_w, ymax = hi_w), alpha = 0.15, colour = NA) + geom_line() + geom_point() +
  facet_wrap(~ policy, scales = "free_y") +
  scale_colour_manual(values = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"), name = "deltaA (sent days)") +
  scale_fill_manual(values = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"), name = "deltaA (sent days)") +
  labs(x = "Assumed tilt delta0 (0 = missing at random)", y = "Difference in steps per day",
       title = "HeartSteps: policy contrasts by STAR-DR (95% wild cluster bootstrap)",
       subtitle = "Grey band = distribution-free worst-case bounds (no assumption on how R and O are related)") +
  theme_bw() + theme(legend.position = "bottom")
figure("fig2_realdata_tipping", p2, w = 14, h = 5.5)
curve3 <- res_real[grepl("^STAR-EM", estimator) & grepl("-", policy)]
curve3[, dA_lab := factor(dA)]
p2b <- ggplot(curve3, aes(delta, est, colour = dA_lab, fill = dA_lab)) +
  geom_rect(data = bl, aes(xmin = -Inf, xmax = Inf, ymin = lw, ymax = up), inherit.aes = FALSE,
            fill = "grey90", colour = NA) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lo_w, ymax = hi_w), alpha = 0.15, colour = NA) + geom_line() + geom_point() +
  facet_wrap(~ policy, scales = "free_y") +
  scale_colour_manual(values = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"), name = "deltaA (sent days)") +
  scale_fill_manual(values = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"), name = "deltaA (sent days)") +
  labs(x = "Assumed tilt delta0 (0 = missing at random)", y = "Difference in steps per day",
       title = "HeartSteps: policy contrasts by STAR-EM (95% pooled wild-bootstrap intervals, no distribution assumed)",
       subtitle = "Grey band = distribution-free worst-case bounds") +
  theme_bw() + theme(legend.position = "bottom")
figure("fig2b_realdata_STAR-EM", p2b, w = 14, h = 5.5)


say("real-data mode done in ", round(difftime(Sys.time(), T_START, units = "mins"), 1), " min")
}

## ---------------------------------------------------------------- sim
if (MODE == "sim") {
P <- readRDS(prep_file); panel <- P$panel; qc <- P$qc; pool <- P$pool; mm <- P$mm; STAREM_K <- P$STAREM_K
cost_rows <- list(); em_rows <- list()
if (PILOT) {
  ## timing pilot: 2 replications of the heaviest N = 37 scenario (3 tilts) and of the N = 150 scenario
  pt <- rbindlist(lapply(c(3L, 8L), function(cnd) rbindlist(lapply(1:2, function(r) {
    t0 <- proc.time()[["elapsed"]]; x <- run_rep(cnd, r)
    hwm <- tryCatch(as.numeric(gsub("[^0-9]", "", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) / 1024^2,
                    error = function(e) NA_real_)                      # peak resident memory so far (GB)
    say(sprintf("pilot: scenario %d rep %d done in %.0f s; peak memory so far %.2f GB", cnd, r,
                proc.time()[["elapsed"]] - t0, hwm))
    data.table(cond = cnd, rep = r, ok = !is.null(x), secs = proc.time()[["elapsed"]] - t0, peak_gb = hwm) }))))
  print(pt)
  s37 <- mean(pt[cond == 3]$secs); s150 <- mean(pt[cond == 8]$secs)
  ## the pilot runs on 2 cores; N = 37 tasks run on 1 core, so their time is inflated by 1.3
  ## (only the forests use 2 cores; the boosting is single-core). SU = cores x hours.
  su37  <- 7 * REPS * s37 * 1.3 / 3600
  su150 <- REPS_BIG * s150 * 2 / 3600
  txt <- c(sprintf("pilot (2 cores): %.0f s per replication, heaviest N = 37 scenario; %.0f s, N = %d scenario",
                   s37, s150, N_BIG),
           sprintf("projected SU: 7 scenarios x %d reps (1 core) = %.0f;  N = %d scenario x %d reps (2 cores) = %.0f;  total = %.0f",
                   REPS, su37, N_BIG, REPS_BIG, su150, su37 + su150),
           sprintf("N = %d scenario: about %.1f SU per replication -> reps that fit X SU = X / %.2f", N_BIG, s150 * 2 / 3600, s150 * 2 / 3600),
           sprintf("(the heaviest scenario is used for all 7, so the N = 37 figure is an upper estimate)"),
           sprintf("time limit check: a N = 37 task (%d reps) needs about %.1f h; a N = %d task (%d reps) about %.1f h",
                   CHUNK, CHUNK * s37 * 1.3 / 3600, N_BIG, CHUNK_BIG, CHUNK_BIG * s150 / 3600))
  txt <- c(txt, sprintf("peak memory: N = 37 scenario %.2f GB; N = %d scenario %.2f GB (Anvil shared gives about 1.9 GB per core)",
                        max(pt[cond == 3]$peak_gb), N_BIG, max(pt$peak_gb)))
  writeLines(txt, file.path(out_dir, "pilot_timing.txt")); say(txt)
} else {
  tk <- tasks[task == TASK]
  if (!nrow(tk)) { say("task ", TASK, " is beyond the ", nrow(tasks), " tasks -- nothing to do"); quit(save = "no") }
  outf <- file.path(out_dir, sprintf("sim_%03d.rds", TASK))
  if (file.exists(outf)) { say("already done: ", outf); quit(save = "no") }
  say(sprintf("task %d: scenario %d (%s), replications %d-%d", TASK, tk$cond, grid[cond == tk$cond]$scenario, tk$r1, tk$r2))
  sim <- lapply(tk$r1:tk$r2, function(r) { x <- run_rep(tk$cond, r)
    if (!is.null(x)) say(sprintf("  rep %d done (%.0f s)", r, x$secs[1])); x })
  saveRDS(list(sim = rbindlist(sim), cost = rbindlist(cost_rows), em = rbindlist(em_rows, fill = TRUE)), outf)
  say("saved ", outf, " in ", round(difftime(Sys.time(), T_START, units = "mins"), 1), " min")
}
}

## ---------------------------------------------------------------- combine
if (MODE == "combine") {
P <- readRDS(prep_file); truth_all <- P$truth_all
fl <- list.files(out_dir, "^sim_[0-9]+[.]rds$", full.names = TRUE)
say("combining ", length(fl), " finished tasks of ", nrow(tasks))
R <- lapply(fl, readRDS)
sim <- rbindlist(lapply(R, `[[`, "sim")); cost_rows <- lapply(R, `[[`, "cost"); em_rows <- lapply(R, `[[`, "em")
done <- sim[, .(reps_done = uniqueN(rep)), by = cond][grid[, .(cond, scenario)], on = "cond"]
done[is.na(reps_done), reps_done := 0L]; print(done); fwrite(done, file.path(out_dir, "table0_replications_done.csv"))
NREP <- max(done$reps_done)
fwrite(sim, file.path(out_dir, "sim_raw.csv"))
## paired comparison with the complete-data oracle in the SAME replication (removes the
## Monte Carlo noise shared by all estimators): mean(estimate - oracle) and its MC SE
orc <- sim[estimator == "Oracle-PDIS", .(cond, rep, policy, orc = est)]
gap <- merge(sim[!startsWith(estimator, "Bound")], orc, by = c("cond", "rep", "policy"))[,
  .(gap_vs_oracle = mean(est - orc), gap_mcse = sd(est - orc) / sqrt(.N)), by = .(cond, estimator, delta, policy)]
## computational cost per estimate (one simulated trial of N = 37 participants)
cost <- rbindlist(cost_rows)[, .(mean_seconds = round(mean(secs), 2), max_seconds = round(max(secs), 2),
                                 peak_memory_mb = round(max(peak_mb), 1)), by = component]
print(cost); fwrite(cost, file.path(out_dir, "table7_compute_cost_sim.csv"))
emd <- rbindlist(em_rows)[, .(fits = .N,
                              cells_mean_iterations = round(mean(iterations), 1), cells_pct_monotone = 100 * mean(monotone_ascent),
                              starem_pct_converged = 100 * mean(starem_converged),
                              starem_mean_iterations = round(mean(starem_iterations), 1),
                              starem_max_iterations = max(starem_iterations),
                              starem_pct_increasing_steps = round(mean(starem_pct_increasing), 1),
                              starem_mean_start_spread = round(mean(starem_start_spread), 3)),
                          by = .(cond, d0, dA, starem_K)]
emd <- merge(grid[, .(cond, scenario)], emd, by = "cond")
print(emd); fwrite(emd, file.path(out_dir, "table8_EM_convergence_sim.csv"))

sim <- merge(merge(sim, grid, by = "cond"), truth_all[, .(world, policy, truth)], by = c("world", "policy"))
perf <- sim[, {
  n <- .N; err <- est - truth
  cw <- mean(lo_w <= truth & truth <= hi_w); cb <- mean(lo <= truth & truth <= hi)
  .(nsim = n, bias = mean(err), bias_mcse = sd(est) / sqrt(n), emp_se = sd(est),
    rmse = sqrt(mean(err^2)), cover_wild = cw, cover_boot = cb,     # cover_boot: pairs bootstrap (DR) / Rubin t (MI, comparison)
    cover_mcse = sqrt(cw * (1 - cw) / n), width_wild = mean(hi_w - lo_w))
}, by = .(cond, scenario, rate, d0_true, dA_true, shape, estimator, delta, policy)]
fwrite(perf[!startsWith(estimator, "Bound")], file.path(out_dir, "table5_sim_performance.csv"))

## bounds: does [lower, upper] contain the truth, and how wide is it?
## (the bounds are themselves ESTIMATED from N participants, so their point values can miss the
##  truth by sampling error; the confidence version uses lower-bound lo and upper-bound hi)
bsim <- dcast(sim[startsWith(estimator, "Bound")], cond + rep + policy + truth ~ estimator,
              value.var = c("est", "lo_w", "hi_w"))
bsim <- bsim[, .(nsim = .N,
                 cover_bounds_point = mean(`est_Bound-lower` <= truth & truth <= `est_Bound-upper`),
                 cover_bounds_95ci  = mean(`lo_w_Bound-lower` <= truth & truth <= `hi_w_Bound-upper`),
                 mean_width = mean(`est_Bound-upper` - `est_Bound-lower`)), by = .(cond, policy)]
bsim <- merge(grid[, .(cond, scenario)], bsim, by = "cond")
print(bsim); fwrite(bsim, file.path(out_dir, "table6_sim_bounds.csv"))

pc <- perf[policy == "always-never" & !startsWith(estimator, "Bound")]
pc[, Estimator := ifelse(grepl("^STAR|tilt", estimator), paste0(estimator, " (delta0=", delta, ")"), estimator)]
pc[, Condition := factor(scenario, levels = grid$scenario)]
p3 <- ggplot(pc, aes(bias, Estimator)) + geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(aes(x = bias - 2 * bias_mcse, xend = bias + 2 * bias_mcse, yend = Estimator), na.rm = TRUE) +
  geom_point() + facet_wrap(~ Condition, ncol = 3, scales = "free_x") +
  labs(x = "Bias in the always-vs-never contrast (steps per day), with 2 MC SE",
       y = NULL, title = sprintf("Simulation (up to %d replications per condition)", NREP)) + theme_bw()
figure("fig3_sim_bias", p3, w = 15, h = 13)
p4 <- ggplot(pc, aes(cover_wild, Estimator)) + geom_vline(xintercept = 0.95, linetype = 2) +
  geom_point() + facet_wrap(~ Condition, ncol = 3) + xlim(0, 1) +
  labs(x = "Coverage of the 95% wild cluster bootstrap interval", y = NULL) + theme_bw()
figure("fig4_sim_coverage", p4, w = 15, h = 13)


perf2 <- merge(perf[!startsWith(estimator, "Bound")], gap, by = c("cond", "estimator", "delta", "policy"), all.x = TRUE)
fwrite(perf2, file.path(out_dir, "table5_sim_performance.csv"))
say("combine done in ", round(difftime(Sys.time(), T_START, units = "mins"), 1), " min")
}
