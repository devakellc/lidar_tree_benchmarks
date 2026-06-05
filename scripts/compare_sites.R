#!/usr/bin/env Rscript
# Cross-site structure-gradient comparison: pools each site's sweep and reports
# the density-sensitivity of detection by crown class for SJER (open oak) ->
# SOAP (mixed conifer) -> TEAK (red fir). Writes a combined CSV + figure.
d  <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
sites <- strsplit(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else
                  "SJER,SOAP,TEAK", ",")[[1]]
classes <- c("dominant","codominant","intermediate","suppressed")

pool_rung <- function(r, rl) {
  s <- r[r$rung == rl & r$chm_res == 0.5 & r$vwf_a == 0.10, ]
  if (!nrow(s)) return(NULL)
  s$tp_core <- round(s$precision * s$n_det)
  o <- data.frame(rung = rl, frdens = round(mean(s$frdens), 2),
                  n_ref = sum(s$n_ref),
                  recall = sum(s$TP)/sum(s$n_ref),
                  precision = sum(s$tp_core, na.rm=TRUE)/sum(s$n_det))
  o$F1 <- 2*o$recall*o$precision/(o$recall+o$precision)
  for (cl in classes) {
    nc <- sum(s[[paste0("n_",cl)]], na.rm=TRUE)
    tp <- sum(ifelse(s[[paste0("n_",cl)]]>0, round(s[[paste0("rec_",cl)]]*s[[paste0("n_",cl)]]),0), na.rm=TRUE)
    o[[paste0("rec_",cl)]] <- if (nc) tp/nc else NA
  }
  # overstory vs understory rollup
  no <- sum(s$n_dominant+s$n_codominant); nu <- sum(s$n_intermediate+s$n_suppressed)
  to <- sum(ifelse(s$n_dominant>0,round(s$rec_dominant*s$n_dominant),0),
            ifelse(s$n_codominant>0,round(s$rec_codominant*s$n_codominant),0), na.rm=TRUE)
  tu <- sum(ifelse(s$n_intermediate>0,round(s$rec_intermediate*s$n_intermediate),0),
            ifelse(s$n_suppressed>0,round(s$rec_suppressed*s$n_suppressed),0), na.rm=TRUE)
  o$rec_overstory  <- if (no) to/no else NA
  o$rec_understory <- if (nu) tu/nu else NA
  o$n_overstory <- no; o$n_understory <- nu
  o
}

all <- list()
for (st in sites) {
  f <- file.path(d, "neon", st, "sweep_results.csv")
  if (!file.exists(f)) { cat("skip", st, "(no results)\n"); next }
  r <- read.csv(f, stringsAsFactors = FALSE)
  rr <- do.call(rbind, lapply(c("native","8","4","2","1"), function(rl) pool_rung(r, rl)))
  rr <- cbind(site = st, rr)
  all[[st]] <- rr
  cat(sprintf("\n=== %s (n_ref/rung=%d; overstory n=%d understory n=%d) ===\n",
              st, rr$n_ref[1], rr$n_overstory[1], rr$n_understory[1]))
  print(rr[, c("rung","frdens","recall","precision","F1",
               "rec_overstory","rec_understory")], row.names=FALSE, digits=2)
}
comb <- do.call(rbind, all)
write.csv(comb, file.path(d, "neon", "cross_site_summary.csv"), row.names = FALSE)

## combined figure: overstory recall vs density, one line per site
fig <- file.path(d, "neon", "figs"); dir.create(fig, showWarnings = FALSE, recursive = TRUE)
png(file.path(fig, "structure_gradient.png"), 1050, 760, res = 130)
par(mar = c(4.2,4.2,2.5,1)); cols <- c(SJER="#1b9e77", SOAP="#d95f02", TEAK="#7570b3")
plot(NA, xlim = range(comb$frdens), ylim = c(0,1), log = "x",
     xlab = "first-return density (pulses/m^2)", ylab = "recall",
     main = "Overstory (solid) vs understory (dashed) recall by site")
for (st in names(all)) { s <- all[[st]]; o <- order(s$frdens)
  lines(s$frdens[o], s$rec_overstory[o],  type="o", pch=19, col=cols[st], lwd=2)
  uu <- s$rec_understory[o]; uu[s$n_understory[o] < 5] <- NA   # hide noisy understory (SJER n=2)
  lines(s$frdens[o], uu, type="o", pch=1, col=cols[st], lwd=2, lty=2) }
legend("topright", names(all), col=cols[names(all)], lwd=2, pch=19, bty="n")
dev.off()
cat(sprintf("\ncombined -> %s ; figure -> %s\n",
            file.path(d,"neon","cross_site_summary.csv"), fig))
