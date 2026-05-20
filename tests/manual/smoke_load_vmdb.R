#!/usr/bin/env Rscript
# Smoke-test load_vmdb_rates() / load_vmdb_magnitudes() against the live
# imo-vmdb server backed by data/vmdb.db. Direct SQLite queries on the
# same DB serve as the ground truth.
#
# Prerequisites:
#   - imo-vmdb server is running at http://127.0.0.1:8000
#   - vismeteor was installed from current source (devtools::install)
#   - Environment variable VMDB_DB_PATH points to the SQLite file
#     (most convenient: put it in a project-local .Renviron, which
#     RStudio loads automatically when the project is opened)

suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(vismeteor)
})

DB_PATH <- Sys.getenv("VMDB_DB_PATH", unset = NA_character_)
if (is.na(DB_PATH) || !nzchar(DB_PATH)) {
    stop(
        "VMDB_DB_PATH is not set. ",
        "Add it to ~/.Renviron or a project-local .Renviron, e.g.:\n",
        "  VMDB_DB_PATH=/path/to/imo-vmdb/data/vmdb.db\n",
        "or set it once in the R console: ",
        "Sys.setenv(VMDB_DB_PATH = \"...\")",
        call. = FALSE
    )
}
if (!file.exists(path.expand(DB_PATH))) {
    stop(
        "VMDB_DB_PATH does not point to an existing file: ", DB_PATH,
        call. = FALSE
    )
}
BASE_URL <- "http://127.0.0.1:8000/api/v1"

con <- dbConnect(SQLite(), path.expand(DB_PATH))

results <- list()
record <- function(name, ok, n_r = NA_integer_, n_sql = NA_integer_, note = "") {
    res <- data.frame(
        fall = name, ok = ok, n_r = n_r, n_sql = n_sql, anmerkung = note,
        stringsAsFactors = FALSE
    )
    results[[length(results) + 1L]] <<- res
    cat(sprintf(
        "[%s] %-32s  n_r=%-6s n_sql=%-6s  %s\n",
        if (ok) "OK " else "FAIL", name,
        format(n_r), format(n_sql), note
    ))
}

sql_count <- function(query) {
    as.integer(dbGetQuery(con, query)$n)
}

# Define one test case. Body is evaluated; if it errors, FAIL is recorded
# and the run continues. The body must end with a record(...) call on success.
# `endpoint` is "rates" or "magnitudes"; the recorded name is "<endpoint>/<name>".
case <- function(endpoint, name, body) {
    full_name <- paste0(endpoint, "/", name)
    tryCatch(body, error = function(e) {
        record(
            full_name, FALSE, NA_integer_, NA_integer_,
            sprintf("ERROR: %s", conditionMessage(e))
        )
    })
}

# ---- rates/shower_period: PER, 12.-13. August 2015 -----------------------
case("rates", "shower_period", {
    res <- load_vmdb_rates(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE shower = 'PER'
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    sample_id <- as.integer(as.character(res$observations$rate_id[1]))
    sql_row <- dbGetQuery(con, sprintf(
        "SELECT freq, lim_magn, t_eff, shower FROM rate WHERE id = %d", sample_id
    ))
    r_row <- res$observations[1, ]
    stopifnot(r_row$freq == sql_row$freq)
    stopifnot(abs(r_row$lim_magn - sql_row$lim_magn) < 1e-9)
    stopifnot(abs(r_row$t_eff - sql_row$t_eff) < 1e-9)
    stopifnot(as.character(r_row$shower) == sql_row$shower)
    record("rates/shower_period", TRUE, n_r, n_sql, "shower=PER, period")
})

# ---- rates/sporadic_mapping: NA -> 'SPO' -> shower IS NULL ---------------
case("rates", "sporadic_mapping", {
    res <- load_vmdb_rates(BASE_URL,
        shower = NA,
        period = c("2015-08-12", "2015-08-13")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE shower IS NULL
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(!any(is.na(res$observations$shower)))
    stopifnot(all(as.character(res$observations$shower) == "SPO"))
    record("rates/sporadic_mapping", TRUE, n_r, n_sql, "sporadic mapping NA -> SPO")
})

# ---- rates/sl_filter -----------------------------------------------------
case("rates", "sl_filter", {
    res <- load_vmdb_rates(BASE_URL,
        sl = c(135, 145),
        period = c("2015-01-01", "2018-12-31")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE sl_start >= 135 AND sl_end <= 145
        AND period_start >= '2015-01-01T00:00:00'
        AND period_end   <= '2018-12-31T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(all(res$observations$sl_start >= 135))
    stopifnot(all(res$observations$sl_end <= 145))
    record("rates/sl_filter", TRUE, n_r, n_sql, "sl filter")
})

# ---- rates/lim_magn ------------------------------------------------------
case("rates", "lim_magn", {
    res <- load_vmdb_rates(BASE_URL,
        lim_magn = c(6.0, 6.5),
        shower = "PER",
        period = c("2015-08-01", "2015-08-31")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE shower = 'PER'
        AND lim_magn BETWEEN 6.0 AND 6.5
        AND period_start >= '2015-08-01T00:00:00'
        AND period_end   <= '2015-08-31T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(all(res$observations$lim_magn >= 6.0 - 1e-9))
    stopifnot(all(res$observations$lim_magn <= 6.5 + 1e-9))
    record("rates/lim_magn", TRUE, n_r, n_sql, "lim_magn range")
})

# ---- rates/rate_id_multi -------------------------------------------------
case("rates", "rate_id_multi", {
    sample_rate_ids <- dbGetQuery(
        con,
        "SELECT id FROM rate ORDER BY id LIMIT 3"
    )$id
    res <- load_vmdb_rates(BASE_URL, rate_id = sample_rate_ids)
    n_r <- nrow(res$observations)
    stopifnot(n_r == 3L)
    stopifnot(setequal(
        as.integer(as.character(res$observations$rate_id)),
        sample_rate_ids
    ))
    record("rates/rate_id_multi", TRUE, n_r, 3L, "rate_id multi")
})

# ---- rates/session_id ----------------------------------------------------
case("rates", "session_id", {
    sample_session_id <- dbGetQuery(
        con,
        "SELECT session_id FROM rate GROUP BY session_id
         ORDER BY count(*) DESC LIMIT 1"
    )$session_id
    n_sql <- sql_count(sprintf(
        "SELECT count(*) AS n FROM rate WHERE session_id = %d", sample_session_id
    ))
    res <- load_vmdb_rates(BASE_URL, session_id = sample_session_id)
    n_r <- nrow(res$observations)
    stopifnot(n_r == n_sql)
    stopifnot(all(as.integer(as.character(res$observations$session_id))
    == sample_session_id))
    record("rates/session_id", TRUE, n_r, n_sql, "session_id")
})

# ---- rates/sun_moon_alt: sun_alt_max + moon_alt_max ----------------------
case("rates", "sun_moon_alt", {
    res <- load_vmdb_rates(BASE_URL,
        sun_alt_max = -10, moon_alt_max = 0,
        shower = "PER",
        period = c("2015-08-01", "2015-08-31")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE shower = 'PER'
        AND sun_alt <= -10 AND moon_alt <= 0
        AND period_start >= '2015-08-01T00:00:00'
        AND period_end   <= '2015-08-31T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(all(res$observations$sun_alt <= -10 + 1e-9))
    stopifnot(all(res$observations$moon_alt <= 0 + 1e-9))
    record("rates/sun_moon_alt", TRUE, n_r, n_sql, "sun_alt_max + moon_alt_max")
})

# ---- rates/with_sessions -------------------------------------------------
case("rates", "with_sessions", {
    res <- load_vmdb_rates(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13"),
        with_sessions = TRUE
    )
    stopifnot(is.data.frame(res$sessions))
    obs_sess <- unique(as.integer(as.character(res$observations$session_id)))
    sess_ids <- as.integer(as.character(res$sessions$session_id))
    stopifnot(setequal(obs_sess, sess_ids))
    stopifnot(is.factor(res$sessions$country))
    record(
        "rates/with_sessions", TRUE, nrow(res$observations), length(sess_ids),
        "with_sessions: session set-equal"
    )
})

# ---- rates/with_magnitudes (xtabs) ---------------------------------------
case("rates", "with_magnitudes", {
    res <- load_vmdb_rates(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13"),
        with_magnitudes = TRUE
    )
    stopifnot(methods::is(res$magnitudes, "table"))
    mag_ids_r <- as.integer(rownames(res$magnitudes))
    obs_mag_ids <- unique(stats::na.omit(
        as.integer(as.character(res$observations$magn_id))
    ))
    stopifnot(setequal(obs_mag_ids, mag_ids_r))
    mid <- mag_ids_r[1]
    r_total <- sum(res$magnitudes[as.character(mid), ])
    sql_total <- dbGetQuery(con, sprintf(
        "SELECT sum(freq) AS s FROM magnitude_detail WHERE id = %d", mid
    ))$s
    stopifnot(abs(r_total - sql_total) < 1e-9)
    record(
        "rates/with_magnitudes", TRUE, length(mag_ids_r), length(obs_mag_ids),
        sprintf(
            "magnitude_details xtabs (sample mid=%d sum=%.1f)",
            mid, r_total
        )
    )
})

# ---- rates/magn_solo: magn_solo type + value match -----------------------
case("rates", "magn_solo", {
    res <- load_vmdb_rates(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13")
    )
    stopifnot(is.logical(res$observations$magn_solo))
    idx <- which(!is.na(res$observations$magn_id))[1]
    sample_id <- as.integer(as.character(res$observations$rate_id[idx]))
    r_solo <- res$observations$magn_solo[idx]
    sql_solo <- dbGetQuery(con, sprintf(
        "SELECT magn_solo FROM rate WHERE id = %d", sample_id
    ))$magn_solo
    stopifnot(as.logical(r_solo) == as.logical(sql_solo))
    record(
        "rates/magn_solo", TRUE, sum(!is.na(res$observations$magn_solo)), NA_integer_,
        sprintf("magn_solo logical, sample id=%d -> %s", sample_id, r_solo)
    )
})

# ---- rates/per_plus_sporadic: PER + sporadic kombiniert ------------------
case("rates", "per_plus_sporadic", {
    res <- load_vmdb_rates(BASE_URL,
        shower = c("PER", NA),
        period = c("2015-08-12", "2015-08-13")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE (shower = 'PER' OR shower IS NULL)
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    levels_present <- unique(as.character(res$observations$shower))
    stopifnot(setequal(levels_present, c("PER", "SPO")))
    stopifnot(!any(is.na(res$observations$shower)))
    record("rates/per_plus_sporadic", TRUE, n_r, n_sql, "rates PER + sporadic (SPO present)")
})

# ---- rates/empty: leeres Ergebnis ----------------------------------------
case("rates", "empty", {
    res <- load_vmdb_rates(BASE_URL,
        shower = "PER",
        period = c("2020-01-01", "2020-01-02")
    )
    stopifnot(is.data.frame(res$observations))
    stopifnot(nrow(res$observations) == 0)
    stopifnot(is.null(res$sessions))
    stopifnot(is.null(res$magnitudes))
    record("rates/empty", TRUE, 0L, 0L, "empty result, sessions+magnitudes NULL")
})

# ---- rates/multi_shower: Multi-Shower + Sporadic -------------------------
case("rates", "multi_shower", {
    res <- load_vmdb_rates(BASE_URL,
        shower = c("PER", "GEM", NA),
        period = c("2015-08-12", "2015-08-13")
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE (shower IN ('PER','GEM') OR shower IS NULL)
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    shower_levels <- unique(as.character(res$observations$shower))
    stopifnot(all(shower_levels %in% c("PER", "GEM", "SPO")))
    record(
        "rates/multi_shower", TRUE, n_r, n_sql,
        sprintf(
            "multi-shower levels: %s",
            paste(sort(shower_levels), collapse = ",")
        )
    )
})

# ---- rates/period_types: period Typvarianten -----------------------------
case("rates", "period_types", {
    pwin_char <- c("2015-08-12", "2015-08-13")
    pwin_posix <- as.POSIXct(c("2015-08-12 00:00:00", "2015-08-13 23:59:59"),
        tz = "UTC"
    )
    pwin_date <- as.Date(c("2015-08-12", "2015-08-13"))
    r_char <- load_vmdb_rates(BASE_URL, shower = "PER", period = pwin_char)
    r_posix <- load_vmdb_rates(BASE_URL, shower = "PER", period = pwin_posix)
    r_date <- load_vmdb_rates(BASE_URL, shower = "PER", period = pwin_date)
    nrows <- c(
        nrow(r_char$observations),
        nrow(r_posix$observations),
        nrow(r_date$observations)
    )
    stopifnot(length(unique(nrows)) == 1L)
    record(
        "rates/period_types", TRUE, nrows[1], nrows[1],
        sprintf("period char/POSIXct/Date all -> %d rows", nrows[1])
    )
})

# ---- rates/large_page: große Page mit beiden Includes --------------------
case("rates", "large_page", {
    t0 <- Sys.time()
    res <- load_vmdb_rates(BASE_URL,
        period = c("2015-08-01", "2015-08-31"),
        with_sessions = TRUE,
        with_magnitudes = TRUE
    )
    dt <- as.numeric(Sys.time() - t0, units = "secs")
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM rate
      WHERE period_start >= '2015-08-01T00:00:00'
        AND period_end   <= '2015-08-31T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(n_r > 1000L)
    stopifnot(is.data.frame(res$sessions))
    stopifnot(methods::is(res$magnitudes, "table"))
    record(
        "rates/large_page", TRUE, n_r, n_sql,
        sprintf("large page + both includes (%.2fs)", dt)
    )
})

# ---- magnitudes/shower_period --------------------------------------------
case("magnitudes", "shower_period", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13"),
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM magnitude
      WHERE shower = 'PER'
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(all(as.character(res$observations$shower) == "PER"))
    record("magnitudes/shower_period", TRUE, n_r, n_sql, "magnitudes shower+period")
})

# ---- magnitudes/magn_id_multi --------------------------------------------
case("magnitudes", "magn_id_multi", {
    sample_magn_ids <- dbGetQuery(
        con,
        "SELECT id FROM magnitude ORDER BY id LIMIT 2"
    )$id
    res <- load_vmdb_magnitudes(BASE_URL,
        magn_id = sample_magn_ids,
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    stopifnot(n_r == 2L)
    stopifnot(setequal(
        as.integer(as.character(res$observations$magn_id)),
        sample_magn_ids
    ))
    record("magnitudes/magn_id_multi", TRUE, n_r, 2L, "magn_id multi")
})

# ---- magnitudes/with_sessions --------------------------------------------
case("magnitudes", "with_sessions", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13"),
        with_sessions = TRUE,
        with_magnitudes = FALSE
    )
    stopifnot(is.data.frame(res$sessions))
    obs_sess <- unique(as.integer(as.character(res$observations$session_id)))
    sess_ids <- as.integer(as.character(res$sessions$session_id))
    stopifnot(setequal(obs_sess, sess_ids))
    record(
        "magnitudes/with_sessions", TRUE, nrow(res$observations), length(sess_ids),
        "magnitudes + sessions set-equal"
    )
})

# ---- magnitudes/without_histogram: with_magnitudes = FALSE ---------------
case("magnitudes", "without_histogram", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13"),
        with_magnitudes = FALSE
    )
    stopifnot(nrow(res$observations) > 0)
    stopifnot(is.null(res$magnitudes))
    record(
        "magnitudes/without_histogram", TRUE, nrow(res$observations), NA_integer_,
        "with_magnitudes=FALSE -> magnitudes NULL"
    )
})

# ---- magnitudes/with_histogram: with_magnitudes = TRUE (default) ---------
case("magnitudes", "with_histogram", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = "PER",
        period = c("2015-08-12", "2015-08-13")
    )
    stopifnot(methods::is(res$magnitudes, "table"))
    mag_ids_r <- as.integer(rownames(res$magnitudes))
    obs_mag_ids <- as.integer(as.character(res$observations$magn_id))
    stopifnot(setequal(mag_ids_r, obs_mag_ids))
    mid <- mag_ids_r[1]
    r_total <- sum(res$magnitudes[as.character(mid), ])
    sql_total <- dbGetQuery(con, sprintf(
        "SELECT freq FROM magnitude WHERE id = %d", mid
    ))$freq
    stopifnot(abs(r_total - sql_total) < 1e-9)
    record(
        "magnitudes/with_histogram", TRUE, length(mag_ids_r), length(obs_mag_ids),
        sprintf(
            "xtabs sums match magnitude.freq (mid=%d, %.1f)",
            mid, r_total
        )
    )
})

# ---- magnitudes/per_plus_sporadic ----------------------------------------
case("magnitudes", "per_plus_sporadic", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = c("PER", NA),
        period = c("2015-08-12", "2015-08-13"),
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM magnitude
      WHERE (shower = 'PER' OR shower IS NULL)
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    levels_present <- unique(as.character(res$observations$shower))
    stopifnot(setequal(levels_present, c("PER", "SPO")))
    stopifnot(!any(is.na(res$observations$shower)))
    record("magnitudes/per_plus_sporadic", TRUE, n_r, n_sql, "magnitudes PER + sporadic (SPO present)")
})

# ---- magnitudes/sl_filter ------------------------------------------------
case("magnitudes", "sl_filter", {
    res <- load_vmdb_magnitudes(BASE_URL,
        sl = c(135, 145),
        period = c("2015-01-01", "2018-12-31"),
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM magnitude
      WHERE sl_start >= 135 AND sl_end <= 145
        AND period_start >= '2015-01-01T00:00:00'
        AND period_end   <= '2018-12-31T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(all(res$observations$sl_start >= 135))
    stopifnot(all(res$observations$sl_end <= 145))
    record("magnitudes/sl_filter", TRUE, n_r, n_sql, "magnitudes sl filter")
})

# ---- magnitudes/lim_magn -------------------------------------------------
case("magnitudes", "lim_magn", {
    res <- load_vmdb_magnitudes(BASE_URL,
        lim_magn = c(6.0, 6.5),
        shower = "PER",
        period = c("2015-08-01", "2015-08-31"),
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM magnitude
      WHERE shower = 'PER'
        AND lim_magn BETWEEN 6.0 AND 6.5
        AND period_start >= '2015-08-01T00:00:00'
        AND period_end   <= '2015-08-31T23:59:59'")
    stopifnot(n_r == n_sql)
    stopifnot(all(res$observations$lim_magn >= 6.0 - 1e-9))
    stopifnot(all(res$observations$lim_magn <= 6.5 + 1e-9))
    record("magnitudes/lim_magn", TRUE, n_r, n_sql, "magnitudes lim_magn range")
})

# ---- magnitudes/session_id -----------------------------------------------
case("magnitudes", "session_id", {
    sample_session_id <- dbGetQuery(
        con,
        "SELECT session_id FROM magnitude GROUP BY session_id
         ORDER BY count(*) DESC LIMIT 1"
    )$session_id
    n_sql <- sql_count(sprintf(
        "SELECT count(*) AS n FROM magnitude WHERE session_id = %d",
        sample_session_id
    ))
    res <- load_vmdb_magnitudes(BASE_URL,
        session_id = sample_session_id,
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    stopifnot(n_r == n_sql)
    stopifnot(all(as.integer(as.character(res$observations$session_id))
    == sample_session_id))
    record("magnitudes/session_id", TRUE, n_r, n_sql, "magnitudes session_id")
})

# ---- magnitudes/empty: leeres Ergebnis -----------------------------------
case("magnitudes", "empty", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = "PER",
        period = c("2020-01-01", "2020-01-02")
    )
    stopifnot(is.data.frame(res$observations))
    stopifnot(nrow(res$observations) == 0)
    stopifnot(is.null(res$sessions))
    stopifnot(is.null(res$magnitudes))
    record("magnitudes/empty", TRUE, 0L, 0L, "empty magnitudes")
})

# ---- magnitudes/multi_shower: Multi-Shower + Sporadic --------------------
case("magnitudes", "multi_shower", {
    res <- load_vmdb_magnitudes(BASE_URL,
        shower = c("PER", "GEM", NA),
        period = c("2015-08-12", "2015-08-13"),
        with_magnitudes = FALSE
    )
    n_r <- nrow(res$observations)
    n_sql <- sql_count("
      SELECT count(*) AS n FROM magnitude
      WHERE (shower IN ('PER','GEM') OR shower IS NULL)
        AND period_start >= '2015-08-12T00:00:00'
        AND period_end   <= '2015-08-13T23:59:59'")
    stopifnot(n_r == n_sql)
    record("magnitudes/multi_shower", TRUE, n_r, n_sql, "magnitudes multi-shower + sporadic")
})

# ---- Summary -------------------------------------------------------------
cat("\n=== Summary ===\n")
summary_df <- do.call(rbind, results)
print(summary_df, row.names = FALSE)
cat(sprintf(
    "\n%d/%d Fälle grün.\n",
    sum(summary_df$ok), nrow(summary_df)
))

dbDisconnect(con)
