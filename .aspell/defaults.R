## see https://dirk.eddelbuettel.com/blog/2017/08/10/#008_aspell_cran_incoming;
## https://github.com/r-lib/usethis/issues/1466
Rd_files <- vignettes <- R_files <- description <-
    list(encoding = "UTF-8",
         language = "en",
         dictionaries = c("en_stats", "mbest"))
## update word list
if (FALSE) {
    saveRDS(c("Schmaus", "Zhang"), file = "mbest.rds")
}

