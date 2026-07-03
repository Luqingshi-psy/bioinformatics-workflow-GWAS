
library(ieugwasr)
library(TwoSampleMR)
library(dplyr)
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyOTE3Mjk1NDA0QHFxLmNvbSIsImlhdCI6MTc3ODgzOTYwMSwiZXhwIjoxNzgwMDQ5MjAxfQ.tFR-_bUCtI_d4f3XVcj873cP5ngYIRTuOPo-Xr0HSEvme2X_3CGuBhNq3LGZZAWgGXD9wVphqq6gMVb0SdrOxev3JgygK7qyPko7eRXWKJvZVC1JYE6-dQFhyDW5kWVG2AFHMmzvrBXhXUt10iiTVC0hDwRvtZ42kyk5BHTpz35h3cRg8dKtBPiyB3XoOaqmvm9x8a2kT3XtJUi3TMCeY8--yD2yF6vlFEzfKUTbm1ijmb_GjJBB0SVbSP0nhD1pFo6eyn5_J3QzLlodcb_GVZ_aP0Bq13irdvXjbOkFaJSTZD41yZPugaoTuDHTJTfAwLmGHxWWDhDQZ8KKhBwMmg")

ao <- available_outcomes()


# 1. fetch GWAS list（id and trait only）
#    Note：first run downloads all data，please be patient
all_gwas <- gwasinfo()
trait_df <- all_gwas %>% select(id, trait)

keywords <- c(
  "Indoxyl sulfate", "3-indoxylsulfate", 
  "Propionic acid", "Propanoate", "Propionate",
  "Trimethylamine N-oxide", "TMAO",
  "2-hydroxybutyrate", "alpha-hydroxybutyrate"
)

results <- list()

for (kw in keywords) {
  matched <- trait_df[grepl(kw, trait_df$trait, ignore.case = TRUE), ]
  if (nrow(matched) > 0) {
    matched$keyword <- kw
    results[[kw]] <- matched
  } else {
    results[[kw]] <- data.frame(id = NA, trait = NA, keyword = kw)
  }
}

final <- bind_rows(results) %>% select(keyword, trait, id)

print(final)

write.csv(final, "metabolite_ieugwas_ids.csv", row.names = FALSE)


