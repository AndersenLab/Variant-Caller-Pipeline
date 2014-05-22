# Checking Strains

f <- group_by(concordance_results[[1]], Query) %.%
  filter(max(Concordance)==Concordance)