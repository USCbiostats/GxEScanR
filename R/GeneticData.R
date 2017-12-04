FamilyFileInfo <- function(familyFile) {
  df <- read.table(familyFile, stringsAsFactors = FALSE)
  if (ncol(df) != 6) {
    print("Family file does not have 6 columns")
    return (0);
  }
  return (nrow(df))
}