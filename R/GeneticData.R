FamilyFileCheck <- function(familyFile) {
  df <- read.table(familyFile, stringsAsFactors = FALSE)
  if (ncol(df) != 6) {
    print("Family file does not have 6 columns")
    return (0);
  }
  return (nrow(df))
}

MapFileCheck <- function(mapFile) {
  df <- read.table(mapFile, stringsAsFactors = FALSE)
  if (ncol(df) != 6) {
    print("Map file does not have 6 columns")
    return (0);
  }
  if (is.numeric(df[,3]) == FALSE | is.numeric(df[,4]) == FALSE) {
    print("Third or fourth column of map file is not numeric")
    return (0)
  }
  return (nrow(df))
}