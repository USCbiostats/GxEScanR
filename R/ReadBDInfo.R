ReadBDInfo <- function(bdInfo) {
  filename <- bdInfo$filename
  format <- integer(2)
  format[1] <- bdInfo$format
  format[2] <- bdInfo$version
  numSub <- integer(1)
  numSub <- bdInfo$NumSamples
  numSNPs <- integer(1)
  numSNPs <- bdInfo$NumSNPs
  bufferSize <- integer(1)
  bufferSize <- 100000001L
  buffer <- integer(bufferSize / 4)
  sections <- integer(1)
  filesize <- file.info(filename)$size
  locations <- numeric(numSNPs + 1)
  sections <- GetLocations(bdInfo$Indices, locations, filesize, bufferSize)
  snpSection <- integer(numSNPs)
  fileLocation <- numeric(sections + 1)
  snpLocation <- integer(numSNPs)
  GetSections(locations, snpSection, fileLocation, snpLocation, bufferSize)
  currectSection <- integer(1)
  currentSection <- -1L
  dosage <- numeric(numSub)
  p0 <- numeric(numSub)
  p1 <- numeric(numSub)
  p2 <- numeric(numSub)
  return (list(filename = filename,
               format = format,
               numSub = numSub,
               numSNPs = numSNPs,
               bufferSize = bufferSize,
               buffer = buffer,
               sections = sections,
               snpSection = snpSection,
               fileLocation = fileLocation,
               snpLocation = snpLocation,
               currentSection = currentSection,
               dosage = dosage,
               p0 = p0,
               p1 = p1,
               p2 = p2))
}
