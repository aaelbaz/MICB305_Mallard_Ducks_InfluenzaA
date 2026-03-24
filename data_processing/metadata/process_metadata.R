library(tidyverse)

run_info <- read_csv("duck_flu_sra_metadata.csv")
run_info <- rename(run_info, SampleID = Identifying_Code)
metadata <- read_csv("duck_flu_paper_metadata.csv")

merged <- merge(run_info, metadata, by=c("SampleID"))

# Check that data merged correctly (values expected to match actually match)
bool <- FALSE
for (row in 1:nrow(merged)){
  if (merged$Influenza.x[row] != merged$Influenza.y[row]){
    bool <- TRUE
  }
  if (merged$SampleID[row] != merged$`Library Name`[row]){
    bool <- TRUE
  }
}
bool
# Still false, so should be good.

colnames(merged)

# Wrangle metadata
metadata_trimmed_conservative <- merged |> select(SampleID, Run, Bases, BioSample, Bytes,
                                                  Collection_Date, Experiment, geo_loc_name_country_continent,
                                                  geo_loc_name, HASubType_NASubType, HOST, Influenza.x,
                                                  isolation_source, lat_lon, LibraryLayout,
                                                  LibrarySelection, LibrarySource, Organism, `SRA Study`,
                                                  sex, IAV_level, `Ct value`, Genotype, HASubType, NASubType)

colnames(metadata_trimmed_conservative)

metadata_trimmed <- metadata_trimmed_conservative |> select(-geo_loc_name_country_continent,
                                                            -LibraryLayout, -Bytes, -Organism,
                                                            -LibrarySource, -BioSample,
                                                            -Experiment, -`SRA Study`,
                                                            -LibrarySelection, -Organism, -lat_lon,
                                                            -HOST, -isolation_source)


# Further Wrangling into tidy form that can be used in QIIME2 Pipelines
metadata_trimmed <- metadata_trimmed |> rename(sample_name = Run,
                                               sample_identifier = SampleID,
                                               num_bases = Bases,
                                               collection_date = Collection_Date,
                                               influenza_status = Influenza.x,
                                               influenza_level = IAV_level,
                                               ct_value = `Ct value`,
                                               genotype = Genotype,
                                               ha_subtype = HASubType,
                                               na_subtype = NASubType,
                                               HA_NA = HASubType_NASubType)

metadata_trimmed <- separate_wider_delim(metadata_trimmed, geo_loc_name, ":", names = c("country", "specific"))
metadata_trimmed <- metadata_trimmed |> mutate(specific = replace(specific, specific == " California\\, Grizzly Island Wildlife Area SolanoCounty", " California\\, Solano County\\, Grizzly Island Wildlife Area\\, Suisun Marsh"))
## ^^ I looked it up and Grizzly Island Wildlife Area is within Suisun Marsh (https://wildlife.ca.gov/Lands/Places-to-Visit/Grizzly-Island-WA)
metadata_trimmed <- separate_wider_delim(metadata_trimmed, specific, "\\,", names = c("state", "county", "precise_location", "region"))

metadata_trimmed <- metadata_trimmed |> mutate(country = str_squish(country))
metadata_trimmed <- metadata_trimmed |> mutate(state = str_squish(state))
metadata_trimmed <- metadata_trimmed |> mutate(region = str_squish(region))
metadata_trimmed <- metadata_trimmed |> mutate(precise_location = str_squish(precise_location))

## Reordering Columns to work with QIIME2
metadata_trimmed <- relocate(metadata_trimmed, sample_name, .before = sample_identifier)
metadata_trimmed <- relocate(metadata_trimmed, HA_NA, .after = na_subtype)


# Export & Re-Import to Check Everything is as Expected
write_tsv(metadata_trimmed, "duck_flu_metadata_wrangled.tsv")

wrangled_metadata <- read_tsv("duck_flu_metadata_wrangled.tsv")

