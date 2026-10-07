# TurkStat province teaching snapshot

Source: Turkish Statistical Institute (TurkStat / TÜİK), Geographical Statistics
Portal, https://cip.tuik.gov.tr/, retrieved 19 September 2026.
Four statistical indicators refer to 2025; raw JSON responses also retain the
other available years. Province GeoJSON is downloaded directly from the portal:
https://cip.tuik.gov.tr/assets/geometri/nuts3.json.
Its authoritative boundary vintage is unspecified. No basemap tiles are included.
No TRmaps or other prepared boundary/data package supplies these files.

Publicly provided TurkStat data may be reused with source attribution:
https://www.tuik.gov.tr/Kurumsal/Yasal_Uyari.
Attribute the data to TurkStat, separately from citations to gwrs and methods.
These source data remain attributed to TurkStat, not to the package author.

See dictionary.csv for units, indicator IDs and metadata URLs, and
snapshot-manifest.json for source URLs, retrieval times, sizes and SHA-256.
Raw JSON is unmodified. The CSV selects 2025 by label, parses decimal commas,
and joins on province codes 01-81. No observations are imputed or excluded.
The raw GeoJSON contains four invalid geometries; the vignette explicitly
repairs an in-memory projected copy and reports the geometry audit.

Run vignette("gwrs-turkey-provinces", package = "gwrs") for the complete example.
This small educational data set is not a big-data benchmark or a causal study.
