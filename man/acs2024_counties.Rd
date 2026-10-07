\name{acs2024_counties}
\alias{acs2024_counties}
\docType{data}
\title{County analysis data from the 2020--2024 ACS five-year Data Profiles}
\description{The 3,107-county analysis table used in \emph{Nonconvex Geographically Weighted Regression for Local Variable Selection}. The 91 columns retain the study's estimates, 90 percent margins of error, controlled-MOE flags, identifiers, projected coordinates and saved spatial folds.}
\usage{data(acs2024_counties)}
\format{A data frame with 3,107 rows and 91 columns:
\describe{
  \item{\code{fips}}{Five-character county FIPS identifier; leading zeros retained. .}
  \item{\code{county_name}}{Census county or county-equivalent name. .}
  \item{\code{state_fips}}{Two-character state FIPS identifier. .}
  \item{\code{state_abbreviation}}{Two-letter postal abbreviation. .}
  \item{\code{state_name}}{State name. .}
  \item{\code{county_fips}}{Three-character within-state county FIPS identifier. .}
  \item{\code{easting_m}}{Easting of a point on the county surface after projection to NAD83 / Conus Albers (EPSG:5070); not longitude or latitude. metres.}
  \item{\code{northing_m}}{Northing of a point on the county surface after projection to NAD83 / Conus Albers (EPSG:5070); not longitude or latitude. metres.}
  \item{\code{spatial_fold}}{Saved spatial fold from k-means on standardized projected coordinates; seed 20260907, nstart 100. Fold labels must be reused for the paper analysis. 1 to 5.}
  \item{\code{median_household_income_2024}}{Original county median household income estimate, before the natural logarithm. 2024 inflation-adjusted US dollars.}
  \item{\code{log_median_household_income_2024}}{Natural logarithm of median_household_income_2024. Estimate!!INCOME AND BENEFITS (IN 2024 INFLATION-ADJUSTED DOLLARS)!!Total households!!Median household income (dollars) natural log of dollar value.}
  \item{\code{bachelors_or_higher_pct}}{Percent!!EDUCATIONAL ATTAINMENT!!Population 25 years and over!!Bachelor's degree or higher percent.}
  \item{\code{disability_pct}}{Percent!!DISABILITY STATUS OF THE CIVILIAN NONINSTITUTIONALIZED POPULATION!!Total Civilian Noninstitutionalized Population!!With a disability percent.}
  \item{\code{foreign_born_pct}}{Percent!!PLACE OF BIRTH!!Total population!!Foreign-born percent.}
  \item{\code{limited_english_pct}}{Percent!!LANGUAGE SPOKEN AT HOME!!Population 5 years and over!!Language other than English!!Speak English less than "very well" percent.}
  \item{\code{broadband_subscription_pct}}{Percent!!COMPUTERS AND INTERNET USE!!Total households!!With a broadband Internet subscription percent.}
  \item{\code{under_18_pct}}{Percent!!SEX AND AGE!!Total population!!Under 18 years percent.}
  \item{\code{age_65_plus_pct}}{Percent!!SEX AND AGE!!Total population!!65 years and over percent.}
  \item{\code{hispanic_latino_pct}}{Percent!!HISPANIC OR LATINO AND RACE!!Total population!!Hispanic or Latino (of any race) percent.}
  \item{\code{black_non_hispanic_pct}}{Percent!!HISPANIC OR LATINO AND RACE!!Total population!!Not Hispanic or Latino!!Black or African American alone percent.}
  \item{\code{labor_force_participation_pct}}{Percent!!EMPLOYMENT STATUS!!Population 16 years and over!!In labor force percent.}
  \item{\code{unemployment_pct}}{Percent!!EMPLOYMENT STATUS!!Population 16 years and over!!In labor force!!Civilian labor force!!Unemployed percent.}
  \item{\code{management_business_science_arts_occupation_pct}}{Percent!!OCCUPATION!!Civilian employed population 16 years and over!!Management, business, science, and arts occupations percent.}
  \item{\code{service_occupation_pct}}{Percent!!OCCUPATION!!Civilian employed population 16 years and over!!Service occupations percent.}
  \item{\code{natural_resources_construction_maintenance_occupation_pct}}{Percent!!OCCUPATION!!Civilian employed population 16 years and over!!Natural resources, construction, and maintenance occupations percent.}
  \item{\code{production_transport_material_moving_occupation_pct}}{Percent!!OCCUPATION!!Civilian employed population 16 years and over!!Production, transportation, and material moving occupations percent.}
  \item{\code{agriculture_forestry_fishing_mining_industry_pct}}{Percent!!INDUSTRY!!Civilian employed population 16 years and over!!Agriculture, forestry, fishing and hunting, and mining percent.}
  \item{\code{manufacturing_industry_pct}}{Percent!!INDUSTRY!!Civilian employed population 16 years and over!!Manufacturing percent.}
  \item{\code{professional_scientific_management_admin_waste_industry_pct}}{Percent!!INDUSTRY!!Civilian employed population 16 years and over!!Professional, scientific, and management, and administrative and waste management services percent.}
  \item{\code{self_employed_nonincorporated_pct}}{Percent!!CLASS OF WORKER!!Civilian employed population 16 years and over!!Self-employed in own not incorporated business workers percent.}
  \item{\code{work_from_home_pct}}{Percent!!COMMUTING TO WORK!!Workers 16 years and over!!Worked from home percent.}
  \item{\code{mean_commute_minutes}}{Estimate!!COMMUTING TO WORK!!Workers 16 years and over!!Mean travel time to work (minutes) minutes.}
  \item{\code{vacant_housing_pct}}{Percent!!HOUSING OCCUPANCY!!Total housing units!!Vacant housing units percent.}
  \item{\code{mobile_home_pct}}{Percent!!UNITS IN STRUCTURE!!Total housing units!!Mobile home percent.}
  \item{\code{renter_occupied_pct}}{Percent!!HOUSING TENURE!!Occupied housing units!!Renter-occupied percent.}
  \item{\code{no_vehicle_pct}}{Percent!!VEHICLES AVAILABLE!!Occupied housing units!!No vehicles available percent.}
  \item{\code{severe_crowding_pct}}{Percent!!OCCUPANTS PER ROOM!!Occupied housing units!!1.51 or more percent.}
  \item{\code{log_median_household_income_2024__moe90}}{Census 90\% margin of error for log_median_household_income_2024. Despite the log-response prefix, this is the original income MOE in dollars; it has NOT been log-transformed. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. 2024 inflation-adjusted US dollars.}
  \item{\code{bachelors_or_higher_pct__moe90}}{Census 90\% margin of error for bachelors_or_higher_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{disability_pct__moe90}}{Census 90\% margin of error for disability_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{foreign_born_pct__moe90}}{Census 90\% margin of error for foreign_born_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{limited_english_pct__moe90}}{Census 90\% margin of error for limited_english_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{broadband_subscription_pct__moe90}}{Census 90\% margin of error for broadband_subscription_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{under_18_pct__moe90}}{Census 90\% margin of error for under_18_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{age_65_plus_pct__moe90}}{Census 90\% margin of error for age_65_plus_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{hispanic_latino_pct__moe90}}{Census 90\% margin of error for hispanic_latino_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{black_non_hispanic_pct__moe90}}{Census 90\% margin of error for black_non_hispanic_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{labor_force_participation_pct__moe90}}{Census 90\% margin of error for labor_force_participation_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{unemployment_pct__moe90}}{Census 90\% margin of error for unemployment_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{management_business_science_arts_occupation_pct__moe90}}{Census 90\% margin of error for management_business_science_arts_occupation_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{service_occupation_pct__moe90}}{Census 90\% margin of error for service_occupation_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{natural_resources_construction_maintenance_occupation_pct__moe90}}{Census 90\% margin of error for natural_resources_construction_maintenance_occupation_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{production_transport_material_moving_occupation_pct__moe90}}{Census 90\% margin of error for production_transport_material_moving_occupation_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{agriculture_forestry_fishing_mining_industry_pct__moe90}}{Census 90\% margin of error for agriculture_forestry_fishing_mining_industry_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{manufacturing_industry_pct__moe90}}{Census 90\% margin of error for manufacturing_industry_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{professional_scientific_management_admin_waste_industry_pct__moe90}}{Census 90\% margin of error for professional_scientific_management_admin_waste_industry_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{self_employed_nonincorporated_pct__moe90}}{Census 90\% margin of error for self_employed_nonincorporated_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{work_from_home_pct__moe90}}{Census 90\% margin of error for work_from_home_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{mean_commute_minutes__moe90}}{Census 90\% margin of error for mean_commute_minutes. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. minutes.}
  \item{\code{vacant_housing_pct__moe90}}{Census 90\% margin of error for vacant_housing_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{mobile_home_pct__moe90}}{Census 90\% margin of error for mobile_home_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{renter_occupied_pct__moe90}}{Census 90\% margin of error for renter_occupied_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{no_vehicle_pct__moe90}}{Census 90\% margin of error for no_vehicle_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{severe_crowding_pct__moe90}}{Census 90\% margin of error for severe_crowding_pct. Controlled MOE sentinel -555555555 is represented by zero with its companion flag. percentage points.}
  \item{\code{log_median_household_income_2024__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{bachelors_or_higher_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{disability_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{foreign_born_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{limited_english_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{broadband_subscription_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{under_18_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{age_65_plus_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{hispanic_latino_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{black_non_hispanic_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{labor_force_participation_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{unemployment_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{management_business_science_arts_occupation_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{service_occupation_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{natural_resources_construction_maintenance_occupation_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{production_transport_material_moving_occupation_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{agriculture_forestry_fishing_mining_industry_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{manufacturing_industry_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{professional_scientific_management_admin_waste_industry_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{self_employed_nonincorporated_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{work_from_home_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{mean_commute_minutes__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{vacant_housing_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{mobile_home_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{renter_occupied_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{no_vehicle_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
  \item{\code{severe_crowding_pct__moe90_controlled}}{TRUE when the corresponding Census MOE was a controlled-estimate sentinel (-555555555) and was stored as zero; FALSE otherwise. TRUE/FALSE.}
}}
\details{
The domain is the contiguous 48 states and the District of Columbia. Of 3,109 eligible geographic records, De Baca County, New Mexico (35011), lacks a median household income estimate, and Loving County, Texas (48301), lacks a mean commute estimate. These two records were excluded under the study's complete-case estimate policy. Alaska, Hawaii and territories are outside the domain. No estimate is imputed, winsorized or standardized in this data object.

The response is the natural logarithm of income in 2024 inflation-adjusted dollars. Predictors remain in their original units. Percentage estimates use the populations stated in the Census labels; they do not all share the same denominator. County associations do not identify individual-level or causal effects.

A Census controlled-estimate MOE sentinel (-555555555) is stored as zero and marked by a logical companion column. Other unavailable MOEs were treated as missing without automatic sample exclusion; this released table contains no missing values. In particular, \code{log_median_household_income_2024__moe90} retains the original dollar-scale income MOE. It is not an MOE on the log scale and must not be used as one.

The spatial folds use projected representative points, seed 20260907 and 100 k-means starts. Coordinates are metres in EPSG:5070 and are not geographic longitude/latitude. Polygon geometry and the common diagnostic graph are provided in the separate reproduction bundle, not in this data frame.

The machine-readable 91-column dictionary and provenance record are installed under \code{system.file("extdata", package = "gwrs")}. Preparation code and the public source snapshot are included in the repository's reproduction bundle.
}
\source{U.S. Census Bureau, 2020--2024 American Community Survey five-year Data Profiles, downloaded 7 September 2026: \url{https://api.census.gov/data/2024/acs/acs5/profile.html}. County and state cartographic boundaries: \url{https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html}.}
\references{
U.S. Census Bureau. 2020--2024 American Community Survey five-year Data Profiles.

Yuzbasi, B. (2026). \emph{Nonconvex Geographically Weighted Regression for Local Variable Selection}. Manuscript.

For work using this prepared dataset, please acknowledge the Census source, the accompanying paper and \code{gwrs}. Use \code{citation("gwrs")} for package citation information.
}
\examples{
data(acs2024_counties)
dim(acs2024_counties)
head(acs2024_counties[, c("fips", "county_name",
                         "median_household_income_2024", "spatial_fold")])
dictionary <- read.csv(system.file("extdata", "acs2024_counties_dictionary.csv",
                                   package = "gwrs"))
head(dictionary)
}
\keyword{datasets}
