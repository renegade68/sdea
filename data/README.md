# Example data and redistribution

## Synthetic benchmark

The 9,040 firms are generated from seed 119. The four reference scores per firm in `tests/fixtures/full_dea_9040.csv` and the measurements in `tests/validation_9040.csv` are project-generated results. These files can accompany a public release without permission from an empirical data provider. The project can specify separate reuse terms for its own contributions.

## Swedish pharmacies

Daniel Santín deposited version 1 of the replication package for the study by Juan Aparicio and Daniel Santín on 15 August 2025. The deposit includes `pharmacy.xlsx` under CC BY 4.0 and credits the original study by Färe, Grosskopf, Lindgren and Roos (1992). The license permits redistribution and adaptation, including commercial use. A public Stata extract should credit the deposit and original study, link to the dataset and license, retain supplied notices, and describe the conversion. [Versioned deposit](https://doi.org/10.17632/dx95hwwcwz.1) · [Original study](https://doi.org/10.1007/BF00158770) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

We compared the deposited workbook (SHA-256 `ebfe070d874974755fddba7f83f2499b6d05a1138fa18473a7a61855f6071104`) with the archived Stata file identified below. We assigned the year from sheets R80–R89, used the explicit ID for 1980 and the 1-based pharmacy row within each later sheet, and compared the eight named inputs and outputs. All 3,360 values across the 420 pharmacy-years matched exactly. The archived file also contains `tl`. All eleven archived columns match the project’s earlier pharmacy CSV, but the licensed workbook contains no `tl` field and we found no redistribution terms for that additional column. The empirical example needs only the ten verified columns (year, pharmacy ID and eight measures). A Stata extract built from those workbook columns can be distributed with the attribution above and a note that the annual sheets were stacked and the identifiers added or renamed.

## U.S. banks

The FDIC catalog marks the official FFIEC Call Report bulk distributions, including the four-period balance sheet and income statement files, as U.S. public domain and unrestricted. Its distribution record expressly permits copying, publication and adaptation. A panel derived solely from these public fields can be published with the source, reporting periods and transformation details. [FDIC distribution and rights record](https://catalog.data.gov/dataset/ffiec-call-reports) · [FFIEC bulk download](https://cdr.ffiec.gov/public/PWS/DownloadBulkData.aspx)

The archived Stata file has 227,488 distinct bank-quarter keys and 17 fields. The source series consists of annual FFIEC bulk files covering four reporting quarters each. In a check against the official series, the March 2005 Bank of Hancock County row agrees with FDIC BankFind on employees (25), assets ($79.758 million) and deposits ($65.568 million). The full mapping of source fields, units, derived columns and sample filters has yet to be reproduced for this copy. Retain it externally until that mapping establishes the provenance of all fields, or publish a fresh extract from the official distribution with its preparation script. [FDIC BankFind observation](https://api.fdic.gov/banks/financials?filters=CERT%3A10057%20AND%20REPDTE%3A20050331&fields=CERT%2CNAME%2CASSET%2CDEP%2CNUMEMP%2CREPDTE&format=json)

## Running the empirical examples

After following the installation steps in the [package README](../README.md#installation), run the command below from the sdea folder. Replace the placeholder with the directory containing both Stata files. The empty second argument selects the core installed alongside sdea. The examples read the data in place.

```stata
do examples/empirical.do "." "" "<data-directory>"
```

| File | Coverage and source | Variables used in the example |
| --- | --- | --- |
| `pharmacy_yearly_1980_1989.dta` | Swedish pharmacy panel, 42 pharmacy groups observed annually for 10 years (420 rows). | `dmu` identifies the group and `year` the year. The deposit defines `x1` and `x2` as annual pharmacist and technical-staff hours, `x3` as equipment services and `x4` as building services. It describes `y1`–`y3` as counts of hospital drug deliveries, outpatient prescriptions and medical appliances, and `y4` as over-the-counter goods at 1980 prices. |
| `us_banks_quarterly_2005_2012.dta` | U.S. bank panel derived from FFIEC Call Reports, 2005 Q1 to 2012 Q4 (227,488 bank-quarter observations). | `id` identifies the bank (IDRSSD). `year` and `quarter` select the reporting period. The example uses deposits and full-time equivalent employees as inputs, and assets and loans as outputs. |

The examples compare the 1989 pharmacy cross-section and the first 200 eligible banks by ID in 2005 Q1 with the full-sample DEA solver. Eligibility requires positive, observed values for all four bank model variables. Reconstructing the processed bank panel from the raw Call Reports would also require its extraction, transformation and sample selection procedure.

## SHA-256 of processed files

```text
64d0286a8af4b676ba93c050926290e25b6b7b230c7b4538364151d52bed58d5  pharmacy_yearly_1980_1989.dta
ce8f599e6cf9d2b873beb04927422c84775b9886271db494f096d63d38e42b2e  us_banks_quarterly_2005_2012.dta
```
