# External data for the empirical examples

Supply the directory containing the two processed Stata files as the third argument to `examples/empirical.do`. The examples read the files in place. Redistribution terms for these processed copies have yet to be established.

| File | Coverage and source | Variables used in the example |
| --- | --- | --- |
| `pharmacy_yearly_1980_1989.dta` | Swedish pharmacy panel attributed to Färe et al. (1992), 42 pharmacy groups observed annually for 10 years (420 rows). | `dmu` identifies the group and `year` the year. Inputs `x1`–`x4` are pharmacist labor, technical-staff labor, building services and equipment services. Outputs `y1`–`y4` are hospital drug deliveries, outpatient prescriptions, medical appliances and over-the-counter goods, in thousands of Swedish kronor. |
| `us_banks_quarterly_2005_2012.dta` | U.S. bank panel derived from FFIEC Call Reports, 2005 Q1 to 2012 Q4 (227,488 bank-quarter observations). | `id` identifies the bank (IDRSSD). `year` and `quarter` select the reporting period. The example uses deposits and full-time equivalent employees as inputs, and assets and loans as outputs. |

The examples compare the 1989 pharmacy cross-section and the first 200 eligible banks by ID in 2005 Q1 with the full-sample DEA solver. Eligibility requires positive, observed values for all four bank model variables. Reconstructing the processed bank panel from the raw Call Reports would also require its extraction, transformation and sample selection procedure.

## SHA-256 of processed files

```text
64d0286a8af4b676ba93c050926290e25b6b7b230c7b4538364151d52bed58d5  pharmacy_yearly_1980_1989.dta
ce8f599e6cf9d2b873beb04927422c84775b9886271db494f096d63d38e42b2e  us_banks_quarterly_2005_2012.dta
```
