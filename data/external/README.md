# Norwegian population mortality

`ssb_07902_2025.csv` contains the public 2025 period life-table probabilities
of death between exact ages x and x+1, by sex, from Statistics Norway table
[07902](https://www.ssb.no/en/statbank/table/07902). Retrieved 23 September 2026
from the full **Life tables** table on the official
[Deaths statistics page](https://www.ssb.no/en/befolkning/fodte-og-dode/statistikk/dode).
Table 07902 was updated 12 March 2026. The published probabilities per 1,000
were divided by 1,000 to obtain `qx`; ages 0 through 106 and both sexes are retained.
The source's three-decimal precision per 1,000 is preserved. These are public
population statistics, not trial data.

The extrapolation report uses piecewise constant annual hazards
`-log(1-qx)` at attained age, holding 2025 mortality rates fixed through follow-up.
It does not forecast mortality improvements or change economic model inputs.
