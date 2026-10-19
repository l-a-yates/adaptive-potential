# Data_Productivity.csv

Productivity of female *Drosophila* from a paternal half-sib breeding design,
measured in six environments (3 temperatures x 2 diets). One row per female.

- 6906 females, 143 sires, 688 dams (full-sib families), 11 treatment-by-set cells
- Built from the per-block files in `raw/` by `raw/Analysis_CombineData.R`
  (productivity joined to development-vial data)

## Columns

| column | type | description |
|---|---|---|
| `Row` | int | Row ID in the combined source file (unique, 1-6906) |
| `TREAT` | chr | Treatment: temperature x diet, e.g. `T14DA` |
| `TEMP` | chr | Temperature: `T14`, `T25`, `T28` (°C) |
| `DIET` | chr | Diet: `DA` standard, `DB` low-calorie |
| `Box_Development` | chr | Box holding the development vial |
| `VIAL_Development` | int | Vial in which the female developed |
| `Offspring_Development` | int | Number of offspring emerging from that development vial (a measure of larval density) |
| `SIRE` | int | Sire ID; 138 of the 143 sires appear in all six treatments, the other 5 in one |
| `DAM` | chr | Dam `a`-`e` within sire; a dam is identified by `SIRE` + `DAM` |
| `REP` | int | Replicate female within a development vial (usually 1 or 2) |
| `SET` | int | Experimental set (block) within treatment, 1-3 |
| `TOTAL` | int | Productivity: total offspring produced over the laying period |

## Notes

- Laying periods differ by treatment (days): T14DA 11, T14DB 15, T25DA 6,
  T25DB 8, T28DA 6, T28DB 9. The analysis divides by these to get
  productivity per day.
- Zeros (no offspring) are common under stress: 33.6, 63.0, 3.6, 23.2, 3.1
  and 34.9 % for T14DA ... T28DB.
- Treatment-by-set cells: T14DA 1-2, T14DB 1-3, T25DA 1, T25DB 1-2,
  T28DA 1, T28DB 1-2.
