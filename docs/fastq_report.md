# BGI1 Issues

A subset of the indices for __BGI1 RET1__ were mixed up somehow. The following table summarizes the problems.


| BGI1 Calls it | BGI1 Index | BGI2+3 Index  | Actual Strain | Comment            |
|---------------|------------|---------------|---------------|--------------------|
| AF16          | AATTCAAA   | AGGTACCA      | JU258         | Swapped with JU258 |
| JU258         | AGGTACCA   | AATTCAAA      | AF16          | Swapped with AF16  |
| DL238         | ACCAACTA   | AGTCAGAA      | CB4856_CGC    |                    |
| CB4856_CGC    | TGACGTCA   | ACCAACTA      | HK104         |                    |
| HK104         | AGTCAGAA   | TGACGTCA      | DL238         |                    |
| MY23          | TCGCAGGA   | CTGCGACA      | N2_CGC        |                    |
| N2_CGC        | CATCCGGA   | TCGCAGGA      | QX1430        |                    |
| QX1430        | CTGCGACA   | CATCCGGA      | MY23          |                    |

Some strains failed to sequence because of improperly specified indices.               

* WN2021
* CB4857_CGC

# BGI2 Issues

BGI mislabed libraries. __RET2__ and __RET3__ were sequenced in lane __0004__ and __0005__. However, the indices used in RET2 and RET3 are identicle to the sets used in RET6 and RET7. BGI therefore lanes __0004__ and __0005__ as RET6 and 7 even though they were truly RET2 and RET3.

For lane 1, the problem goes much deeper. The set of indices used in RET1 __are not__ identicle to __any__ other library, but if you take the union of index sets used in RET4 and RET5, you mostly cover the indices found in RET1 (although you don't capture all of them!). Apparently, this is what BGI did - meaning we were able to rescue a fair amount of the sample. The problem is that they demultiplexed using the indices 

| Lane	 | BGI Calls (RET) |	Actual (RET) |
| --- | --- | --- | --- |
| 0003	 | 4 + 5 |	1 |
| 0004	 | 6	| 2 |
| 0005	 | 7	| 3 |
| 0006	 | 4	| 4 |
| 0007	 | 5	| 5 |
| 0008	 | 6	| 6 |
| 0009	 | 7	| 7 |


| BGI Reported | BGI Called (RET) | Actual Library (RET) | Supposed Strain (Index Mixups) | Actual Strain (index Fixed) |
|--------------|------------------|----------------------|--------------------------------|-----------------------------|
| DL200        | 4                | 1                    | CB4852                         |                             |
| ED3073       | 4                | 1                    | JU778                          |                             |
| JU1400       | 4                | 1                    | JU360                          |                             |
| JU1896       | 4                | 1                    | ED3049                         |                             |
| LSJ1         | 4                | 1                    | JU1580                         |                             |
| NIC198       | 4                | 1                    | QG557                          |                             |
| NIC231       | 4                | 1                    | CX11315                        |                             |
| NIC4         | 4                | 1                    | NIC232                         |                             |
| QG557        | 4                | 1                    | CB4856 / CB4857                |                             |
| QX1212       | 4                | 1                    | CB4856 / CB4857                |                             |
| CB4855_CGC   | 5                | 1                    | MY23                           |                             |
| ED3052       | 5                | 1                    | AF16                           |                             |
| JU397        | 5                | 1                    | N2_CGC                         |                             |
| KR314        | 5                | 1                    | QX1430                         |                             |
| QG536        | 5                | 1                    | DL238                          |                             |
| QX2265       | 5                | 1                    | CB4856_CGC                     |                             |
| WN2011       | 5                | 1                    | JU258                          |                             |
| WN2018       | 5                | 1                    | HK104                          |                             |
| CB4853_CGC   | 6                | 2                    | JT11398                        |                             |
| CB4857_UK    | 6                | 2                    | QX1793                         |                             |
| CX11262      | 6                | 2                    | ED3017                         |                             |
| ED3012       | 6                | 2                    | CX11285                        |                             |
| EG4349       | 6                | 2                    | NIC2                           |                             |
| GXW1         | 6                | 2                    | JU363                          |                             |
| JU1212       | 6                | 2                    | AB4                            |                             |
| JU1395       | 6                | 2                    | JU782                          |                             |
| JU1440       | 6                | 2                    | EG4725                         |                             |
| JU1581       | 6                | 2                    | AB1                            |                             |
| JU2001       | 6                | 2                    | MY1                            |                             |
| JU393        | 6                | 2                    | CX11276                        |                             |
| JU406        | 6                | 2                    | QX1791                         |                             |
| JU792        | 6                | 2                    | PB306                          |                             |
| N2_HRH       | 6                | 2                    | ED3046                         |                             |
| NIC199       | 6                | 2                    | JU775                          |                             |
| NIC236       | 6                | 2                    | MY16                           |                             |
| PB303        | 6                | 2                    | JU346                          |                             |
| PX179        | 6                | 2                    | JU1200                         |                             |
| QG538        | 6                | 2                    | DL226                          |                             |
| QX1214       | 6                | 2                    | QG556                          |                             |
| QX2266       | 6                | 2                    | QX1792                         |                             |
| WN2001       | 6                | 2                    | NIC196                         |                             |
| WN2017       | 6                | 2                    | QX1794                         |                             |
| CB4854       | 7                | 3                    | EG4347                         |                             |
| CB4932       | 7                | 3                    | JU2007                         |                             |
| CX11271      | 7                | 3                    | NIC197                         |                             |
| ED3040       | 7                | 3                    | LKC34                          |                             |
| ED3077       | 7                | 3                    | JU1586                         |                             |
| EG4946       | 7                | 3                    | NIC207                         |                             |
| JU1530       | 7                | 3                    | QX1215                         |                             |
| JU1652       | 7                | 3                    | JU1246                         |                             |
| JU323        | 7                | 3                    | JU830                          |                             |
| JU367        | 7                | 3                    | CX11314                        |                             |
| JU561        | 7                | 3                    | JU440                          |                             |
| JU774        | 7                | 3                    | QX1211                         |                             |
| JU847        | 7                | 3                    | ED3011                         |                             |
| MY18         | 7                | 3                    | CB4853_UK                      |                             |
| NIC200       | 7                | 3                    | NIC1                           |                             |
| PS2025       | 7                | 3                    | JU1491                         |                             |
| QG537        | 7                | 3                    | QX2268                         |                             |
| QG558        | 7                | 3                    | NIC3                           |                             |
| QX1216       | 7                | 3                    | JU315                          |                             |
| WN2010       | 7                | 3                    | CX11254                        |                             |
| WN2013       | 7                | 3                    | CB4851_CGC                     |                             |
| WN2016       | 7                | 3                    | MY10                           |                             |
| WN2019       | 7                | 3                    | CB4858_UK                      |                             |
| WN2020       | 7                | 3                    | JU1088                         |                             |