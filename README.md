# RbaseDB
## SNP effects based on annotation

Genetic (DNA) sequence is annotated with effects of a change in each nucleotide base. This allows to query: what if at position X there would be a change to [A|G|C|T]?


# Design

###  BASE ONTOLOGY

### Concept by Mark Blaxter, Perl code by John Davey translated to R
### by Emanuel Heitlinger:

- A: phase 1 any change is nonsynonymous
- B: phase 2 any change is nonsynonymous
- C: phase 3 any change is nonsynonymous
- D: phase 1 change to CT is nonsynonymous
- E: phase 2 change to CT is nonsynonymous
- F: phase 3 change to CT is nonsynonymous
- G: phase 1 change to AG is nonsynonymous
- H: phase 2 change to AG is nonsynonymous
- I: phase 3 change to AG is nonsense
- K: phase 1 change to GT is nonsynonymous
- L: phase 2 change to A is nonsense, to anything else is nonsynonymous
- J: phase 3 change to G is nonsynonymous
- M: phase 3 change to G is nonsense, to A is nonsynonymous
- N: phase 3 any change synonymous
- O: phase 1 change to T nonsense, others nonsynonymous
- P: phase 3 change to AG is nonsynonymous
- Q: phase 1 change to T nonsense, to G nonsynonymous
- R: phase 2 change to AG nonsense, others nonsynonymous
- S: phase 3 change to A nonsense, others nonsynonymous
- T: phase 3 change to A nonsense, G nonsynonymous

-W: all changes are unknown # EH added 08/23/2011

|      |  a       |    g       |    c       |    t
|------|----------|------------|------------|----------
| a    |aaa K OBF |  aga R QBF |  aca T ABN |  ata I ABJ
|      |aag K OBF |  agg R KBF |  acg T ABN |  atg M ABC
|      |aac N ABP |  agc S ABP |  acc T ABN |  atc I ABJ
|      |aat N ABP |  agt S ABP |  act T ABN |  att I ABJ
|      |          |            |            |
| g    |gaa E OBF |  gga G OBN |  gca A ABN |  gta V ABN
|      |gag E OBF |  ggg G ABN |  gcg A ABN |  gtg V ABN
|      |gac D ABP |  ggc G ABN |  gcc A ABN |  gtc V ABN
|      |gat D ABP |  ggt G ABN |  gct A ABN |  gtt V ABN
|      |          |            |            |
| c    |caa Q OBF |  cga R QBN |  cca P ABN |  cta L GBN
|      |cag Q OBF |  cgg R KBN |  ccg P ABN |  ctg L GBN
|      |cac H ABP |  cgc R ABN |  ccc P ABN |  ctc L ABN
|      |cat H ABP |  cgt R ABN |  cct P ABN |  ctt L ABN
|      |          |            |            |
| t    |taa * AEF |  tga * AEC |  tca S ARN |  tta L GRF
|      |tag * ABF |  tgg W ALS |  tcg S ALN |  ttg L GLF
|      |tac Y ABI |  tgc C ABT |  tcc S ABN |  ttc F ABP
|      |tat Y ABI |  tgt C ABT |  tct S ABN |  ttt F ABP




