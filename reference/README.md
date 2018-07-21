
## References

### Reference models

**mBAL_g1g4combo_logRegModel_MODEL_RNA_rank_ispatho**
This model was trained on mBAL respiratory samples using only the input features
- RNA
- rank 
- is_pathogen

The performance on the mBAL study was:
Derivation LOOCV AUC = 0.96 +/- 0.03
Derivation LOPO-CV AUC = 0.91 (0.78 - 0.99)
Validation AUC = 0.96 (0.91 - 1.00) 


**mBAL_g1g4combo_logRegModel_MODEL_RNA_rank_ispatho_isvirus**
This model was trained on mBAL respiratory samples using only the input features
- RNA
- rank 
- is_pathogen
- is_virus

The performance on the mBAL study was:
Derivation LOOCV AUC = 0.95 +/- 0.04
Derivation LOPO-CV AUC = 0.90 (0.76 - 0.99)
Validation AUC = 0.97 (0.93 - 1.00) 


### Feature-calling references

**pathogens_bangladesh_official.txt** - this is the list of reference pathogens used for the Bangladesh study.

**viruses_bangldesh_uniq.txt** - this is a list of all viruses identified in the Bangladesh cohort.


**known_resrpiatory_pathogens.txt** - this is the reference list of respiratory pathogens described in the manuscript.


