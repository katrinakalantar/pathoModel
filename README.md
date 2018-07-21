
## README

Katrina Kalantar
July 2018

Code to run the pathogen v. commensal microbe logistic regression model.

```
python [IDseq report filename] [reference known pathogens filename] [model]
```

For example, running from the github root directory, we can predict the pathogen output for sample 0058 as follows:

```
python3 ./scripts/predict_pathogens.py ./data/070518/rapid-response-007_reports/chrf_rna_0058_s58.csv ./reference/pathogens_bangladesh_official.txt ./reference/mBAL_g1g4combo_logRegModel_MODEL_RNA_rank_ispatho
```

The expected result is a matrix of all pathogens:

```
   RNAvalue  ranks                     microbe  microbe_genus  pathogenic_red   score
0  1.021189      0  Mycobacterium tuberculosis           1763            True   0.458405
1  0.000000      1           Torque teno virus     -200687329           False   0.010700
```
The "score" column contains the logistic regression predicted probability that this is an pathogen. In benchmark studies p > .46 has been a reasonable threshold for identifying single pathogen cases. To capture coinfections, a more lenient threshold of 0.2 may be useful.

