## Updates 4/20/17

1. Add a new function mhglm_ml, which can fit multilevel hiererchical models of two levels or more ('_ml' indicates multilevel).
2. Define mhglm_ml's returned object as class 'mhglm_ml'. Implement various class functions for it. See changes in mhglm.R file.
3. In test/testthat/test-mhglm.multilevel.R, write new tests for:
- whether the results of mhglm and mhglm_ml on two levels model are same (the returned object has different structure, but fixef, ranef, covariance matrix estimates are same).
- whether mhglm_ml works on three-levels model
- whether different formats of formulas are translated as desired
4. In linalg.R, updated a few lines to make pseudo.inverse more robust to numerical errors

## TODO 4/20/17
1. Test out the formula '~(1+Days||Subject)' for both mhglm and mhglm_ml, where '||' indicates diagonal covariance matrix. 
2. mhglm_ml do not have parallel capability. Need to ignore that flag. 
3. mhglm_ml breaks when standardize = TRUE in the control. The default is set to standardize = FALSE. Need to fix. 
4. Test mhglm_ml on more than three levels.
5. Test mhglm_ml on some extrem datasets (small dataset, very skewed dataset, etc).


