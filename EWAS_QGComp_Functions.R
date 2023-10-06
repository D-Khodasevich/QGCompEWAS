# EWAS QGComp Application
# Author: Dennis Khodasevich


# ewas_qgcomp: function for iteratively running quantile g-computation across dna methylation array data 
### pheno: matrix containing exposure and covariate data
### meth: dna methylation matrix with CpG sites as rows and Individual Samples as columns
####### rownames should reflect CpG name, colnames should reflect sample ID
####### order of rows in pheno must match order of rows in meth
### mix_comp: vector of mixture components
### covars: vector of covariates to adjust for
### mval_conversion: whether to perform an m-value conversion on the DNAm data. Default is FALSE, set to TRUE if using a beta value matrix
### qval: value of q for quantile g computation function. Default is 4
### output_type: whether to output only the mixture summary output or include weight summaries as well. 
###### basic: only output is the summary dataframe including probeID, psi, se, CI, t-value, and p-value
###### full: basic output plus complete mixture weight summaries (note: higher run time for greater mixture components)

ewas_qgcomp <- function(pheno, meth, mix_comp, covars, 
                        mval_conversion=FALSE, qval=4, output_type = "basic", ...) {
  # create empty storage for results  
  cpg_list <- rownames(meth)
  nprobes <- nrow(meth)
  results <- data.frame(probeID = rownames(meth), 
                        beta = rep(0, times=nprobes), 
                        se = rep(0, times=nprobes), 
                        lowCI = rep(0, times=nprobes), 
                        upCI = rep(0, times=nprobes), 
                        pval = rep(0, times=nprobes), 
                        tval = rep(0, times=nprobes))
  
  if(output_type == "full") {
    pos_weights <- data.frame(probeID = rep(0, times=nprobes*length(mix_comp)), 
                              component = rep(0, times=nprobes*length(mix_comp)), 
                              weight = rep(999, times=nprobes*length(mix_comp)), 
                              direction = rep("Positive", times=nprobes*length(mix_comp)))
    neg_weights <- data.frame(probeID = rep(0, times=nprobes*length(mix_comp)), 
                              component = rep(0, times=nprobes*length(mix_comp)), 
                              weight = rep(999, times=nprobes*length(mix_comp)), 
                              direction = rep("Negative", times=nprobes*length(mix_comp)))
    pos_weights$probeID <- as.character(pos_weights$probeID)
    pos_weights$component <- as.character(pos_weights$component)
    neg_weights$probeID <- as.character(neg_weights$probeID)
    neg_weights$component <- as.character(neg_weights$component)
  } # weight dataframe initialization
  
  options(stringsAsFactors = FALSE, scipen = 999)
  
  pos_limit = 0
  neg_limit = 0
  
  for(i in 1:length(cpg_list)) {
    cpg_i <- (meth[i, ])
    current <- pheno
    cur_cpg <- cpg_list[i]
    
    current$cpg_i <- paste0(cpg_i)
    current$cpg_i <- as.numeric(current$cpg_i)
    
    if(mval_conversion == TRUE) {
      current$cpg_i <- log2((current$cpg_i)/(1-current$cpg_i))
    } # m value conversion
    
    qc.fit <- qgcomp.noboot(cpg_i ~.,dat=current[,c(mix_comp, covars, 'cpg_i')],    
                            expnms=mix_comp, q=qval, family=gaussian())
    
    tstat <- qc.fit[["tstat"]][[2]] 
    sums <- quiet(summary(qc.fit)$coefficients  )                      
    sums <- sums[2, ]
    PSI = sums[1]
    SE = sums[2]
    LOW = sums[3]
    UP = sums[4]
    PVAL = sums[5]
    TSTAT = tstat 
    
    i = as.integer(i)
    
    data.table::set(results, i, 2L, PSI)
    data.table::set(results, i, 3L, SE)
    data.table::set(results, i, 4L, LOW)
    data.table::set(results, i, 5L, UP)
    data.table::set(results, i, 6L, PVAL)
    data.table::set(results, i, 7L, TSTAT)
    
    if(output_type == "full") {
      pos <- as.data.frame(qc.fit$pos.weights)
      metp <- rownames(pos)
      neg <- as.data.frame(qc.fit$neg.weights)
      metn <- rownames(neg)
      
      if(length(metp > 0)) {
        
        for(j in 1:nrow(pos)) {
          addon = pos_limit + j
          poscur = pos[j,1]
          metcur = metp[j]
          addon = as.integer(addon)
          
          data.table::set(pos_weights, addon, 3L, poscur)
          data.table::set(pos_weights, addon, 2L, metcur)
          data.table::set(pos_weights, addon, 1L, cur_cpg)
        }
        pos_limit = pos_limit + j} else{
          pos_limit = pos_limit + 1
          pos_limit = as.integer(pos_limit)
          data.table::set(pos_weights, pos_limit, 3L, value=777)
          data.table::set(pos_weights, pos_limit, 2L, value=777)
          data.table::set(pos_weights, pos_limit, 1L, cur_cpg)
          
        }
      
      if(length(metn > 0)) {
        
        for(f in 1:nrow(neg)) {
          addon_n = neg_limit + f
          negcur = neg[f,1]
          negmetcur = metn[f]
          addon_n = as.integer(addon_n)
          
          data.table::set(neg_weights, addon_n, 3L, negcur)
          data.table::set(neg_weights, addon_n, 2L, negmetcur)
          data.table::set(neg_weights, addon_n, 1L, cur_cpg)
        }
        neg_limit = neg_limit + f} else {
          neg_limit = neg_limit + 1
          neg_limit = as.integer(neg_limit)
          data.table::set(neg_weights, neg_limit, 3L, value=777)
          data.table::set(neg_weights, neg_limit, 2L, value=777)
          data.table::set(neg_weights, neg_limit, 1L, cur_cpg)
        }
    } # saving weight summaries
    
    if(i %% 5000 == 0) {cat("Progress:", 100*i/nprobes, "%\n")} # progress bar, made for 450k sized data, adjust as necessary
  } # iterative qgcomp function
  
  # trim weights 
  if(output_type == "full") {
    pos_weights <- pos_weights %>% 
      filter(weight != 999)
    pos_weights[, 2:3][pos_weights[, 2:3] == 777] <- NA
    neg_weights <- neg_weights %>% 
      filter(weight != 999)
    neg_weights[, 2:3][neg_weights[, 2:3] == 777] <- NA
    ewas_qgcomp_fit <- list(results, pos_weights, neg_weights)
    nms <- c("results", "pos_weights", "neg_weights")
    names(ewas_qgcomp_fit) <- nms
    return(ewas_qgcomp_fit)
  } else {
    return(results)
  }
}

# output guide
### results: standard qgcomp summary measures: "probeID" "beta"    "se"      "lowCI"   "upCI"    "pval"    "tval"   
### pos_weights: positive component weight summaries
### neg_weights: negative component weight summaries




# quiet function
### dependency for ewas_qgcomp function, quiets summary output
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 



# ewas_qgcomp.boot: function for iteratively running quantile g-computation across dna methylation array data 
### pheno: matrix containing exposure and covariate data
### meth: dna methylation matrix with CpG sites as rows and Individual Samples as columns
####### rownames should reflect CpG name, colnames should reflect sample ID
####### order of rows in pheno must match order of rows in meth
### mval_conversion: whether to perform an m-value conversion on the DNAm data. Default is FALSE, set to TRUE if using a beta value matrix
### qval: value of q for quantile g computation function. Default is 4
### output_type: whether to output only the mixture summary output or include weight summaries as well. 
###### basic: only output is the summary dataframe including probeID, psi, se, CI, t-value, and p-value
###### full: basic output plus complete mixture weight summaries (note: higher run time for greater mixture components)

ewas_qgcomp.boot <- function(pheno, meth, mix_comp, covars, 
                             mval_conversion=FALSE, qval=4, output_type = "basic", 
                             Bval=1000, seedval=123,...) {
  # create empty storage for results  
  cpg_list <- rownames(meth)
  nprobes <- nrow(meth)
  results <- data.frame(probeID = rownames(meth), 
                        beta = rep(0, times=nprobes), 
                        se = rep(0, times=nprobes), 
                        lowCI = rep(0, times=nprobes), 
                        upCI = rep(0, times=nprobes), 
                        pval = rep(0, times=nprobes), 
                        tval = rep(0, times=nprobes))
  
  if(output_type == "full") {
    component_effects <- data.frame(probeID = rep(0, times=nprobes*length(mix_comp)), 
                                    component = rep(0, times=nprobes*length(mix_comp)), 
                                    weight = rep(999, times=nprobes*length(mix_comp)))
    component_effects$probeID <- as.character(component_effects$probeID)
    component_effects$component <- as.character(component_effects$component)
  } # weight dataframe initialization
  
  options(stringsAsFactors = FALSE, scipen = 999)
  pos_limit = 0
  
  for(i in 1:length(cpg_list)) {
    cpg_i <- (meth[i, ])
    current <- pheno
    cur_cpg <- cpg_list[i]
    
    current$cpg_i <- paste0(cpg_i)
    current$cpg_i <- as.numeric(current$cpg_i)
    
    if(mval_conversion == TRUE) {
      current$cpg_i <- log2((current$cpg_i)/(1-current$cpg_i))
    } # m value conversion
    
    qc.fit <- qgcomp.boot(cpg_i ~.,dat=current[,c(mix_comp, covars, 'cpg_i')],    
                          expnms=mix_comp, q=qval, family=gaussian(), B=Bval, seed=seedval)
    
    tstat <- qc.fit[["tstat"]][[2]] 
    sums <- quiet(summary(qc.fit)$coefficients  )                      
    sums <- sums[2, ]
    PSI = sums[1]
    SE = sums[2]
    LOW = sums[3]
    UP = sums[4]
    PVAL = sums[5]
    TSTAT = tstat 
    
    i = as.integer(i)
    
    data.table::set(results, i, 2L, PSI)
    data.table::set(results, i, 3L, SE)
    data.table::set(results, i, 4L, LOW)
    data.table::set(results, i, 5L, UP)
    data.table::set(results, i, 6L, PVAL)
    data.table::set(results, i, 7L, TSTAT)
    
    if(output_type == "full") {
      pos <- as.data.frame(qc.fit$fit$coefficients)
      pos$metabolite <- rownames(pos)
      co_length = length(mix_comp) + 1 
      pos <- as.data.frame(pos[c(2:co_length), ])
      metp <- rownames(pos)
      
      for(j in 1:length(mix_comp)) {
        addon = pos_limit + j
        poscur = pos[j,1]
        metcur = metp[j]
        addon = as.integer(addon)
        
        data.table::set(component_effects, addon, 3L, poscur)
        data.table::set(component_effects, addon, 2L, metcur)
        data.table::set(component_effects, addon, 1L, cur_cpg)
      }
      pos_limit = pos_limit + j
    }
    
    if(i %% 5 == 0) {cat("Progress:", 100*i/nprobes, "%\n")} # progress bar
  } # iterative qgcomp function
  
  if(output_type == "full") {
    ewas_qgcomp_fit.boot <- list(results, component_effects)
    nms <- c("results", "component_effects")
    names(ewas_qgcomp_fit.boot) <- nms
    return(ewas_qgcomp_fit.boot)
  } else {
    return(results)
  }}
