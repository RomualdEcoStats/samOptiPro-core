## ---- Internal functions (1) ----

customize_samplers <- function(conf, low_ess=NULL, strategy=list(sampler_for_low_ess='slice', block_correlated=TRUE, block_size=3)){ if(!is.null(low_ess)) for(nm in low_ess){ if(conf$isSamplerAssigned(nm)) conf$removeSampler(nm); conf$addSampler(nm, type='slice') } ; conf }


