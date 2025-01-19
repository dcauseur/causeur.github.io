
library(FrF2)

### Plan complet pour facteurs à deux modalités

plan3 <- FrF2(nruns=8,nfactors=3,factor.names=paste0("F",1:3),
                    replications=1)
print(plan3)

### Plan 2^{4-1}

plan4 <- FrF2(nruns=8,nfactors=4,factor.names=paste0("F",1:4))
plan4
design.info(plan4)

### Plan 2^{5-2}

plan5 <- FrF2(nruns=8,nfactors=5,factor.names=paste0("F",1:5))
plan5
design.info(plan5)