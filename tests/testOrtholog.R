library(ribiosAnnotation)

## find out all non-human orthologues of MAPK14
## querydb("SELECT RO_GENE_ID1 as HG, RO_GENE_ID2 as NHG FROM GTI_ORTHOLOGS WHERE RO_GENE_ID1='1432' AND TAX_ID2!='9606'", "bia")

## find out human ortholog of mouse Mapk14
## querydb("SELECT RO_GENE_ID1 as HG, RO_GENE_ID2 as NHG FROM GTI_ORTHOLOGS WHERE RO_GENE_ID2='26416' AND TAX_ID1='9606'", "bia")

mouse.geneids <- c("20466", "240411", "59052", "72503", "15951", "170743")
system.time(mouse.ho <- humanOrthologs(mouse.geneids))
system.time(mouse.ho.uniq <- humanUniqOrtholog(mouse.geneids))
stopifnot(identical(length(mouse.geneids), length(mouse.ho)))
stopifnot(identical(length(mouse.geneids), length(mouse.ho.uniq)))
stopifnot(is.list(mouse.ho))
stopifnot(is.vector(mouse.ho.uniq))


human.geneids <- c("10013", "10014", "10016", "1001", "25895")
system.time(nho <- nonhumanOrthologs(human.geneids))
system.time(nho.mouse <- nonhumanOrthologs(human.geneids, taxid=c(10090)))
system.time(nho.mouse.rat <- nonhumanOrthologs(human.geneids, taxid=c(10090, 10116))) ## same results as setting TAXID=NULL when only rat and mouse orthologs exist.
system.time(nho.mouse.uniq <- nonhumanUniqOrtholog(human.geneids, taxid=c(10090)))

stopifnot(identical(length(human.geneids), length(nho)))
stopifnot(identical(length(human.geneids), length(nho.mouse)))
stopifnot(identical(length(human.geneids), length(nho.mouse.rat)))
stopifnot(is.vector(nho.mouse.uniq))
## NA in uniqOrtholog should be NULLs in ortholog lists
stopifnot(identical(sum(is.na(nho.mouse.uniq)),
                    sum(sapply(nho.mouse, length)==0)))
