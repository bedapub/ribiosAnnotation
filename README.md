*ribiosAnnotation*: Annotate genes, mRNAs, and proteins with ribios and Roche Bioinformatics Infrastructure
===

## What is *ribiosAnnotation*?

*ribiosUtils* is a R package that performs feature annotation, for instance genes, mRNAs, and proteins, using either the Roche Bioinformatics infrastructure or the Bioconductor environment. It is distributed under the GPL-3 license.

## Installation

Run following commands in the R console.

```{R}
library(devtools)
devtools::install_github("bedapub/ribiosAnnotation")
```

## Annotation Backends

ribiosAnnotation supports two annotation backends:

- `mongodb` (default): online mode using the existing MongoDB instance using the Roche Bioinformatics infrastructure.
- `bioconductor`: offline mode using Bioconductor OrgDb packages and the babelgene package for ortholog mapping.

You can choose the backend per call:

```r
annotateGeneIDs(c(1, 2, 3), backend = "mongodb")
annotateGeneIDs(c(1, 2, 3), backend = "bioconductor")
```

Or set a global backend for the current R session:

```r
ribiosAnnotation::setAnnotationBackend("bioconductor")
ribiosAnnotation::getAnnotationBackend()

# Calls without explicit backend now use the global setting
annotateGeneSymbols(c("AKT1", "MAPK14"))
```

Or set a global default for the session with an environment variable:

```r
Sys.setenv(RIBIOS_ANNOTATION_BACKEND = "bioconductor")
annotateGeneSymbols(c("AKT1", "MAPK14"))
```

For `backend = "bioconductor"`, install:

- `AnnotationDbi`
- `org.Hs.eg.db` (human)
- `org.Mm.eg.db` (mouse)
- `org.Rn.eg.db` (rat)
- `org.Cf.eg.db` (dog)
- `org.Ss.eg.db` (pig)
- `org.Pt.eg.db` (chimpanzee)
- `org.Mmu.eg.db` (rhesus macaque)
- `babelgene` (for non-human to human orthologue mapping)

## CI Environment Variables

In CI, backend behavior is controlled with these environment variables:

- `RIBIOS_ANNOTATION_BACKEND`: default backend for function calls without an explicit `backend` argument (recommended CI default: `bioconductor`).
- `RIBIOS_TEST_BACKENDS`: comma-separated list of backends to exercise in tests.
	- `bioconductor`: run Bioconductor-only checks (works without MongoDB).
	- `bioconductor,mongodb`: also enable MongoDB-dependent tests.
- `RIBIOS_ANNOTATION_SECRETS_JSON`: path to the MongoDB secrets JSON file when MongoDB tests are enabled.

Typical setup in GitHub Actions:

- Standard CI (no MongoDB): set `RIBIOS_ANNOTATION_BACKEND=bioconductor` and `RIBIOS_TEST_BACKENDS=bioconductor`.
- Optional MongoDB CI job: additionally provide a valid secrets file and set `RIBIOS_TEST_BACKENDS=bioconductor,mongodb`.

## Contact

[Jitao David Zhang](mailto:jitao_david.zhang@roche.com) maintains and develops *ribiosUtils* and other ribios packages in memory of Clemens Broger, a pioneer of bioinformatics and cheminformatics in drug discovery, a man true to himself. Jitao David Zhang thanks Balazs Banfai, Marco Berrera, and Roland Schmucki for their help and input.
