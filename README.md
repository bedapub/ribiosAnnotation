*ribiosAnnotation*: Annotate genes, mRNAs, and proteins with ribios and Roche Bioinformatics Infrastructure
===

## What is *ribiosAnnotation*?

*ribiosUtils* is a R package that performs feature annotation, for instance genes, mRNAs, and proteins, using Roche Bioinformatics infrastructure. It is distributed under the GPL-3 license.

## Installation

Run following commands in the R console.

```{R}
library(devtools)
devtools::install_github("bedapub/ribiosAnnotation")
```

The user has to specify a secrets file in JSON format (the template is accessible here [inst/secrets/secrets-template.json](inst/secrets/secrets-template.json), with invalid credentials). The secrets file can be specified with the `RIBIOS_ANNOTATION_SECRETS_JSON` environmental variable. If not, ribiosAnnotation looks for the file at `~/.credentials/ribiosAnnotation-secrets.json`.

## Contact

[Jitao David Zhang](mailto:jitao_david.zhang@roche.com) maintains and develops *ribiosUtils* and other ribios packages in memory of Clemens Broger, a pioneer of bioinformatics and cheminformatics in drug discovery, a man true to himself. Jitao David Zhang thanks Balazs Banfai, Marco Berrera, and Roland Schmucki for their help and input.
