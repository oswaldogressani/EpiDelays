EpiDelays: A Software for Estimation of Epidemiological Delays
================
Oswaldo Gressani (<oswaldo.gressani@uhasselt.be>)

**EpiDelays** is a package that can be used to fit epidemiological delay
distributions to interval-censored data. The `parfitml` routine can be
used to fit classic delay distributions (`gamma`, `lognormal`,
`weibull`, `gaussian`) via maximum likelihood. The `nonparfit` routine
uses the nonparametric methodology developed by Gressani and Hens (2025)
to obtain estimates of key epidemiological delay features without
imposing any parametric assumptions. For both routines, the
nonparametric bootstrap can be implemented to compute standard errors
and confidence intervals for often reported epidemiological delay
features.

This package is currently under construction and additional features
will be added in future releases.

#### Associated literature

1.  Gressani, O. and Hens, N. (2025). Nonparametric serial interval
    estimation with uniform mixtures. *PLoS Computational Biology*,
    **21**(8):e1013338. <https://doi.org/10.1371/journal.pcbi.1013338>

#### Package version

This is version 0.0.2 - “Parametric wave 1”.<br> Release date:
2025-12-14 (December 14, 2025).

#### Authors and contributors

1.  Oswaldo Gressani (author, maintainer) <br>
2.  Dongxuan Chen (contributor)

#### Acknowledgments

This project is supported by the VERDI project (101045989) and the
ESCAPE project (101095619), funded by the European Union. Views and
opinions expressed are however those of the authors only and do not
necessarily reflect those of the European Union or European Health and
Digital Executive Agency (HADEA). Neither the European Union nor the
granting authority can be held responsible for them. This project is
also supported by the BE-PIN project (contract nr. TD/231/BE-PIN) funded
by BELSPO (Belgian Science Policy Office) as part of the POST-COVID
programme.

#### License

Copyright (C) 2025-2026 Oswaldo Gressani. All rights reserved.
