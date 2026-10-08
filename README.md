EpiDelays: A Software for Estimation of Epidemiological Delays
================
Oswaldo Gressani (<oswaldo.gressani@uhasselt.be>)

**EpiDelays** is a package that can be used to fit epidemiological delay
distributions to singly and doubly interval-censored data. The
`parfitml` routine can be used to fit parametric delay distributions
(`gamma`, `lognormal`, `weibull`, `gaussian`, `skewnorm`) via maximum
likelihood. The `nonparfit` routine uses the nonparametric methodology
developed by Gressani and Hens (2025) to obtain estimates of key
epidemiological delay features without imposing any parametric
assumptions. The `nonparfit` routine relies on the nonparametric
bootstrap to compute estimated standard errors and confidence interval
of delay features. For the `parfitml` routine, quantification of
uncertainty can be done via the nonparametric bootstrap. For large
samples, an option in EpiDelays is to rely on the asymptotic normality
of the maximum likelihood estimator. In that case, approximate
confidence intervals are computed by simulating draws from the Gaussian
distribution following from asymptotic MLE theory (Mandel, 2013). This
is usually a faster alternative than the nonparametric bootstrap
whenever large sample sizes are involved.

This package is currently under construction and additional features
will be added in future releases.

#### Associated literature

1.  Gressani, O. and Hens, N. (2025). Nonparametric serial interval
    estimation with uniform mixtures. *PLoS Computational Biology*,
    **21**(8):e1013338. <https://doi.org/10.1371/journal.pcbi.1013338>

2.  Mandel, M. (2013). Simulation-based confidence intervals for
    functions with complicated derivatives. *The American Statistician*,
    **67**(2), 76-81.

#### Package version

This is version 0.0.3 - “Parametric wave 2”.<br> Release date:
2026-10-08 (October 8, 2026).

#### Authors and contributors

1.  Oswaldo Gressani (author, maintainer) <br>
2.  Dongxuan Chen (contributor)

#### Acknowledgments

This project is supported by the Research Foundation - Flanders (FWO)
through a Senior Postdoctoral Fellowship awarded to Oswaldo Gressani
(grant no. 1236527N). This project was supported by the VERDI project
(101045989) and the ESCAPE project (101095619), funded by the European
Union. Views and opinions expressed are however those of the authors
only and do not necessarily reflect those of the European Union or
European Health and Digital Executive Agency (HADEA). Neither the
European Union nor the granting authority can be held responsible for
them. We also acknowledge the BE-PIN project (contract nr.
TD/231/BE-PIN) funded by BELSPO (Belgian Science Policy Office) as part
of the POST-COVID programme.

#### License

Copyright (C) 2025-2026 Oswaldo Gressani. All rights reserved.
