# Ising study

|L|T|fault|energy/spin|abs magnetization|susceptibility|heat capacity|
|---|---|---|---|---|---|---|
|2|1.5|bias|-2|1|0|0|
|2|1.5|coupled|-1.94167|0.980313|2.27323|0.201199|
|2|1.5|healthy|-1.93944|0.97941|2.05587|0.208088|
|2|2.26919|bias|-2|1|0|0|
|2|2.26919|coupled|-1.69389|0.89849|1.49649|0.404797|
|2|2.26919|healthy|-1.69424|0.898646|1.5006|0.406456|
|2|3.5|bias|-2|1|0|0|
|2|3.5|coupled|-1.23354|0.744931|0.773535|0.319175|
|2|3.5|healthy|-1.22722|0.744809|0.771775|0.324808|
|3|1.5|bias|-2|1|0|0|
|3|1.5|coupled|-1.96435|0.991088|0.0113862|0.121453|
|3|1.5|healthy|-1.94727|0.98473|2.0343|0.212119|
|3|2.26919|bias|-2|1|0|0|
|3|2.26919|coupled|-1.7158|0.914823|0.0957247|0.386629|
|3|2.26919|healthy|-1.61611|0.872716|3.0636|0.615367|
|3|3.5|bias|-2|1|0|0|
|3|3.5|coupled|-1.23369|0.745625|0.122313|0.233494|
|3|3.5|healthy|-0.922593|0.610887|1.19945|0.401777|
|4|1.5|bias|-1.875|0.875|0|0|
|4|1.5|coupled|-1.95811|0.98918|0.0149951|0.143543|
|4|1.5|healthy|-1.94801|0.9851|0.216916|0.217774|
|4|2.26919|bias|-2|1|0|0|
|4|2.26919|coupled|-1.69727|0.906688|0.101652|0.378681|
|4|2.26919|healthy|-1.56336|0.84296|4.75122|0.782853|
|4|3.5|bias|-2|1|0|0|
|4|3.5|coupled|-1.17595|0.692821|0.189947|0.22135|
|4|3.5|healthy|-0.772969|0.488429|1.4534|0.370603|
|8|1.5|bias|-1.9375|0.90625|0|0|
|8|1.5|coupled|-1.87948|0.905735|0.0530506|0.247264|
|8|1.5|healthy|-1.94857|0.985678|0.0295384|0.211425|
|8|2.26919|bias|-2|1|0|0|
|8|2.26919|coupled|-1.62403|0.821636|0.145988|0.290439|
|8|2.26919|healthy|-1.47824|0.766589|12.7513|1.16513|
|8|3.5|bias|-2|1|0|0|
|8|3.5|coupled|-1.25417|0.636361|0.149171|0.102643|
|8|3.5|healthy|-0.665985|0.252739|1.72043|0.260828|

Random-scan Metropolis, J=kB=1, periodic square lattice. L2 has parallel periodic bonds; both enumerator and kernel use that convention. Exact enumeration covers L2,L3,L4 (up to 65536 configurations). L8 has no exact finite-size reference here. Infinite-lattice Tc is only a regime marker. Susceptibility uses (mean M²-mean M²)/(NT), explicitly not the absolute-M convention. Heat capacity uses energy variance/(NT²). Finite-chain plug-in nonlinear estimates retain sampling bias.

32 independent replica chains per cell; healthy and faulty arms share preallocated draws within a seed group. Unit=chain, not sweep. 1200 sweeps,300 burn-in; no mixing guarantee. Eight Wolff chains per L4/L8/temperature use1000 cluster updates,200 burn-in. Update counts are not work-matched sweeps, so no speedup claim is justified. Traces retain cluster sizes.

`analysis.json` contains effect CIs and Holm-adjusted paired energy tests. The CIs are pointwise, unadjusted across cells; only Holm p-values control the named family. Repeated temperatures/lattices are separate tests, not independent confirmation of a universal law.

![Ising](figures/ising.png)

Literature: [Ferrenberg et al.1992](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.3382), [Wolff1989](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.62.361/fulltext). RNG/algorithm interactions and critical slowing down predate this work.
