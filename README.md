# causality-toolbox
Functions for calculating information theory measures from time-series, with the purpose of distinguishing the cause and the effect variables on a physical process.

This version is written in Matlab/Octave.

* cami.m: calculates, for bivariate systems, the Causal Mutual Information (CaMI) in both directions, the Mutual Information, the Transfer Entropy in both directions and the Directionality Index

* camir.m: calculates, for bivariate sytems, the Causal Mutual Information Rate (CaMIR) in both directions, the Mutual Information Rate and the Transfer Entropy Rate in both directions

* totals.m: calculates, for multivariate systems, the Total Correlation (generalization of Mutual Information) and the Joint Entropy of the entire system

Details of each program and usage example is available by typing in the command window: help [function-name]

-----------------
(C) Arthur Valencio* and Dr Murilo S. Baptista

ICSMB, University of Aberdeen

* Support: CNPq (Brazil)

This package is available as is, without any warranty. Use it at your own risk.

-------------------
If useful, please cite:

Arthur Valencio and Murilo S. Baptista. Causality Toolbox: functions for calculating information theory measures from time-series. Open source codes for Matlab. 2018.  


* Bibtex entry:

@misc{cami,
author={Valencio, Arthur and Baptista, Murilo da Silva},
title={Causality Toolbox: functions for calculating information theory measures from time-series},
note={Open source codes for Matlab. Available at \url{https://github.com/artvalencio/causality-toolbox/}.},
year={2018}
}
