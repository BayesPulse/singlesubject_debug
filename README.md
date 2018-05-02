
# singlesubject\_debug

There is a persistent and evasive bug in the single-subject model code.  This
repo is an attempt to root out the issue.

We are using the C source code from Matt's thesis work as the codebase for
debugging.  The pulsatile package and the new libpulsatile library/package is
based off of this same codebase.

Two symptoms exist so far.  The first is 'sticky points' in the chains of mean
width and the SD of pulse widths for essentially any dataset.  The second is
random noise/lack of convergence  in the posterior chains of pulse locations,
specifically in noisier datasets.



