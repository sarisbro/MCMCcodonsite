This method is described in Aris-Brosou, S. 2006. Identifying sites under positive selection with uncertain parameter estimates. Genome. 49:767-776. 

- mcmc_041105_M2a.pl: use this to run M2a

- mcmc_041105_Mx.pl: all other models (M0, M3, M7, M8) can be run with this code

- mcmc_041105_Mx.exe: a pre-compiled version of the above for MS-DOS (and Windows I suppose)

- simulmcmc_041208.pl: script to perform simulations

Note that you will need codeml (Ziheng Yang's PAML 3.14 prior to September 2004, available from this page: http://abacus.gene.ucl.ac.uk/software/pamlOld.html) in your path. 

To run the MCMC sampler, format the codeml control file as you would do to run a maximum likelihood analysis; the settings of the sampler can be changed from the first lines of the scripts.
