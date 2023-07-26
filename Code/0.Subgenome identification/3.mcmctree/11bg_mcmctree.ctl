          seed = -1
       seqfile = 11bg_alignment.txt
      treefile = 11bg.trees
      outfile = 11bg_time

         ndata = 1
       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<1.0'  * safe constraint on root age, used if no fossil for root.

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0    * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 5  * gamma prior for overall rates for genes
  sigma2_gamma = 1 10   * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.08 0.08  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 200000
      sampfreq = 100
       nsample = 20000

*** Note: Make your window wider (100 columns) before running the program.
