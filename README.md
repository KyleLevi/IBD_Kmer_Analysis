# IBD_Kmer_Analysis

## Python Implementaion - Thanks Rob!
The python kmerbias should function the same as the perl kmerbias, both depend on the output from GetNucFrequency_PerSeq_varK.pl (GCF_000403175.4mers is an example output file)

```
python3 kmerbias.py -f GCF_000403175.4mers -k GGCC,GCGC -v > pyout
perl kmerbias.pl GCF_000403175.4mers "GGCC,GCGC"  > plout
```

plout and pyout are not identical because of the floating point math, but they are close to several significant digits.


## Perl Version

The code should be able to run the code like this:

```perl GetNucFrequency_PerSeq_varK.pl ./GCF_000403175.fna 4 > GCF_000403175.4mers```

and then you can count the bias:

```perl kmerbias.pl GCF_000403175.4mers "GGCC,GGCG,GGGC"```


The first output file is GCF_00403175.4mers

The second command should output:

| Matrix | NZ_KE159482.1    | NZ_KE159483.1     | NZ_KE159484.1     | NZ_KE159485.1     | NZ_KE159486.1     | NZ_KE159487.1     | NZ_KE159488.1     | NZ_KE159489.1     | NZ_KE159490.1     | NZ_KE159491.1     | NZ_KE159492.1     | 
|--------|------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------| 
| GGCC   | 0.72379167166026 | 0.564725231578111 | 1.06033886652903  | 0.700487472590855 | 0.722849965243236 | 0.744502198607545 | 0.718905768185197 | 0.627581109687921 | 0.638351277899619 | 0.907762632169645 | 0                 | 
| GGCG   | 1.02371445359977 | 1.01361527445047  | 0.880629901562083 | 1.01214386205571  | 0.99189604408644  | 1.00392096757395  | 1.01391615871228  | 1.03649903331046  | 1.03819590211593  | 0.875969559120189 | 0.934879508443214 | 
| GGGC   | 1.09172008969991 | 1.11655103110103  | 0.956758304197849 | 1.10617995726091  | 1.08601140376983  | 1.09210234817395  | 1.10214747927039  | 1.15141594439126  | 1.21693684689459  | 1.09103255691336  | 1.22430819416776  | 




