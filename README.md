# liddell
Perform Liddell's exact test on matched pairs.<br/>
Paired proportions have traditionally been compared using McNemar's test
but an exact alternative due to Liddell is preferable. The usual reason
for using approximate methods, when exact methods are available, is that
the former are "quick" if sometimes "dirty". 
There is no such justification in the present circumstances.
Liddell FDK. Simplified exact analysis of case-referent studies; matched
pairs; dichotomous exposure. Journal of Epidemiology and Community Health 1983;37:82-84.

Syntax: 	liddell(x,alpha)
     
    Inputs:
          X - 2x2 data matrix 
          ALPHA (default 0.05) 
    Outputs:
          - Point estimate of relative risk (R')
          - Exact confidence interval
          - F statistics with p-value
          - Power

  Example:
In the following example, a researcher attempts to determine if a drug
has an effect on a particular disease. 

                     Drug
                 +         -
            --------------------
        +   |   101   |   59   |
            |-------------------   Placebo
        -   |   121   |   33   |
            --------------------
                                      

  x=[101 59; 121 33];

  Calling on Matlab the function: 
            liddell(x)

  Answer is:

Liddell's exact test<br/>
Point estimate of relative risk (R') =  2.0508<br/>
Exact 95% confidence interval = 1.4904 to 2.8493<br/>
F = 2.0167 p-value (two-side)= 0.00000443<br/>
alpha = 0.0500  Zb = 2.7566  Power (2-tails) = 0.0058<br/>

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2008) Liddell's test: Perform Liddell's exact test on
matched pairs
http://www.mathworks.com/matlabcentral/fileexchange/22024
