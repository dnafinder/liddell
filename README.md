[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/liddell)

📘 Overview
liddell is a MATLAB implementation of Liddell's exact test for matched pairs, an exact alternative to McNemar's test for paired binary outcomes. It is designed for case-referent or paired designs where each subject is observed under two conditions (e.g., drug vs placebo, pre vs post), and each observation is classified as a dichotomous outcome (e.g., disease vs no disease, positive vs negative).

The test uses the 2×2 matched-pairs table:

             Condition 2
                  +      -
Condition 1  +   a      b
            -   c      d

The information about treatment effect is carried by the discordant pairs b and c. Liddell's exact test provides an exact inference on the relative risk, its confidence interval, and an associated F statistic with a two-sided p-value.

✨ Features
- Exact test for matched pairs (paired binary data)
- Calculates:
  • Maximum likelihood estimate of relative risk R  
  • Exact confidence interval for R  
  • F statistic and two-sided p-value  
  • Approximate two-sided power
- Accepts a simple 2×2 table as input
- Offers both printed output and an optional structured output (STATS)
- Suitable for epidemiological case-referent studies and matched-pair designs

📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/liddell

2. Add the folder to your MATLAB path:
      addpath('path_to_liddell')

3. Verify that MATLAB can find the function:
      which liddell

⚙️ Requirements
- MATLAB (any recent version)
- No additional toolboxes are strictly required (only core MATLAB functions are used)

📈 Usage
Basic example (drug vs placebo):

    % 2x2 matched-pairs table
    %               Drug
    %             +     -
    % Placebo +  101   59
    %         -  121   33
    x = [101 59; 121 33];

    % Run Liddell's exact test with default alpha = 0.05
    STATS = liddell(x);

When called without an output:

    liddell(x);

the function prints:

    Liddell's exact test
    Maximum likelihood estimate of relative risk (R) =  ...
    Exact 95% confidence interval = ... to ...
    F = ...  p-value (two-sided) = ...
    alpha = 0.0500  Zb = ...  Power (two-sided) = ...

Custom significance level:

    alpha = 0.01;
    STATS = liddell(x, alpha);

🔢 Inputs
liddell(x)
liddell(x, alpha)

- x     : 2×2 numeric matrix of nonnegative integers:
          x(1,1) = a (both positive),
          x(1,2) = b (positive in condition 1 only),
          x(2,1) = c (positive in condition 2 only),
          x(2,2) = d (both negative).

- alpha : Optional significance level in (0,1), default 0.05.
          Used to build the confidence interval and compute Z_alpha for
          power approximation.

📤 Outputs
If called with an output argument, LIDDELL returns a structure STATS with fields:

- STATS.Rhat       : maximum likelihood estimate of relative risk R
- STATS.CI         : [Rl Ru], exact (1 - alpha) confidence interval for R
- STATS.F          : F statistic for the hypothesis H0: R = 1
- STATS.pvalue     : two-sided p-value from the F distribution
- STATS.alpha      : significance level used
- STATS.Zb         : Z_beta used in power calculation
- STATS.power      : approximate two-sided power
- STATS.N          : total number of pairs (sum of all cells in x)
- STATS.table      : original 2×2 table x
- STATS.discordant : [r s], the sorted discordant counts

🧠 Interpretation
- Rhat > 1 suggests that the outcome is more frequent in one condition than the other (depending on how the table is defined).
- The exact confidence interval [Rl Ru] describes the plausible range of relative risk at level (1 - alpha).
- A small p-value (e.g., p < 0.05) indicates evidence against the null hypothesis that R = 1 (no difference in risk between the two conditions).
- The reported power is an approximation and should be interpreted cautiously; it is not a full design-stage power analysis.

📌 Notes
- Liddell's exact test is particularly recommended when the number of discordant pairs is small or when the assumptions behind asymptotic approximations
