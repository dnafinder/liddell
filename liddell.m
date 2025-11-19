function STATS = liddell(x, varargin)
%LIDDELL Liddell's exact test for matched pairs (paired proportions).
%
%   liddell(x)
%   STATS = liddell(x)
%   liddell(x, alpha)
%   STATS = liddell(x, alpha)
%
%   Description:
%       LIDDELL performs Liddell's exact test on matched pairs, which is an
%       exact alternative to McNemar's chi-square test for paired binary
%       outcomes. It is typically used when the same subjects are observed
%       under two conditions (e.g., drug vs placebo) and each observation
%       is classified as "positive" or "negative".
%
%       The input x is a 2x2 table:
%
%                   Condition 2
%                 +            -
%           +   x(1,1)      x(1,2)
%   Cond 1
%           -   x(2,1)      x(2,2)
%
%       The test focuses on the discordant pairs:
%           b = x(1,2) (positive in condition 1, negative in condition 2)
%           c = x(2,1) (negative in condition 1, positive in condition 2)
%
%       Liddell's test provides:
%           - Maximum likelihood estimate of relative risk (R)
%           - Exact confidence interval for R
%           - F statistic with p-value (two-sided)
%           - Approximate power of the test
%
%   Inputs:
%       x     - 2x2 numeric matrix of nonnegative integers, real and finite.
%       alpha - (Optional) Significance level for confidence interval and
%               computation of Z_α (for power). Scalar in (0,1).
%               Default: 0.05
%
%   Outputs:
%       STATS - (Optional) structure with fields:
%                  STATS.Rhat      : maximum likelihood estimate of R
%                  STATS.CI        : [Rl Ru] exact (1 - alpha) CI for R
%                  STATS.F         : F statistic for H0: R = 1
%                  STATS.pvalue    : two-sided p-value
%                  STATS.alpha     : significance level
%                  STATS.Zb        : Z_beta (used in power computation)
%                  STATS.power     : approximate two-sided power
%                  STATS.N         : total number of pairs
%                  STATS.table     : input 2x2 table x
%                  STATS.discordant: [r s] sorted discordant counts
%
%       If no output is requested, LIDDELL prints a formatted summary to
%       the Command Window.
%
%   Example:
%       % Drug vs placebo example:
%       %                  Drug
%       %             +        -
%       %       +   101      59
%       % Placebo
%       %       -   121      33
%       x = [101 59; 121 33];
%       STATS = liddell(x);
%
%   GitHub repository:
%       https://github.com/dnafinder/liddell
%
%   Citation:
%       Cardillo G. (2025). liddell: Liddell's exact test for matched pairs
%       in MATLAB. Available at:
%       https://github.com/dnafinder/liddell
%
%       Original statistical method:
%       Liddell F.D.K. (1983). Simplified exact analysis of case-referent
%       studies; matched pairs; dichotomous exposure. Journal of
%       Epidemiology and Community Health, 37, 82–84.
%
%   License:
%       This function is distributed under the terms specified in the
%       LICENSE file of the liddell repository.
%
%   Author:
%       Giuseppe Cardillo
%       giuseppe.cardillo.75@gmail.com
%
%   Created:
%       2008-01-01 (original concept)
%
%   Updated:
%       2025-11-19 (refactored and documented version)
%
%   Version:
%       1.1.0

% -----------------------------
% Input parsing and validation
% -----------------------------
p = inputParser;

% 2x2 table of nonnegative integers
addRequired(p, 'x', @(v) validateattributes(v, ...
    {'numeric'}, ...
    {'real', 'finite', 'integer', 'nonnegative', 'nonnan', 'size', [2 2]}));

% Significance level alpha
addOptional(p, 'alpha', 0.05, @(v) ...
    validateattributes(v, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnan', '>', 0, '<', 1}));

parse(p, x, varargin{:});
x     = p.Results.x;
alpha = p.Results.alpha;

clear p;

% -----------------------------
% Extract discordant pairs
% -----------------------------
% Discordant pairs are the off-diagonal elements of x.
% We flip the matrix left-right and take the diagonal to obtain:
%   [b; c] where
%       b = x(1,2)
%       c = x(2,1)
ob = diag(fliplr(x));

% Ensure that the first element is the larger discordant count
% so that r >= s (by flipping if necessary).
if ob(1) < ob(2)
    ob = flipud(ob);
end

r = ob(1);
s = ob(2);

% Warn if there are no discordant pairs
if r == 0 && s == 0
    warning('liddell:NoDiscordantPairs', ...
        'There are no discordant pairs (off-diagonal cells are zero). Liddell''s test is not informative.');
end

% -----------------------------
% Relative risk and confidence interval
% -----------------------------
% Maximum likelihood estimate of relative risk R (Mantel & Haenszel, 1959):
%   R = r / s
% Based only on discordant pairs.
if s == 0
    Rhat = Inf;
else
    Rhat = r / s;
end

% Exact (1 - alpha) confidence interval for R using F distribution.
P  = 1 - alpha / 2;
% Upper bound Ru = (r+1)/s * F(P; 2*(r+1), 2*s)
% Lower bound Rl = r/(s+1) / F(P; 2*(s+1), 2*r)
if r > 0 && s > 0
    Ru = (r + 1) / s * finv(P, 2*(r + 1), 2*s);
    Rl = r / (s + 1) / finv(P, 2*(s + 1), 2*r);
else
    % When r == 0 or s == 0, the usual CI formulas may not be valid.
    % Here we return NaNs and issue a warning.
    Ru = NaN;
    Rl = NaN;
    warning('liddell:DegenerateCI', ...
        ['One of the discordant counts is zero. Standard exact CI ', ...
         'formulas for R may not be well-defined.']);
end

% -----------------------------
% Test of H0: R = 1
% -----------------------------
% F statistic:
%   F = r / (s + 1)
Fstat = r / (s + 1);

% Two-sided p-value from F distribution:
pvalue = 2 * (1 - fcdf(Fstat, 2*(s + 1), 2*r));

% -----------------------------
% Approximate power calculation
% -----------------------------
% Compute Za for given alpha:
Za = abs(-sqrt(2) * erfcinv(alpha));

% Total sample size:
N = sum(x(:));

% Proportion of minimum discordant cell:
p = min(ob ./ N);

% Ratio of discordant counts (larger / smaller), avoiding division by zero.
if min(ob) == 0
    pp = Inf;
else
    pp = max(ob(1)/ob(2), ob(2)/ob(1));
end

% If pp is infinite (one discordant cell is zero), power formula is not
% meaningful. We handle this separately.
if ~isfinite(pp)
    Zb   = NaN;
    pwr  = NaN;
    warning('liddell:PowerUndefined', ...
        'Power calculation is undefined when one discordant cell is zero.');
else
    num   = abs(sqrt(N * p * (pp - 1)^2) - sqrt(Za^2 * (pp + 1)));
    denom = sqrt(pp + 1 - p * (pp - 1)^2);
    Zb    = num / denom;
    % Two-sided power approximation using the normal distribution:
    pwr   = (1 - 0.5 * erfc(-Zb / sqrt(2))) * 2;
end

% -----------------------------
% Display results (if no output)
% -----------------------------
if nargout == 0
    disp('Liddell''s exact test');
    fprintf('Maximum likelihood estimate of relative risk (R) =  %0.4f\n', Rhat);
    fprintf('Exact %i%% confidence interval = %0.4f to %0.4f\n', ...
        round((1 - alpha) * 100), Rl, Ru);
    fprintf('F = %0.4f  p-value (two-sided) = %0.8f\n', Fstat, pvalue);
    fprintf('alpha = %0.4f  Zb = %0.4f  Power (two-sided) = %0.4f\n', ...
        alpha, Zb, pwr);
end

% -----------------------------
% Build output structure (if requested)
% -----------------------------
if nargout > 0
    STATS.Rhat       = Rhat;
    STATS.CI         = [Rl, Ru];
    STATS.F          = Fstat;
    STATS.pvalue     = pvalue;
    STATS.alpha      = alpha;
    STATS.Zb         = Zb;
    STATS.power      = pwr;
    STATS.N          = N;
    STATS.table      = x;
    STATS.discordant = [r, s];
end

end
