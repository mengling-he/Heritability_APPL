import warnings
import pandas as pd
from statsmodels.regression.mixed_linear_model import MixedLM
### DO not use this code, use function in R!!!!!!!!!!!!
def estimate_h2_per_time(df, Traits, time_col='Day', group_col='Genotype'):
    """
    Estimate broad-sense heritability (H²) for each metric at each day.

    Parameters:
    - df: pandas DataFrame with your data
    - metrics: list of column names to estimate H² for
    - time_col: column name for day (numeric)
    - group_col: column name for genotype/group

    Returns:
    - h2_df: DataFrame with columns ['Date', 'Metric', 'H2', 'var_genotype', 'var_residual']
    """
    h2_results = []
    dates = sorted(df[time_col].unique())

    for d in dates:
        df_date = df[df[time_col] == d]
        
        for m in Traits:
            if m not in df_date.columns or df_date[m].notna().sum() < 2:
                warnings.warn(f"Not enough data to estimate H2 for metric '{m}' on day {d}. Skipping.")
                H2 = var_genotype = var_residual = None
            else:
                try:
                    # Fit the mixed model
                    md = MixedLM(df_date[m], pd.DataFrame({'Intercept': 1}, index=df_date.index), 
                                 groups=df_date[group_col])
                    mdf = md.fit(reml=True)
                    var_genotype = mdf.cov_re.iloc[0,0]
                    var_residual = mdf.scale
                    H2 = var_genotype / (var_genotype + var_residual)
                except Exception as e:
                    warnings.warn(f"Model failed for Trait '{m}' on day {d}: {e}")
                    H2 = var_genotype = var_residual = None
                    
            h2_results.append({
                time_col: d,
                'Trait': m,
                'H2': H2,
                'var_genotype': var_genotype,
                'var_residual': var_residual
            })
    
    h2_df = pd.DataFrame(h2_results)
    return h2_df