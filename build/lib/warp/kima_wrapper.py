def filter_period(result_df, period, period_prefix='P', Np=4, cond='smaller', replace=False):
    """
    Filter the result DataFrame to include only rows where the period is within 10% of the specified period.

    Args:
        result_df (pd.DataFrame): The DataFrame containing the results.
        period (float): The period to filter by.
        period_prefix (str): The prefix for the period column name.
        Np (int): The number of periods to consider.
        cond (str): The condition for filtering ('smaller' or 'larger').
    """
    import numpy as np
    if not replace:
        result_df = result_df.copy()
    period_cols = [f'{period_prefix}{i+1}' for i in range(Np)]
    qtys = ['P', 'K', 'ecc', 'w', 'phi', 'mp']
    for it, row in result_df.iterrows():
        for idx in range(Np):
            if cond == 'smaller':
                if row[period_cols[idx]] < period:
                    for qty in qtys:
                        col_name = f'{qty}{idx+1}'
                        result_df.at[it, col_name] = np.nan
                    result_df.at[it, 'Np'] -= 1
            elif cond == 'larger':
                if row[period_cols[idx]] > period:
                    for qty in qtys:
                        col_name = f'{qty}{idx+1}'
                        result_df.at[it, col_name] = np.nan
                    result_df.at[it, 'Np'] -= 1

    return result_df
