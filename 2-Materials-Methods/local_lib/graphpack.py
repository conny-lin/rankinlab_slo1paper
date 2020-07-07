import pandas as pd
import numpy as np


# save to excel
def save_excel_graphdata(df_output, iv, MEASURES, savepath):
    with pd.ExcelWriter(savepath) as writer:
        for msr in MEASURES:
            # get data
            df = df_output[[iv, msr]].copy()
            # create summary
            t = df.groupby(iv).agg(['count','mean','sem'])
            # save to excel
            t.to_excel(writer, sheet_name=msr)
