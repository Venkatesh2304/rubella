import pandas as pd 
#cbr
df=pd.read_excel("data/mortality.xlsx")
df["year"] = df["YEAR"].apply(lambda x : x.strftime("%Y"))
df[["year","MR"]].to_excel("mortality.xlsx",index=False)
