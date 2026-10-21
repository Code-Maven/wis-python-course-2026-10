
import pandas as pd

df = pd.read_csv("people.csv")
print(df)
df["BMI"] = 10000*df["Weight"]/df["Height"]**2
print(df)
df["Text"] = pd.cut(
    df["BMI"],
    bins=[0, 18.5, 25, 30, float("inf")],
    labels=["Underweight", "Healthy", "Overweight", "Obesity"],
    right=False
)
print(df)
