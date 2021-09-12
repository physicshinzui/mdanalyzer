import matplotlib.pyplot as plt 
import pandas as pd

def main():
    df1 = pd.read_csv('rmsd_mda.csv', index_col = 1)
    df2 = pd.read_csv('rmsf_mda.csv', index_col = 0)

    fig, axes = plt.subplots(nrows = 1, ncols = 2)
    df1.iloc[:,1].plot(ax=axes[0])
    df2.plot(ax=axes[1])
    
    plt.tight_layout()
    plt.show()

main()