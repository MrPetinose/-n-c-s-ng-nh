import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def visualize_clustering(csv_file='Auto_Result.csv'):
    try:
        df = pd.read_csv(csv_file)
        if df.empty:
            print("File CSV rong!")
            return

        clusters = df['Cluster_ID'].unique()
        palette = sns.color_palette("bright", len(clusters))
        
        peaks = pd.DataFrame()
        if 'Is_Peak' in df.columns:
            peaks = df[df['Is_Peak'] == 'YES']

        sns.set_style("whitegrid")
        ##
        plt.figure("KET QUA PHAN CUM", figsize=(10, 8))
        
        sns.scatterplot(
            data=df, x='X', y='Y', hue='Cluster_ID', 
            palette=palette, c = 'blue',s=60, alpha=0.7, edgecolor='k', legend='full'
        )
        

        if not peaks.empty:
            plt.scatter(peaks['X'], peaks['Y'], c='black', s=100, marker='*', label='Peaks', zorder=10)

        plt.title('Ket qua Phan cum (Spatial Domain)', fontsize=16, fontweight='bold')
        plt.xlabel('Toa do X')
        plt.ylabel('Toa do Y')
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title="Cluster ID")
        plt.tight_layout()

        ##
        plt.figure("BIEU DO QUYET DINH (Rho-Delta)", figsize=(10, 8))

        plt.scatter(df['Rho'], df['Delta'], c='gray', s=60, alpha=0.5, label='Normal Points')

        if not peaks.empty:
            plt.scatter(peaks['Rho'], peaks['Delta'], c='red', s=60, edgecolors='black', label='Density Peaks', zorder=10)
            

        plt.title('Decision Graph (Bieu do Quyet dinh)', fontsize=16, fontweight='bold')
        plt.xlabel('Mat do (Rho)', fontsize=12)
        plt.ylabel('Khoang cach (Delta)', fontsize=12)
        plt.legend()
        
        # Thêm chú thích
        # plt.text(0.95, 0.05, 'Peaks o goc tren-phai', 
        #          horizontalalignment='right', verticalalignment='bottom', 
        #          transform=plt.gca().transAxes, fontsize=10, color='blue', style='italic')
        # plt.tight_layout()

        print(f"Dang hien thi {len(df)} diem:")
        plt.show()

    except FileNotFoundError:
        print(f"LOI: Khong tim thay file '{csv_file}'")
    except Exception as e:
        print(f"Co loi xay ra: {e}")

if __name__ == "__main__":
    visualize_clustering()