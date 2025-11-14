import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

# 엑셀 파일 불러오기
file_path = 'filtered_smiles_dataset_ph_x.xlsx'  # 파일 이름 수정
df = pd.read_excel(file_path)

# 결과 저장 폴더 생성
output_dir = 'regression_plots'
os.makedirs(output_dir, exist_ok=True)

# 분석 대상 열 리스트
targets = ['Est', 'E_singlet', 'E_triplet']

for col in targets:
    if col not in df.columns:
        print(f"⚠️ {col} 열이 존재하지 않습니다. 건너뜁니다.")
        continue
    
    x = df['sigma_para_exp']
    y = df[col]
    
    # 선형회귀
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    
    # 회귀 결과 출력
    print(f"\n📈 {col} vs Hammett")
    print(f"회귀식: {col} = {slope:.5f} * Hammett + {intercept:.5f}")
    print(f"R² = {r_value**2:.4f}, p = {p_value:.3e}")
    
    # 그래프 생성
    plt.figure(figsize=(6,4))
    plt.scatter(x, y, color='blue', label='Data')
    plt.plot(x, slope*x + intercept, color='red',
             label=f'y={slope:.3f}x+{intercept:.3f}\nR²={r_value**2:.2f}')
    plt.xlabel('Hammett parameter (σ)')
    plt.ylabel(col)
    plt.title(f'Hammett vs {col}')
    plt.legend()
    plt.tight_layout()
    
    # 자동 저장
    save_path = os.path.join(output_dir, f'Hammett_vs_{col}.png')
    plt.savefig(save_path, dpi=300)
    plt.close()

print(f"\n✅ 모든 그래프가 '{output_dir}/' 폴더에 저장되었습니다.")
