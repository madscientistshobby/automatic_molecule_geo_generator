import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# ========= 사용자 설정 =========
file_path = 'filtered_smiles_dataset_ph_x.xlsx'
x_col = 'sigma_para_exp'
targets = ['Est', 'E_singlet', 'E_triplet']
output_dir = 'regression_plots_multi'
os.makedirs(output_dir, exist_ok=True)
# ==============================

# 데이터 로드
df = pd.read_excel(file_path)

if x_col not in df.columns:
    raise ValueError(f"{x_col} 열이 존재하지 않습니다.")

# ----- 모델 정의 -----
def linear(x, a, b):
    return a * x + b

def poly2(x, a, b, c):
    return a * x**2 + b * x + c

def poly3(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

def exp_fn(x, a, b, c):
    # a * exp(bx) + c
    return a * np.exp(b * x) + c

def log_fn(x, a, b):
    # a * ln(x) + b (x>0에서만)
    return a * np.log(x) + b

def sigmoid(x, a, b, c, d):
    # a / (1 + exp(-b(x-c))) + d
    return a / (1.0 + np.exp(-b * (x - c))) + d

# (이름, 함수, 도메인 조건)
model_specs = [
    ("Linear",       linear,  None),
    ("Poly2",        poly2,   None),
    ("Poly3",        poly3,   None),
    ("Exponential",  exp_fn,  None),
    ("Logarithmic",  log_fn,  lambda x: x > 0),
    ("Sigmoid",      sigmoid, None),
]

def calc_r2(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1 - ss_res / ss_tot

# ----- 각 target에 대해 피팅 -----
for col in targets:
    if col not in df.columns:
        print(f"⚠️ {col} 열이 존재하지 않습니다. 건너뜁니다.")
        continue

    sub = df[[x_col, col]].dropna()
    x_all = sub[x_col].values.astype(float)
    y_all = sub[col].values.astype(float)

    if len(sub) < 5:
        print(f"⚠️ {col}: 데이터 포인트가 너무 적어서 스킵합니다. (n={len(sub)})")
        continue

    results = []  # (name, popt, r2, x_fit, y_fit, domain_mask_info)

    for name, func, domain_cond in model_specs:
        # 도메인 조건 적용 (예: log)
        if domain_cond is not None:
            mask = domain_cond(x_all)
            if not np.any(mask):
                continue
            x = x_all[mask]
            y = y_all[mask]
        else:
            x = x_all
            y = y_all

        if len(x) < 5:
            continue

        try:
            popt, _ = curve_fit(func, x, y, maxfev=10000)
            y_pred = func(x, *popt)
            r2 = calc_r2(y, y_pred)

            # 해당 모델의 유효 구간에서만 곡선 생성
            x_fit = np.linspace(x.min(), x.max(), 300)
            y_fit = func(x_fit, *popt)

            results.append((name, popt, r2, x_fit, y_fit))
        except Exception:
            continue

    if not results:
        print(f"❌ {col}: 어떤 모델도 피팅에 실패했습니다.")
        continue

    # R² 기준 정렬
    results.sort(key=lambda t: t[2], reverse=True)

    print(f"\n📊 {col} vs {x_col}")
    for name, popt, r2, _, _ in results:
        params_str = ", ".join(f"{p:.3g}" for p in popt)
        print(f"- {name:12s}: R² = {r2:.4f}, params = [{params_str}]")

    best_name, _, best_r2, _, _ = results[0]
    print(f"👉 Best model for {col}: {best_name} (R² = {best_r2:.4f})")

    # ✅ 모든 모델에 대해 개별 플롯 저장
    for name, popt, r2, x_fit, y_fit in results:
        plt.figure(figsize=(6, 4))
        # 원 데이터 전체 산점도
        plt.scatter(x_all, y_all, label='Data', alpha=0.7)
        # 해당 모델 피팅 곡선
        plt.plot(x_fit, y_fit,
                 label=f'{name} fit\nR²={r2:.3f}',
                 linewidth=2)
        plt.xlabel(f'Hammett parameter ({x_col})')
        plt.ylabel(col)
        plt.title(f'{x_col} vs {col} - {name}')
        plt.legend()
        plt.tight_layout()

        fname = f'{x_col}_vs_{col}_{name}.png'
        save_path = os.path.join(output_dir, fname)
        plt.savefig(save_path, dpi=300)
        plt.close()

print(f"\n✅ 모든 타겟/모델 피팅 및 플롯 저장 완료: '{output_dir}/' 폴더 확인")
