import numpy as np
import matplotlib.pyplot as plt

lam = 3.1
ap = 11.0
bp = 11.0

theta_deg = np.linspace(0, 90, 901)
theta_rad = np.radians(theta_deg)

u_E = (np.pi * bp / lam) * np.sin(theta_rad)
u_E[u_E == 0] = 1e-10
F_CE = np.sin(u_E) / u_E
F_1E = (1 + np.cos(theta_rad)) / 2
F_E = np.abs(F_CE * F_1E)
F_E_norm = F_E / np.max(F_E)

denom_H = 1 - (2 * ap * np.sin(theta_rad) / lam)**2
denom_H[np.abs(denom_H) < 1e-5] = 1e-5
F_CH = np.cos((np.pi * ap / lam) * np.sin(theta_rad)) / denom_H
F_1H = (1 + np.cos(theta_rad)) / 2
F_H = np.abs(F_CH * F_1H)
F_H_norm = F_H / np.max(F_H)

F_1e = F_1E

idx_07_E = np.where(F_E_norm <= 0.707)[0][0]
idx_07_H = np.where(F_H_norm <= 0.707)[0][0]

theta_null_deg = np.degrees(np.arcsin(lam / bp))
idx_null = np.argmin(np.abs(theta_deg - theta_null_deg))

theta_max1_deg = np.degrees(np.arcsin(1.5 * lam / bp))
idx_max1_E = np.argmin(np.abs(theta_deg - theta_max1_deg))
idx_max1_H = np.argmin(np.abs(theta_deg - theta_max1_deg))

plt.figure(figsize=(10, 6))

plt.plot(theta_deg, F_1e, label='$F_{1e}(\\theta)$', color='grey', linestyle=':', linewidth=1)
plt.plot(theta_deg, F_H_norm, label='$F_H(\\theta)$', color='blue', linewidth=2)
plt.plot(theta_deg, F_E_norm, label='$F_E(\\theta)$', color='red', linewidth=2)

plt.axhline(0.707, color='green', linestyle='--', label='Рівень 0.707 (-3 дБ)')

plt.plot(theta_deg[idx_07_H], F_H_norm[idx_07_H], marker='o', color='blue')
plt.plot(theta_deg[idx_07_E], F_E_norm[idx_07_E], marker='o', color='red')
plt.axvline(theta_deg[idx_07_H], color='blue', linestyle=':')
plt.axvline(theta_deg[idx_07_E], color='red', linestyle=':')

plt.plot(theta_deg[idx_null], F_E_norm[idx_null], marker='o', color='red', label='Нулі випромінювання')
plt.plot(theta_deg[idx_null], F_H_norm[idx_null], marker='o', color='red')

plt.plot(theta_deg[idx_max1_H], F_H_norm[idx_max1_H], marker='o', color='blue', label='Макс. $F_H$')
plt.plot(theta_deg[idx_max1_E], F_E_norm[idx_max1_E], marker='o', color='green', label='Макс. $F_E$')

plt.title('Рисунок 2 – Діаграма спрямованості антени у площині Е та Н (Варіант 7)', fontsize=12, fontweight='bold')
plt.xlabel('Кут $\\theta$, град', fontsize=10)
plt.ylabel('Нормована амплітуда $F(\\theta)$', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=9)
plt.xlim(0, 90)
plt.ylim(-0.02, 1.05)

plt.tight_layout()
plt.savefig('horn_antenna_var7', dpi=300)
plt.show()