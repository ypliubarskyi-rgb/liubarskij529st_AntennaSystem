import numpy as np
import matplotlib.pyplot as plt
import math


f = [400, 420, 440, 460, 480, 520, 540, 560, 580, 600, 620, 640, 630, 610, 600, 590, 580, 570, 560, 550, 540, 530, 520, 510, 500]
n = [6, 7, 8, 9, 12, 11, 10, 9, 8, 7, 6, 5, 12, 11, 10, 9, 8, 7, 6, 6, 7, 7, 8, 8, 9]


variant = 7

print(f"=== ВАРІАНТ {variant} ===")
print(f"Частота f = {f[variant - 1]:.0f} (МГц)")
print(f"Число елементів антени N = {n[variant - 1]:.0f}")

lambd = 299792458 / (f[variant - 1] * 10 ** 6)
print(f"λ = {lambd:.3f} (м)")

d = 0.25 * lambd
print(f"dcp = {d:.3f} (м)")

k = (2 * np.pi) / lambd
print(f"k = {k:.3f} (рад/м)")


F1E = []
FC = []
FE = []
steps = []

max_x_FC, max_y_FC = [], []
max_x_FE, max_y_FE = [], []
min_x_FC, min_y_FC = [], []
min_x_FE, min_y_FE = [], []


for teta in np.arange(0.01, np.pi / 2, 0.0001):
    mn1 = abs(np.cos((np.pi / 2) * np.sin(teta)) / np.cos(teta))


    num = np.sin((n[variant - 1] * k * d * (1 - np.cos(teta))) / 2)
    den = n[variant - 1] * np.sin((k * d * (1 - np.cos(teta))) / 2)
    mn2 = abs(num / den)

    mn3 = mn1 * mn2

    F1E.append(mn1)
    FC.append(mn2)
    FE.append(mn3)
    deg_val = math.degrees(teta)
    steps.append(deg_val)

steps = np.array(steps)
FC = np.array(FC)
FE = np.array(FE)
F1E = np.array(F1E)


idx_0707_H = np.argmin(np.abs(FC - 0.707))
idx_0707_E = np.argmin(np.abs(FE - 0.707))

theta_0707_H = steps[idx_0707_H]
theta_0707_E = steps[idx_0707_E]

fS1 = FC[idx_0707_H]
fS2 = FE[idx_0707_E]

SGP1 = 2 * theta_0707_H
SGP2 = 2 * theta_0707_E


for i in range(1, len(FC) - 1):
    if FC[i] > FC[i - 1] and FC[i] > FC[i + 1] and FC[i] < 0.99:
        max_x_FC.append(steps[i])
        max_y_FC.append(FC[i])
    if FE[i] > FE[i - 1] and FE[i] > FE[i + 1] and FE[i] < 0.99:
        max_x_FE.append(steps[i])
        max_y_FE.append(FE[i])

    if FC[i] < FC[i - 1] and FC[i] < FC[i + 1] and FC[i] < 0.05:
        min_x_FC.append(steps[i])
        min_y_FC.append(FC[i])
    if FE[i] < FE[i - 1] and FE[i] < FE[i + 1] and FE[i] < 0.05:
        min_x_FE.append(steps[i])
        min_y_FE.append(FE[i])

print("\n=== РЕЗУЛЬТАТИ АНАЛІЗУ ===")
print("Ширина головної пелюстки в площині H (2θ_0.7) = " + str(round(SGP1, 2)) + '°')
print("Ширина головної пелюстки в площині E (2θ_0.7) = " + str(round(SGP2, 2)) + '°')


len_value = [len(max_x_FC), len(max_x_FE), len(min_x_FC), len(min_x_FE)]
print("\nТабл. 1 - Аналіз ДС Директорної антени в площині Н та Е (Варіант 7)")
print("-----------------------------------------")
print("| № |θminH|θminE|θmaxH|FH(θ)|θmaxE|FE(θ)|")
print("-----------------------------------------")
for i in range(0, max(len_value)):
    v1 = f" {i + 1:.0f} "
    v2 = f"{min_x_FC[i]:.2f}" if i < len(min_x_FC) else '  -  '
    v3 = f"{min_x_FE[i]:.2f}" if i < len(min_x_FE) else '  -  '
    v4 = f"{max_x_FC[i]:.2f}" if i < len(max_x_FC) else '  -  '
    h4 = f"{max_y_FC[i]:.3f}" if i < len(max_y_FC) else '  -  '
    v5 = f"{max_x_FE[i]:.2f}" if i < len(max_x_FE) else '  -  '
    h5 = f"{max_y_FE[i]:.3f}" if i < len(max_y_FE) else '  -  '
    print(f"|{v1}|{v2}|{v3}|{v4}|{h4}|{v5}|{h5}|")
print("-----------------------------------------")


fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
ax.plot(steps, F1E, linewidth=1, linestyle=':', color='gray', label="$F_{1e}(\\theta)$")
ax.plot(steps, FC, linewidth=1.8, color='blue', label="$F_{H}(\\theta)$")
ax.plot(steps, FE, linewidth=1.8, color='red', label="$F_{E}(\\theta)$")


ax.axhline(0.707, color='green', linestyle='--', linewidth=1, label="Рівень 0.707 (-3 дБ)")


ax.axvline(theta_0707_E, color='red', linestyle=':', linewidth=1.5)
ax.plot(theta_0707_E, fS2, 'ro', markersize=6, zorder=5)

ax.axvline(theta_0707_H, color='blue', linestyle=':', linewidth=1.5)
ax.plot(theta_0707_H, fS1, 'bo', markersize=6, zorder=5)


ax.plot(max_x_FC, max_y_FC, "bo", markersize=6, label="Макс. $F_H$")
ax.plot(max_x_FE, max_y_FE, "go", markersize=6, label="Макс. $F_E$")


ax.plot(min_x_FC, min_y_FC, "ro", markersize=6, label="Нулі випромінювання")

ax.set_xlabel('Кут $\\theta$, град', fontsize=10)
ax.set_ylabel('Нормована амплітуда $F(\\theta)$', fontsize=10)
ax.set_xlim(0, 90)
ax.set_ylim(-0.02, 1.05)
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax.set_title("Рисунок 2 – Діаграма спрямованості антени у площині Е та Н (Варіант 7)", fontsize=11, fontweight='bold')
ax.legend(loc="upper right", fontsize=8)
plt.tight_layout()
fig.savefig("Figure_2_Rectangular_Var7.png", dpi=300)
plt.show()


fig2 = plt.figure(figsize=(7, 7), dpi=300)
ax_polar = fig2.add_subplot(111, projection='polar')

theta_rad = np.radians(steps)

ax_polar.plot(theta_rad, FC, linewidth=1.8, color="blue", label="$F_H(\\theta)$")
ax_polar.plot(theta_rad, FE, linewidth=1.8, color="red", label="$F_E(\\theta)$")


ax_polar.scatter(np.radians(max_x_FC), max_y_FC, color="black", s=30, label="Макс. FH", zorder=5)
ax_polar.scatter(np.radians(max_x_FE), max_y_FE, color="gray", s=30, label="Макс. FE", zorder=5)
ax_polar.scatter(np.radians(min_x_FC), min_y_FC, color="red", s=30, label="Нулі", zorder=5)

ax_polar.set_theta_zero_location("N")
ax_polar.set_theta_direction(-1)
ax_polar.set_thetamin(0)
ax_polar.set_thetamax(90)

ax_polar.set_title("Рисунок 3 – Діаграма спрямованості антени\nу площині Е та Н у полярному вигляді (Варіант 7)",
                   fontsize=10, fontweight='bold', pad=20)
ax_polar.set_rlabel_position(90)
ax_polar.grid(True, linestyle='--', linewidth=0.5)
ax_polar.legend(loc="lower left", bbox_to_anchor=(-0.1, -0.15), fontsize=8)

plt.tight_layout()
fig2.savefig("Figure_3_Polar_Var7.png", dpi=300)
plt.show()