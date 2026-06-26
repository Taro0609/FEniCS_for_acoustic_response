import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os

# 周波数リスト（10 Hz ～ 1000 Hz, 50 Hz刻み）
frequencies = np.arange(10, 1000 + 1, 50)

pressure_max = []
displacement_max = []

for idx, f in enumerate(frequencies):
    pres_filename = f"pressure_p0_{idx:06d}_p0_000000.vtu"
    disp_filename = f"displacement_p0_{idx:06d}_p0_000000.vtu"

    p_abs_max = np.nan
    u_max = np.nan

    # 音圧読み込み
    try:
        pressure_grid = pv.read(pres_filename)
        if "acoustic_pressure_real" not in pressure_grid.point_data:
            raise KeyError("Field 'acoustic_pressure_real' not found")
        p_vals = pressure_grid.point_data["acoustic_pressure_real"]
        p_abs_max = np.max(np.abs(p_vals))
        pressure_max.append(p_abs_max)
    except Exception as e:
        print(f"❌ Error reading pressure file {pres_filename}: {e}")
        pressure_max.append(np.nan)

    # 変位読み込み
    try:
        disp_grid = pv.read(disp_filename)
        if "f_real" not in disp_grid.point_data:
            raise KeyError("Field 'f_real' not found in displacement file")
        u_vals = disp_grid.point_data["f_real"]  # 実数ベクトル
        u_norms = np.linalg.norm(u_vals, axis=1)
        u_max = np.max(u_norms)
        displacement_max.append(u_max)
    except Exception as e:
        print(f"❌ Error reading displacement file {disp_filename}: {e}")
        displacement_max.append(np.nan)

    if not np.isnan(p_abs_max) and not np.isnan(u_max):
        print(f"✅ {f} Hz: Pressure max = {p_abs_max:.2e}, Displacement max = {u_max:.2e}")
    else:
        print(f"⚠️ Skipped {f} Hz due to read error.")

# 特定座標を定義（例：前耳骨など）
target_point = np.array([7.55542314, 17.2636752, 27.42938731])

def find_closest_point_index(mesh_points, target):
    distances = np.linalg.norm(mesh_points - target, axis=1)
    return np.argmin(distances)

pressure_at_target = []
displacement_at_target = []

for idx, f in enumerate(frequencies):
    pres_filename = f"pressure_p0_{idx:06d}_p0_000000.vtu"
    disp_filename = f"displacement_p0_{idx:06d}_p0_000000.vtu"

    p_val = np.nan
    u_val = np.nan

    try:
        # Pressure
        pressure_grid = pv.read(pres_filename)
        p_coords = pressure_grid.points
        p_idx = find_closest_point_index(p_coords, target_point)
        p_val = pressure_grid.point_data["acoustic_pressure_real"][p_idx]
        pressure_at_target.append(p_val)
    except Exception as e:
        print(f"❌ Error reading pressure: {e}")
        pressure_at_target.append(np.nan)

    try:
        # Displacement
        disp_grid = pv.read(disp_filename)
        u_coords = disp_grid.points
        u_idx = find_closest_point_index(u_coords, target_point)
        u_vector = disp_grid.point_data["f_real"][u_idx]
        u_norm = np.linalg.norm(u_vector)
        displacement_at_target.append(u_norm)
    except Exception as e:
        print(f"❌ Error reading displacement: {e}")
        displacement_at_target.append(np.nan)

    print(f"✅ {f} Hz: P={p_val:.2e}, U={u_norm:.2e}")




# プロット
fig, ax1 = plt.subplots()

color1 = 'tab:red'
ax1.set_xlabel('Frequency [Hz]')
ax1.set_ylabel('Max Acoustic Pressure [Pa]', color=color1)
ax1.plot(frequencies, pressure_max, 'o-', color=color1)
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()
color2 = 'tab:blue'
ax2.set_ylabel('Max Displacement [m]', color=color2)
ax2.plot(frequencies, displacement_max, 's--', color=color2)
ax2.tick_params(axis='y', labelcolor=color2)

plt.title("Frequency Response: Pressure and Displacement")
fig.tight_layout()
plt.grid(True)
plt.savefig("frequency_response_plot.png", dpi=300)
plt.show()



# --- プロット ---
fig, ax1 = plt.subplots()

ax1.set_xlabel("Frequency [Hz]")
ax1.set_ylabel("Acoustic Pressure [Pa]", color='tab:red')
ax1.plot(frequencies, pressure_at_target, 'o-', color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')

ax2 = ax1.twinx()
ax2.set_ylabel("Displacement [m]", color='tab:blue')
ax2.plot(frequencies, displacement_at_target, 's--', color='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')

plt.title("Pressure and Displacement at Specific Point")
plt.grid(True)
plt.tight_layout()
plt.savefig("target_point_response.png", dpi=300)
plt.show()





