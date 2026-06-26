import pyvista as pv
import numpy as np

# メッシュ読み込みと表面抽出（三角形の面だけにする）
mesh = pv.read("output.msh")
surface = mesh.extract_surface()
picked_faces = set()

print("✅ 矩形選択（through=True）を複数回行えます")
print("   → Shift + R → ドラッグ → Qキーでウィンドウを閉じる")

while True:
    plotter = pv.Plotter()

    picked_ids = []

    # 選択時のコールバック関数
    def callback(picked):
        if picked.n_cells == 0:
            print("⚠️ 面が選択されていません")
            return
        ids = picked["vtkOriginalCellIds"] if "vtkOriginalCellIds" in picked.array_names else np.arange(picked.n_cells)
        picked_faces.update(ids)
        picked_ids.append(True)
        print(f"✅ 今回選択された面数: {len(ids)}（累計: {len(picked_faces)}）")

    # メッシュ描画と矩形選択（through=True）
    plotter.add_mesh(surface, show_edges=True, color="lightgray")
    plotter.enable_cell_picking(
        callback=callback,
        through=True,        # ←★ ここを変更
        show=True,
        style="surface"
    )
    plotter.show()

    if not picked_ids:
        print("⚠️ 面が選択されませんでした")

    ans = input("続けて選択しますか？ (y/n): ")
    if ans.lower() != "y":
        break

# ✅ 抽出と保存処理
if picked_faces:
    extracted = surface.extract_cells(list(picked_faces))
    extracted.save("selected_vibration_faces.vtu")
    with open("selected_vibration_face_ids.txt", "w") as f:
        for cid in picked_faces:
            f.write(f"{cid}\n")
    print(f"✅ 最終的な選択面数: {len(picked_faces)}")
    print("保存しました: selected_vibration_faces.vtu / selected_vibration_face_ids.txt")
else:
    print("⚠️ 面は選択されませんでした")
