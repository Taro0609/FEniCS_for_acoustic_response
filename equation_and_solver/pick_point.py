import pyvista as pv

# メッシュの読み込み
mesh = pv.read("output.msh")

# ピック時に呼ばれる関数（引数2つ！）
def on_pick(point, picker):
    print("Picked Point Coordinates:", point)

# プロッター生成
plotter = pv.Plotter()
plotter.add_mesh(mesh, show_edges=True, color="lightgray")

# ピック有効化（use_picker=True が推奨）
plotter.enable_point_picking(callback=on_pick, use_picker=True, show_message=True)

# 表示
plotter.show()
