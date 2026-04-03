import os, json, time
from pathlib import Path
from PIL import Image
import numpy as np
import dearpygui.dearpygui as dpg

REPO_ROOT = "/home/nate/research/Nates-OSPM"
COCKPIT = Path(REPO_ROOT) / "Cockpit"
SHARED  = COCKPIT / "shared"

REFRESH = 10.0

current_galaxy = None
last_update = 0
file_mtimes = {}
texture_meta = {}
image_widgets = set()
zoom_scale = 1.0
widget_dims = {}  # widget_tag -> (base_w, base_h) at zoom 1.0

# =========================
# UTIL
# =========================
def load_latest_galaxy():
    try:
        return (SHARED / "latest_galaxy").read_text().strip()
    except:
        return None

def load_json(path):
    try:
        return json.load(open(path))
    except:
        return {}

def describe_plot(name):
    name = name.lower()
    if "chi2" in name:
        return "χ² surface — lower is better, look for tight minima"
    if "geometry" in name:
        return "parameter projection — clustering = viable region"
    if "hist" in name:
        return "distribution of accepted models"
    if "gp" in name:
        return "GP interpolation — smooth likelihood estimate"
    return "model diagnostic"

def load_image_texture(tag, path):
    if not os.path.exists(path):
        return False

    img = Image.open(path).convert("RGBA")
    w, h = img.size

    arr = np.asarray(img, dtype=np.float32) / 255.0
    flat = arr.flatten().tolist()

    if dpg.does_item_exist(tag):
        if texture_meta.get(tag) != (w, h):
            dpg.delete_item(tag)
            dpg.add_static_texture(w, h, flat, tag=tag, parent="tex_reg")
            texture_meta[tag] = (w, h)
        else:
            dpg.set_value(tag, flat)
    else:
        dpg.add_static_texture(w, h, flat, tag=tag, parent="tex_reg")
        texture_meta[tag] = (w, h)

    return True

# =========================
# POPUP
# =========================
def show_full_image(sender, app_data, user_data):
    path = user_data
    if not path or not os.path.exists(path):
        return

    tag = f"popup_{path}"

    if dpg.does_item_exist(tag):
        dpg.delete_item(tag)

    img = Image.open(path).convert("RGBA")
    w, h = img.size

    arr = np.asarray(img, dtype=np.float32) / 255.0
    flat = arr.flatten().tolist()

    tex_tag = f"{tag}_tex"

    if dpg.does_item_exist(tex_tag):
        dpg.delete_item(tex_tag)

    dpg.add_static_texture(w, h, flat, tag=tex_tag, parent="tex_reg")

    with dpg.window(
        label=os.path.basename(path),
        tag=tag,
        width=min(w+50, 1400),
        height=min(h+50, 1000)
    ):
        dpg.add_image(tex_tag)

# =========================
# ZOOM
# =========================
def zoom_callback(sender, app_data):
    global zoom_scale
    zoom_scale = app_data
    for wt, (bw, bh) in widget_dims.items():
        if dpg.does_item_exist(wt):
            dpg.configure_item(wt, width=int(bw * zoom_scale), height=int(bh * zoom_scale))

# =========================
# RENDER
# =========================
def render_folder(prefix, folder):
    if not folder.exists():
        return

    files = sorted([f for f in os.listdir(folder) if f.endswith(".png")])

    for fname in files:
        path = str(folder / fname)
        tex_tag = f"{prefix}_{fname}"
        widget_tag = f"{prefix}_{fname}_widget"

        mtime = os.path.getmtime(path)

        # update texture if file changed
        if file_mtimes.get(path) != mtime:
            load_image_texture(tex_tag, path)
            file_mtimes[path] = mtime

        if widget_tag not in image_widgets:

            # get real texture size
            w, h = texture_meta.get(tex_tag, (1000, 400))

            # scale to fit nicely in GUI
            max_w = 1000
            max_h = 400

            scale = min(max_w / w, max_h / h)
            base_w = int(w * scale)
            base_h = int(h * scale)
            disp_w = int(base_w * zoom_scale)
            disp_h = int(base_h * zoom_scale)
            widget_dims[widget_tag] = (base_w, base_h)

            with dpg.group(parent=f"{prefix}_container"):

                dpg.add_image(tex_tag, width=disp_w, height=disp_h, tag=widget_tag)

                # caption
                with dpg.group(horizontal=True):
                    dpg.add_text(fname)
                    dpg.add_spacer(width=20)
                    dpg.add_text(describe_plot(fname), color=(150,150,150))

                dpg.add_button(
                    label="Open",
                    callback=show_full_image,
                    user_data=path
                )

                dpg.add_spacer(height=10)  # cleaner than separator

            image_widgets.add(widget_tag)

# =========================
# UPDATE LOOP
# =========================
def update():
    global current_galaxy, last_update

    now = time.time()
    if now - last_update < REFRESH:
        return
    last_update = now

    g = load_latest_galaxy()
    if not g:
        return

    gdir = COCKPIT / g

    changed = False

    if g != current_galaxy:
        dpg.set_value("galaxy_text", f"Galaxy: {g}")
        current_galaxy = g
        changed = True
        # Clear stale widgets and tracking so the new galaxy renders fresh
        image_widgets.clear()
        widget_dims.clear()
        for container_tag in ["raw_container", "interp_container"]:
            children = dpg.get_item_children(container_tag, 1)
            if children:
                for child in children:
                    dpg.delete_item(child)

    for folder in ["raw", "interpreted"]:
        fpath = gdir / folder
        if fpath.exists():
            for f in fpath.glob("*.png"):
                mtime = os.path.getmtime(f)
                if file_mtimes.get(str(f)) != mtime:
                    changed = True
                    break

    if not changed:
        return

    render_folder("raw", gdir / "raw")
    render_folder("interp", gdir / "interpreted")

    state = load_json(gdir / "state/well_state.json")

    notes = "\n".join(state.get("notes", []))
    dpg.set_value("notes_text", notes)

    if "center" in state:
        c = state["center"]
        dpg.set_value(
            "summary_text",
            f"rho_s ~ {c['rho_s']:.3e}   r_s ~ {c['r_s']:.3e}   MBH ~ {c['MBH']:.3e}"
        )

    patch_txt = gdir / "patch/config_patch.txt"
    if patch_txt.exists():
        dpg.set_value("patch_text", patch_txt.read_text())

# =========================
# GUI
# =========================
dpg.create_context()

with dpg.theme() as theme:
    with dpg.theme_component(dpg.mvAll):
        dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (18,18,22))
        dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (28,28,35))
        dpg.add_theme_color(dpg.mvThemeCol_Button, (60,60,80))
        dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (90,90,120))
        dpg.add_theme_color(dpg.mvThemeCol_Text, (220,220,230))
dpg.bind_theme(theme)

with dpg.texture_registry(tag="tex_reg"):
    pass

with dpg.window(label="OSPM Cockpit", width=1800, height=1200):

    with dpg.group(horizontal=True):

        # LEFT PANEL (plots)
        with dpg.child_window(width=400, autosize_y=True, height = 8000):

            dpg.add_text("Galaxy: ---", tag="galaxy_text")
            dpg.add_text("", tag="summary_text")

            dpg.add_slider_float(
                label="Zoom",
                tag="zoom_slider",
                default_value=1.0,
                min_value=0.25,
                max_value=3.0,
                width=300,
                callback=zoom_callback,
            )
            dpg.add_spacer(height=6)

            with dpg.tab_bar():

                with dpg.tab(label="Raw"):
                    with dpg.child_window(tag="raw_container", autosize_x=True, height=2000):
                        pass

                with dpg.tab(label="Interpreted"):
                    with dpg.child_window(tag="interp_container", autosize_x=True, height=2000):
                        pass

        # RIGHT PANEL (info)
        with dpg.child_window(width=-1, autosize_y=True):

            dpg.add_text("Notes")
            dpg.add_input_text(tag="notes_text", multiline=True, width=-1, height=300)

            dpg.add_spacer(height=10)

            dpg.add_text("Patch")
            dpg.add_input_text(tag="patch_text", multiline=True, width=-1, height=250)

dpg.create_viewport(title="LOM Cockpit", width=1600, height=1000)
dpg.setup_dearpygui()
dpg.show_viewport()

while dpg.is_dearpygui_running():
    update()
    dpg.render_dearpygui_frame()

dpg.destroy_context()