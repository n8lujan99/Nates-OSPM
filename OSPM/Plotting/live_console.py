# For use with GPUs only - CPU runs will be too slow for live monitoring
# and will eat up resources that ospm needs to run.
import sys
import time
import re
import glob
import os
import pandas as pd
import numpy as np
from pathlib import Path
from PyQt6 import QtWidgets, QtCore
import pyqtgraph as pg
from pyqtgraph.opengl import GLViewWidget

pg.setConfigOptions(antialias=True)
pg.setConfigOption("background", (10, 10, 15))
pg.setConfigOption("foreground", "w")


# --------------------------------------------------
# Resolve Galaxy
# --------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]
WHICH = REPO_ROOT / "which_galaxy"

if not WHICH.exists():
    raise FileNotFoundError("which_galaxy not found")

GALAXY = WHICH.read_text().strip()
PROFILE_DIR = REPO_ROOT / "Data" / "Gal_Profiles" / GALAXY / "default"


# --------------------------------------------------
# CSV Finder
# --------------------------------------------------

def _timestamp_from_filename(path):
    m = re.search(r"(\d{8})_(\d{6})", os.path.basename(path))
    return int(m.group(1) + m.group(2)) if m else None


def find_latest_daemon_deck():
    c = glob.glob(str(PROFILE_DIR / "*daemon_deck*"))
    if not c:
        raise FileNotFoundError("No daemon deck found")

    ts = [(_timestamp_from_filename(p), p) for p in c]
    ts = [t for t in ts if t[0] is not None]
    return max(ts, key=lambda x: x[0])[1] if ts else max(c, key=os.path.getmtime)


# --------------------------------------------------
# Main Dashboard
# --------------------------------------------------

class LiveConsole(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("OSPM Live Convergence Console")
        self.resize(1600, 900)

        layout = QtWidgets.QGridLayout(self)
        self.setLayout(layout)

        # Plots
        self.plot_chi = pg.PlotWidget(title="χ² vs Run")
        self.plot_land = pg.PlotWidget(title="χ² vs MBH")

        layout.addWidget(self.plot_chi, 0, 0, 1, 2)
        layout.addWidget(self.plot_land, 1, 0, 1, 2)

        # Gauges
        self.accept_gauge = pg.PlotWidget()
        self.stag_bar = pg.PlotWidget(title="Runs Since Improvement")

        layout.addWidget(self.accept_gauge, 0, 2)
        layout.addWidget(self.stag_bar, 1, 2)

        # Stats
        self.stats_label = QtWidgets.QLabel()
        self.stats_label.setStyleSheet("color: white; font-size: 14px;")
        layout.addWidget(self.stats_label, 2, 0, 1, 3)

        # Data State
        self.csv_file = find_latest_daemon_deck()
        self.last_row = 0
        self.df = pd.DataFrame()

        self.best_running = None

        # Timer
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_dashboard)
        self.timer.start(2000)

    # --------------------------------------------------

    def update_dashboard(self):

        try:

            df_full = pd.read_csv(self.csv_file)

            if len(df_full) <= self.last_row:
                return

            new_rows = df_full.iloc[self.last_row:]
            self.last_row = len(df_full)

            if self.df.empty:
                self.df = new_rows.copy()
            else:
                self.df = pd.concat([self.df, new_rows])

            self.df["chi2"] = pd.to_numeric(self.df["chi2"], errors="coerce")
            self.df = self.df.replace([np.inf, -np.inf], np.nan)
            self.df = self.df.dropna(subset=["chi2"])

            if len(self.df) < 10:
                return

            chi = self.df["chi2"].values
            runs = np.arange(len(chi))

            if self.best_running is None:
                self.best_running = np.minimum.accumulate(chi)
            else:
                self.best_running = np.minimum.accumulate(chi)

            best_val = self.best_running[-1]
            best_idx = np.where(self.best_running == best_val)[0][0]
            runs_since_improve = len(chi) - best_idx

            if "status" in self.df:
                acc_rate = np.mean(self.df["status"].tail(100) == "pass")
            else:
                acc_rate = np.nan

            # --------------------------------------------------
            # χ² Plot
            # --------------------------------------------------

            self.plot_chi.clear()
            self.plot_chi.plot(runs, chi, pen=(100, 100, 255, 50))
            self.plot_chi.plot(runs, self.best_running, pen=pg.mkPen("c", width=2))

            # --------------------------------------------------
            # Landscape
            # --------------------------------------------------

            self.plot_land.clear()
            if "MBH" in self.df:
                self.plot_land.plot(
                    self.df["MBH"].values,
                    chi,
                    pen=None,
                    symbol="o",
                    symbolSize=4,
                    symbolBrush=(200, 200, 50, 120),
                )

            # --------------------------------------------------
            # Acceptance Gauge (Arc)
            # --------------------------------------------------

            self.accept_gauge.clear()
            self.accept_gauge.setAspectLocked()
            self.accept_gauge.hideAxis("bottom")
            self.accept_gauge.hideAxis("left")

            theta = np.linspace(np.pi, 2 * np.pi, 200)
            r = 1

            x = r * np.cos(theta)
            y = r * np.sin(theta)

            self.accept_gauge.plot(x, y, pen=(80, 80, 80))

            if not np.isnan(acc_rate):
                theta_val = np.linspace(np.pi, np.pi + np.pi * acc_rate, 200)
                x_val = r * np.cos(theta_val)
                y_val = r * np.sin(theta_val)
                self.accept_gauge.plot(x_val, y_val, pen=pg.mkPen("g", width=4))

            self.accept_gauge.setXRange(-1.2, 1.2)
            self.accept_gauge.setYRange(-1.2, 0.2)

            # --------------------------------------------------
            # Stagnation Bar
            # --------------------------------------------------

            self.stag_bar.clear()
            bg = pg.BarGraphItem(x=[0], height=[runs_since_improve], width=0.5)
            self.stag_bar.addItem(bg)

            # --------------------------------------------------
            # Stats Panel
            # --------------------------------------------------

            self.stats_label.setText(
                f"Galaxy: {GALAXY}   "
                f"Runs: {len(self.df)}   "
                f"Best χ²: {best_val:.4f}   "
                f"Runs Since Improve: {runs_since_improve}   "
                f"Acceptance: {acc_rate:.2f}"
            )

        except Exception as e:
            print("Dashboard error:", e)


# --------------------------------------------------

def run():
    args = parse_args()

    app = QtWidgets.QApplication(sys.argv)
    win = LiveConsole(
        profile_dir=args.profile_dir,
        csv_override=args.csv,
        refresh=args.refresh
    )
    win.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    run()
