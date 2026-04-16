"""Streamlit Dashboard for Wasatch Photonics spectrometer control."""

from __future__ import annotations

import datetime as dt
import time
from pathlib import Path
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy.signal import savgol_filter
from wasatch.WasatchBus import WasatchBus
from wasatch.WasatchDevice import WasatchDevice

# =========================
# Configuration
# =========================
APP_TITLE = "Wasatch Spectrometer Dashboard"
SCRIPT_DIR = Path(__file__).resolve().parent
WAVENUMBER_CSV = SCRIPT_DIR / "wasatch_wavenumbers.csv"
DEFAULT_SAVE_DIR = SCRIPT_DIR / "output"
SPECTRUM_START_IDX = 394
SPECTRUM_END_TRIM = 98
MIN_TRIMMED_SCANS = 26
TRIM_HEAD_SCANS = 20
TRIM_TAIL_SCANS = 5
TEC_SETPOINT_DEGC = 10
LIVE_VIEW_TIMEOUT_SEC = 3600


# =========================
# Data helpers
# =========================
def load_wavenumbers(csv_path: Path) -> np.ndarray:
    """Load the wavenumber axis from a CSV file."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Wavenumber mapping file not found: {csv_path}. "
            "Add wasatch_wavenumbers.csv next to this script."
        )
    return pd.read_csv(csv_path, header=None).iloc[:, 0].to_numpy()


def crop_spectrum(spectrum: Sequence[float]) -> np.ndarray:
    """Crop detector edge pixels to match the calibrated wavenumber range."""
    arr = np.asarray(spectrum)
    return arr[SPECTRUM_START_IDX:-SPECTRUM_END_TRIM]


def remove_fluo_spectra_lowest_point(
    input_spectrum: Sequence[np.ndarray],
    int_width: int,
    need_smoothing: bool,
    spanner: int = 5,
    polydegree: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Estimate and subtract fluorescence background using a lowest-point method."""
    x_vals = np.asarray(input_spectrum[0])
    y_vals = np.asarray(input_spectrum[1]).copy()

    if need_smoothing:
        y_vals = savgol_filter(y_vals, window_length=spanner, polyorder=polydegree)

    bg = y_vals.copy()
    half_width = int(int_width / 2)

    for i in range(half_width + 1, len(y_vals) - half_width - 1):
        lowest = np.inf
        for j in range(1, half_width + 1):
            for k in range(1, half_width + 1):
                temp = ((y_vals[i + k] - y_vals[i - j]) * j / (k + j)) + y_vals[i - j]
                if temp < lowest:
                    lowest = temp
        bg[i] = lowest

    bg_subtracted = y_vals - bg
    return x_vals, bg_subtracted


def trim_scans(spectra: np.ndarray) -> np.ndarray:
    """Discard startup and tail scans when enough scans were collected."""
    if spectra.shape[0] >= MIN_TRIMMED_SCANS:
        return spectra[TRIM_HEAD_SCANS:-TRIM_TAIL_SCANS]
    return spectra


def build_timestamped_filename(base_name: str, acq_time: int, laser_power: int) -> str:
    """Return a unique output filename preserving the original extension."""
    if not base_name.lower().endswith(".csv"):
        base_name += ".csv"
    stem = Path(base_name).stem
    suffix = Path(base_name).suffix
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{stem}_{acq_time}s_{laser_power}lp_{timestamp}{suffix}"


# =========================
# Hardware helpers
# =========================
@st.cache_resource

def init_device() -> Optional[WasatchDevice]:
    """Connect to the first detected Wasatch device and initialize TEC settings."""
    bus = WasatchBus()
    bus.update()
    if not bus.device_ids:
        return None

    device = WasatchDevice(bus.device_ids[0])
    result = device.connect()
    if not result.data:
        return None

    device.change_setting("detector_tec_setpoint_degC", TEC_SETPOINT_DEGC)
    device.change_setting("detector_tec_enable", True)
    return device


def configure_live_settings(device: WasatchDevice, integration_time: int, laser_power: int, laser_enabled: bool) -> None:
    """Apply live acquisition settings to the spectrometer."""
    device.change_setting("integration_time_ms", integration_time)
    device.change_setting("laser_power_perc", laser_power)
    device.change_setting("set_laser_enable", laser_enabled)


def acquire_single_spectrum(device: WasatchDevice) -> Optional[np.ndarray]:
    """Acquire one spectrum and return it as a NumPy array."""
    result = device.acquire_data()
    if result and result.data and result.data.spectrum:
        return np.asarray(result.data.spectrum)
    return None


def acquire_averaged_spectra(device: WasatchDevice, duration_sec: int) -> np.ndarray:
    """Continuously acquire spectra for the requested duration."""
    start_time = time.time()
    spectra = []
    while time.time() - start_time < duration_sec:
        spectrum = acquire_single_spectrum(device)
        if spectrum is not None:
            spectra.append(spectrum)
        time.sleep(0.1)
    return np.asarray(spectra)


# =========================
# Plot helpers
# =========================
def plot_plotly_spectrum(x_vals: np.ndarray, y_vals: np.ndarray, title: str) -> None:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode="lines", name=title))
    fig.update_layout(
        title=title,
        xaxis_title="Wavenumber (cm⁻¹)",
        yaxis_title="Intensity (a.u.)",
        height=400,
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_matplotlib_live(plot_area, x_vals: np.ndarray, y_vals: np.ndarray, title: str, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_vals, y_vals)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Intensity (a.u.)")
    plot_area.pyplot(fig)
    plt.close(fig)


# =========================
# Streamlit app
# =========================
def main() -> None:
    st.set_page_config(page_title=APP_TITLE, layout="wide")
    st.title(APP_TITLE)

    try:
        wavenumbers = load_wavenumbers(WAVENUMBER_CSV)
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    device = init_device()
    if device is None:
        st.error("No spectrometer found or failed to connect.")
        st.stop()

    if "dark_spectrum" not in st.session_state:
        st.session_state.dark_spectrum = None
    if "live_running" not in st.session_state:
        st.session_state.live_running = False
    if "live_stopped" not in st.session_state:
        st.session_state.live_stopped = False

    # Sidebar
    st.sidebar.header("Live Settings")
    integration_time = st.sidebar.slider("Integration Time (ms)", 100, 10000, 100, step=100)
    laser_power = st.sidebar.slider("Laser Power (%)", 0, 100, 5)
    enable_laser = st.sidebar.checkbox("Laser ON", value=False)
    use_dark_measurement_live = st.sidebar.checkbox("Use Dark Measurement")

    st.sidebar.header("Background Subtraction")
    int_width = st.sidebar.slider("Int Width", 3, 25, 12)
    apply_smoothing = st.sidebar.checkbox("Apply Smoothing", value=True)
    spanner = st.sidebar.slider("Smoothing Window", 3, 21, 5, step=2)
    polydegree = st.sidebar.slider("Smoothing Degree", 1, 5, 1)

    configure_live_settings(device, integration_time, laser_power, enable_laser)

    # Acquisition section
    st.subheader("Acquire and Save")
    acq_time = st.number_input("Acquisition Time (sec)", min_value=1, max_value=60, value=5)
    acq_laser_power = st.slider("Acquisition Laser Power (%)", 0, 100, 10)
    custom_folder = st.text_input("Folder to Save In", value=str(DEFAULT_SAVE_DIR))
    custom_filename = st.text_input("Base Filename", value="spectrum.csv")
    acquire_button = st.button("Acquire and Save Spectrum")
    is_dark_acquisition = st.checkbox("Save as Dark Measurement")
    use_dark_measurement_save = st.checkbox("Use Dark Measurement (Save)")

    if acquire_button:
        st.info("Acquiring...")

        output_dir = Path(custom_folder)
        output_dir.mkdir(parents=True, exist_ok=True)

        final_filename = build_timestamped_filename(custom_filename, int(acq_time), int(acq_laser_power))
        output_path = output_dir / final_filename
        base_stem = output_path.stem

        device.change_setting("laser_power_perc", acq_laser_power)
        device.change_setting("set_laser_enable", True)
        spectra = acquire_averaged_spectra(device, int(acq_time))
        device.change_setting("set_laser_enable", False)

        if spectra.size == 0:
            st.error("No valid spectra acquired.")
            st.stop()

        trimmed_spectra = trim_scans(spectra)
        if spectra.shape[0] < MIN_TRIMMED_SCANS:
            st.info("Fewer scans were captured than the trimming threshold, so all scans were used.")

        mean_spectrum = np.sum(trimmed_spectra, axis=0)
        mean_spectrum_cropped = crop_spectrum(mean_spectrum)
        all_scans_cropped = np.asarray([crop_spectrum(scan) for scan in spectra])

        if len(wavenumbers) != len(mean_spectrum_cropped):
            st.error("Wavenumber file length does not match cropped spectrum length.")
            st.stop()

        if is_dark_acquisition:
            st.session_state.dark_spectrum = mean_spectrum_cropped.copy()
            st.success("Dark measurement saved in memory.")

        if use_dark_measurement_save:
            if st.session_state.dark_spectrum is not None:
                mean_spectrum_cropped = mean_spectrum_cropped - st.session_state.dark_spectrum
                st.info("Dark spectrum subtracted from saved measurement.")
            else:
                st.warning("No dark spectrum available. Acquire one first to enable subtraction.")

        df_mean = pd.DataFrame({
            "Wavenumber": wavenumbers,
            "Intensity": mean_spectrum_cropped,
        })
        df_mean.to_csv(output_path, index=False)
        st.success(f"Saved averaged spectrum to: {output_path}")

        scan_cols = [f"Scan_{i + 1:03d}" for i in range(all_scans_cropped.shape[0])]
        df_all_scans = pd.DataFrame(all_scans_cropped.T, columns=scan_cols)
        df_all_scans.insert(0, "Wavenumber", wavenumbers)

        excel_output_path = output_dir / f"{base_stem}_all_scans.xlsx"
        df_all_scans.to_excel(excel_output_path, index=False, sheet_name="all_scans_cropped")
        st.success(f"Saved all individual scans to: {excel_output_path}")

        clean_x, clean_y = remove_fluo_spectra_lowest_point(
            [df_mean["Wavenumber"].to_numpy(), df_mean["Intensity"].to_numpy()],
            int_width=int_width,
            need_smoothing=apply_smoothing,
            spanner=spanner,
            polydegree=polydegree,
        )

        bksub_output_path = output_dir / f"{base_stem}_bksub.csv"
        pd.DataFrame({"Wavenumber": clean_x, "Intensity": clean_y}).to_csv(bksub_output_path, index=False)
        st.success(f"Saved background-subtracted spectrum to: {bksub_output_path}")

        plot_plotly_spectrum(df_mean["Wavenumber"].to_numpy(), df_mean["Intensity"].to_numpy(), "Acquired Averaged Spectrum (Raw)")
        plot_plotly_spectrum(clean_x, clean_y, "Acquired Spectrum (Background-Subtracted)")

    # Live view section
    st.subheader("Real-Time Spectrum")
    refresh_rate = st.slider("Refresh Rate (sec)", 0.1, 5.0, 0.5, step=0.1)

    if st.button("Toggle Live View"):
        st.session_state.live_running = not st.session_state.live_running
        st.session_state.live_stopped = not st.session_state.live_running

    if st.button("Stop Live View"):
        st.session_state.live_running = False
        st.session_state.live_stopped = True

    if st.session_state.live_running:
        st.success("Live view is ON")
    else:
        st.warning("Live view is OFF")

    plot_area_raw = st.empty()
    plot_area_bgsub = st.empty()

    if st.session_state.live_running and not st.session_state.live_stopped:
        stop_time = time.time() + LIVE_VIEW_TIMEOUT_SEC
        while time.time() < stop_time:
            if not st.session_state.live_running or st.session_state.live_stopped:
                break

            spectrum = acquire_single_spectrum(device)
            if spectrum is None:
                plot_area_raw.warning("No spectrum received.")
                time.sleep(refresh_rate)
                continue

            spectrum_raw = crop_spectrum(spectrum)
            if use_dark_measurement_live and st.session_state.dark_spectrum is not None:
                spectrum_raw = spectrum_raw - st.session_state.dark_spectrum

            if len(wavenumbers) == len(spectrum_raw):
                plot_matplotlib_live(
                    plot_area_raw,
                    wavenumbers,
                    spectrum_raw,
                    "Live Raw Spectrum",
                    "Wavenumber (cm⁻¹)",
                )

                clean_x, clean_y = remove_fluo_spectra_lowest_point(
                    [wavenumbers, spectrum_raw],
                    int_width=int_width,
                    need_smoothing=apply_smoothing,
                    spanner=spanner,
                    polydegree=polydegree,
                )
                plot_matplotlib_live(
                    plot_area_bgsub,
                    clean_x,
                    clean_y,
                    "Live Background-Subtracted Spectrum",
                    "Wavenumber (cm⁻¹)",
                )
            else:
                x_vals = np.arange(len(spectrum_raw))
                plot_matplotlib_live(plot_area_raw, x_vals, spectrum_raw, "Live Raw Spectrum", "Pixel Index")

            time.sleep(refresh_rate)


if __name__ == "__main__":
    main()
