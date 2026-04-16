# Wasatch Spectrometer Control Dashboard

A lightweight, real-time control and acquisition dashboard for Wasatch Photonics spectrometers built with Streamlit.

This application was developed to enable full spectrometer control, data acquisition, and processing on systems where the standard ENLIGHTEN software is not practical to run.

---

## Overview

This project provides a browser-based interface for:

- Real-time spectral acquisition
- Laser and integration time control
- Dark measurement capture and subtraction
- Background subtraction using a custom algorithm
- Saving processed and raw spectral data

The dashboard is designed to be simple, responsive, and hardware-focused, enabling reliable operation in constrained computing environments.

---

## Motivation

This dashboard was developed specifically for use on a system running **Windows 10 Enterprise 2016 LTSB**, which has significantly limited storage and RAM capacity.

Due to these constraints, the official Wasatch Photonics software (**ENLIGHTEN**) could not be run reliably. This project replaces that functionality with a lightweight alternative that:

- Minimizes system resource usage  
- Provides direct hardware control  
- Enables reproducible data acquisition workflows  

---

## Features

### Live Control
- Adjustable integration time
- Laser power control
- Laser enable/disable toggle
- Real-time spectral visualization

### Data Acquisition
- Timed acquisition with configurable duration
- Automatic averaging of collected spectra
- Optional dark spectrum subtraction
- Automatic timestamped file naming

### Data Export
- Raw averaged spectra saved as `.csv`
- All individual scans saved as `.xlsx`
- Background-subtracted spectra saved separately

### Signal Processing
- Custom background subtraction algorithm
- Optional Savitzky–Golay smoothing
- Adjustable processing parameters via sidebar controls

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/wasatch-dashboard.git
cd wasatch-dashboard

Here it is properly formatted as clean GitHub-ready Markdown:

````

## 2. Install Dependencies

```bash
pip install -r requirements.txt
````

---

## Required Files

Ensure the following file is present in the project directory:

```
wasatch_wavenumbers.csv
```

This file maps pixel indices to wavenumbers and must match your spectrometer configuration.

---

## Running the Dashboard

```bash
streamlit run wasatch_dashboard.py
```

The application will open in your default web browser.

---

## Usage

### Live Mode

* Toggle live view to stream spectra in real time
* Adjust acquisition parameters in the sidebar

### Acquisition Mode

1. Set acquisition time and laser power
2. Click **Acquire and Save Spectrum**
3. Data will be saved automatically to the specified folder

### Dark Measurement Workflow

1. Acquire a spectrum with the laser off
2. Save as dark measurement
3. Enable dark subtraction for future acquisitions

---

## Output Files

The dashboard generates the following files:

* `*_timestamp.csv` → Averaged spectrum
* `*_all_scans.xlsx` → Individual scans
* `*_bksub.csv` → Background-subtracted spectrum

All files are saved to the specified output directory.

---

## Project Structure

```
.
├── wasatch_dashboard.py
├── wasatch_wavenumbers.csv
├── requirements.txt
├── README.md
└── output/                # Generated data (ignored by git)
```

---

## Notes

* A Wasatch Photonics spectrometer must be connected via USB
* The app assumes a valid wavelength/wavenumber calibration file
* Background subtraction parameters may require tuning depending on sample type


