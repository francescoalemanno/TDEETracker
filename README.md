# TDEE Kalman Tracker

**Minimal, principled, high-precision daily TDEE estimation using Kalman filtering.**

---

## 🚀 Quick Start Guide

This tool allows you to **track your true Total Daily Energy Expenditure (TDEE)** day-by-day, accounting for noise and uncertainty — in a way that simple spreadsheet estimators or moving averages cannot.

1. **Prepare your log file.**
   
   Create a text file (`log.txt`) with daily entries in the following format:

   ```
   yyyy-mm-dd weight_in_kg calories_in_kcal
   yyyy-mm-dd weight_in_kg calories_in_kcal
   ...
   gw target_goal_weight_in_kg
   ```

   Example:
   ```
   2025-04-10 75.4 2300
   2025-04-11 75.3 2250
   2025-04-12 75.1 2400
   gw 72.0
   ```

2. **Build the project.**

   You need a C++17 compatible compiler, [Eigen](https://eigen.tuxfamily.org/), and [CLI11](https://github.com/CLIUtils/CLI11) are included already in the repo.

   Example using `g++`:
   ```bash
   g++ -O2 -std=c++17 main.cpp -o tdee_tracker
   ```

3. **Run it on your data.**
   
   ```bash
   ./tdee_tracker log.txt
   ```

   Example output:
   ```
   2025-04-10 75.40 2300.00 - TDEE: 2410 ± 30  EstW: 75.40 ± 0.05  ΔW7d: -0.10 ± 0.15
   ...
   
   Goal: 72.0 kg | Current: 75.1 kg
   Suggested intake: 2000 cal/day
   Weekly change: -0.50 kg/week
   ```

4. **Optional: Enable smoothing for even more stable estimates**:

   ```bash
   ./tdee_tracker log.txt --smooth
   ```

---

## ⚙️ Details and Customization

You can fine-tune the model using several advanced flags:

| Option                         | Description |
| :----------------------------- | :---------- |
| `-E, --initialTDEE`             | Set initial value for estimated TDEE (default to first day calories). |
| `-S, --smooth`                  | Apply Rauch-Tung-Striebel smoothing (better estimates by incorporating future data). |
| `-C, --calibrate`               | Recalibrate observation noise from your dataset automatically. |
| `-K, --calPerFatKg`             | Set calories per kg of fat mass (default 7700 kcal/kg). |
| `--mw, --rsdObsWeight`          | Measurement noise level for weight (default 0.008, i.e., ~0.8% error). |
| `--mc, --rsdObsCal`             | Measurement noise level for calories (default 0.1, i.e., ~10% error). |
| `--pw, --rsdWeight`             | Process noise for weight drift (default 0.0001). |
| `--pe, --rsdTDEE`               | Process noise for TDEE drift (default 0.01). |
| `--maxDailyChangePct`           | Maximum safe daily weight change (% of weight per day, default: 0.02/31). |
| `--maxTDEEDeltaPct`             | Maximum safe daily TDEE adjustment (default: 25%). |

If no file is provided, the program reads from `stdin` (standard input).

---

## 📚 Mathematical Foundations

This tracker uses a **Kalman filter** based model — a principled statistical tool for **dynamic state estimation under uncertainty**.

The state we model is a **2-dimensional latent variable**:

- **Body weight** (kg)
- **Total Daily Energy Expenditure** (kcal/day)

The system dynamics are:
- **Weight loss/gain** depends on caloric surplus/deficit.
- **TDEE** is allowed to vary slowly over time (to account for metabolic adaptation, lifestyle changes, etc).

At each step:
- **Prediction**: Use physics and physiology to predict weight changes given calories.
- **Update**: Correct predictions using your observed weight, factoring in measurement noise.

We optionally apply **Rauch-Tung-Striebel smoothing**, which refines estimates backward in time using future information — resulting in a much smoother and more reliable trajectory.

Mathematically:

- State transition:
  
  \[
  x_{k} = F_{k-1} x_{k-1} + u_{k-1} + w_{k-1}
  \]

- Observation:

  \[
  z_{k} = H x_{k} + v_{k}
  \]

where:
- \( w_{k} \sim \mathcal{N}(0, Q) \) is the process noise (uncertainty in weight and TDEE evolution).
- \( v_{k} \sim \mathcal{N}(0, R) \) is the measurement noise (uncertainty in scale readings).

---

### 🚀 Why is this far superior to spreadsheets or naive averages?

- **Noise-aware**: It understands that weight measurements are noisy and uncertain, unlike naive "rolling average" methods.
- **Adaptive**: It models day-to-day TDEE variations realistically.
- **Minimal lag**: Unlike moving averages that trail the real changes, Kalman filters can rapidly pick up trends.
- **Smoothing available**: Incorporates future information for retrospective correction (no moving average can do this properly).
- **Quantified uncertainty**: You get confidence intervals ("±" outputs) for both weight and TDEE estimates.
- **No "manual" recalculations**: Goal-based caloric advice is automatically updated based on the real physical model.

Compared to products like **Macrofactor** or **Excel TDEE spreadsheets** (e.g., nSuns-style estimators), this method:

| Feature                          | TDEE Kalman Tracker | Macrofactor | Excel/nSuns Sheet |
| :------------------------------- | :------------------ | :--------- | :--------------- |
| Noise modeling                   | ✅ | ❌ | ❌ |
| Adaptive TDEE drift               | ✅ | ✅ | ❌ |
| Smoothing using future data       | ✅ | ❌ | ❌ |
| Confidence intervals on estimates | ✅ | ❌ | ❌ |
| Fully offline and local           | ✅ | ❌ | ✅ |
| Transparent mathematics           | ✅ | ❌ | ✅ |
| No black-box "trend weight" magic | ✅ | ❌ | ❌ |

---

## 🧠 Final Notes

- Ideal for **physique athletes**, **weight-loss/gain planners**, and anyone serious about **precise, scientific self-tracking**.
- Designed to be minimal, robust, and fully transparent.
- No external dependencies other than Eigen and CLI11.

---

## 📜 License

This project is released under the GPL3 License.
