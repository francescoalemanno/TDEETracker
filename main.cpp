
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "CLI11.hpp"
#include "Eigen/Dense"

struct LogEntry {
    std::chrono::system_clock::time_point date;
    double weight;
    double cals;
    double time;
};

struct Estimate {
    double weight;
    double tdee;
    double sdw;
    double sdtdee;
};

struct Params {
    double calPerFatKg = 7700.0;
    double rsdObsWeight = 0.008;
    double rsdObsCal = 0.1;
    double rsdWeight = 0.0001;
    double rsdTDEE = 0.01;
    double initialTDEE = -1;
    bool smoothTracking = false;
};

std::tm parseDate(const std::string &dateStr) {
    std::tm tm = {};
    std::istringstream ss(dateStr);
    ss >> std::get_time(&tm, "%Y-%m-%d");
    if (ss.fail()) {
        throw std::runtime_error("Failed to parse date: " + dateStr);
    }
    return tm;
}

std::vector<LogEntry> parseLog(std::istream &input, double &goalWeight) {
    std::vector<LogEntry> entries;
    std::string line;
    std::unordered_set<std::string> seen;
    goalWeight = std::numeric_limits<double>::quiet_NaN();

    while (std::getline(input, line)) {
        line.erase(line.find_last_not_of(" \n\r\t") + 1);
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<std::string> fields{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};

        if (fields.size() == 2 && fields[0] == "gw") {
            goalWeight = std::stod(fields[1]);
            if (goalWeight <= 0)
                throw std::runtime_error("Invalid goal weight");
        } else if (fields.size() >= 3) {
            if (seen.count(fields[0]))
                throw std::runtime_error("Duplicate date: " + fields[0]);
            seen.insert(fields[0]);

            std::tm tm = parseDate(fields[0]);
            auto date =
                std::chrono::system_clock::from_time_t(std::mktime(&tm));

            double weight = std::stod(fields[1]);
            double cals = std::stod(fields[2]);
            if (weight <= 0 || cals < 0)
                throw std::runtime_error("Invalid weight or calories");

            entries.push_back({date, weight, cals, 0});
        } else {
            throw std::runtime_error("Malformed line: " + line);
        }
    }

    std::sort(
        entries.begin(), entries.end(),
        [](const LogEntry &a, const LogEntry &b) { return a.date < b.date; });

    auto startDate = entries[0].date;

    for (int i = 0; i < entries.size(); i++) {
        entries[i].time = std::chrono::duration<double, std::ratio<86400>>(
                              entries[i].date - startDate)
                              .count();
    }
    return entries;
}

std::vector<Estimate> kalman(const std::vector<LogEntry> &E, const Params &P) {
    using namespace Eigen;
    if (E.size() < 3) throw std::runtime_error("Need at least 3 entries");

    Matrix2d I = Matrix2d::Identity();
    RowVector2d H(1.0, 0.0);

    Vector2d X(E[0].weight, E[0].cals);
    if (P.initialTDEE>0 && (!std::isnan(P.initialTDEE))) {
        X[1] = P.initialTDEE;
    }
    Matrix2d V;
    V << pow(P.rsdObsWeight * X[0], 2)+pow(P.rsdWeight * X[0], 2),      0, 
         0,               pow(P.rsdTDEE * X[1], 2)+pow(P.rsdObsCal * X[1], 2);

    struct State {
        Vector2d X;
        Matrix2d V;
    };

    std::vector<State> filt = {{X, V}};
    std::vector<State> pred = filt, smooth = filt;

    for (size_t i = 1; i < E.size(); ++i) {
        double dt = E[i].time - E[i - 1].time;
        double avgCals = (E[i - 1].cals + E[i].cals) / 2.0;

        const double Q00 = pow(P.rsdWeight * std::max(E[i].weight, 30.0), 2)  // Weight process noise
                           + pow(dt/P.calPerFatKg*avgCals*P.rsdObsCal,2);     // Calories control noise
        const double Q11 = pow(P.rsdTDEE * std::max(E[i].cals, 1200.0), 2);
        Matrix2d Q;
        Q << Q00, 0, 0,  Q11;

        Matrix<double, 1, 1> R(pow(P.rsdObsWeight * E[i].weight, 2));


        Matrix2d F;
        F << 1.0, -dt / P.calPerFatKg, 0.0, 1.0;

        Matrix2d Vpred = F * V * F.transpose() + Q;
        auto S = H * Vpred * H.transpose() + R;
        auto K = Vpred * H.transpose() * S.inverse();
        Vector2d u(dt * avgCals / P.calPerFatKg, 0.0);
        Vector2d Xpred = F * X + u;

        Vector<double, 1> z(E[i].weight);
        Vector<double, 1> delta = z - H * Xpred;

        V = (I - K * H) * Vpred;
        V = 0.5 * (V + V.transpose());  // Symmetrize for numerical stability
        X = Xpred + K * delta;

        filt.push_back({X, V});
        pred.push_back({Xpred, Vpred});
        smooth.push_back({X, V});
    }

    if (P.smoothTracking) {
        for (int i = E.size() - 2; i >= 0; --i) {
            double dt = E[i + 1].time - E[i].time;

            Matrix2d F;
            F << 1.0, -dt / P.calPerFatKg, 0.0, 1.0;

            Matrix2d G = filt[i].V * F.transpose() * pred[i + 1].V.inverse();

            smooth[i].X += G * (smooth[i + 1].X - pred[i + 1].X);
            smooth[i].V +=
                G * (smooth[i + 1].V - pred[i + 1].V) * G.transpose();
        }
    }

    std::vector<Estimate> final;
    for (const auto &s : smooth) {
        final.push_back({s.X(0), s.X(1), sqrt(s.V(0, 0)), sqrt(s.V(1, 1))});
    }
    return final;
}

void printAdvice(double goal, double currWeight, double currTDEE,
                 double maxTDEEDeltaPct, double maxDailyChangePct,
                 double calPerFatKg) {
    if (std::isnan(goal)) return;

    double deltaMax = std::min(currTDEE * maxTDEEDeltaPct,
                          currWeight * maxDailyChangePct * calPerFatKg);

    double delta =
        std::clamp((goal - currWeight) * calPerFatKg, -deltaMax, deltaMax);
    double suggested = currTDEE + delta;

    std::cout << "\nGoal: " << std::fixed << std::setprecision(1) << goal
              << " kg | Current: " << currWeight << " kg\n"
              << "Suggested intake: " << std::round(suggested) << " cal/day\n"
              << "Weekly change: " << std::fixed << std::setprecision(2)
              << (suggested - currTDEE) / calPerFatKg * 7 << " kg/week\n";
}

// Calibration routine
Params calibrateParameters(std::vector<LogEntry> &entries, const Params &P) {
    if (entries.size() < 10) {
        std::cerr << "Warning: calibration is not possible with less than 10 "
                     "entries.\n";
        return P;
    }
    Params B = P;
    double wrd = 0.0;
    double n = 0.0;
    for (int i = 1; i < entries.size() - 1; i++) {
        double rd1 = (entries[i - 1].weight - entries[i].weight) /
                     (entries[i - 1].time - entries[i].time);
        double rd2 = (entries[i + 1].weight - entries[i].weight) /
                     (entries[i + 1].time - entries[i].time);
        double rdd = (rd1 - rd2) / (entries[i - 1].time - entries[i + 1].time) *
                     std::sqrt(2.0 / 3.0);
        double avd = (entries[i].weight + entries[i - 1].weight +
                      entries[i + 1].weight) /
                     3.0;
        double ratio = rdd / avd;
        wrd += ratio * ratio;
        n += 1.0;
    }
    B.rsdObsWeight = std::sqrt(wrd / n);
    return B;
}

int main(int argc, char **argv) {
    CLI::App app{"TDEE Kalman Tracker"};
    Params params;

    double maxDailyChangePct = 0.02 / 31.0;
    double maxTDEEDeltaPct = 0.25;
    bool calibrate = false;
    std::string filename;

    app.add_option("file", filename, "Data file to read (optional)");
    app.add_option("-E,--initialTDEE", params.initialTDEE,
                 "Initial value for estimated TDEE, by default it is set equal to the first day calories.")
        ->capture_default_str();
    app.add_flag("-S,--smooth", params.smoothTracking,
                 "Apply Rauch-Tung-Striebel smoothing")
        ->capture_default_str();
    app.add_option("-K,--calPerFatKg", params.calPerFatKg)
        ->capture_default_str();
    app.add_flag("-C,--calibrate", calibrate, "Enable calibration routine")
        ->capture_default_str();
    app.add_option("--mw,--rsdObsWeight", params.rsdObsWeight,
                   "Relative measurement noise for Weight (advanced)")
        ->capture_default_str();
    app.add_option("--mc,--rsdObsCal", params.rsdObsCal,
                   "Relative measurement noise for Calories (advanced)")
        ->capture_default_str();
    app.add_option("--pw,--rsdWeight", params.rsdWeight,
                   "Relative process noise for Weight (advanced)")
        ->capture_default_str();
    app.add_option("--pe,--rsdTDEE", params.rsdTDEE,
                   "Relative process noise for TDEE (advanced)")
        ->capture_default_str();
    app.add_option("--maxDailyChangePct", maxDailyChangePct,
                   "maximum daily recommended relative weight change goal.")
        ->capture_default_str();
    app.add_option("--maxTDEEDeltaPct", maxTDEEDeltaPct,
                   "maximum daily recommended relative calories change goal.")
        ->capture_default_str();
    CLI11_PARSE(app, argc, argv);

    std::istream *input;
    std::ifstream file;
    if (!filename.empty()) {
        file.open(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file\n";
            return 1;
        }
        input = &file;
    } else {
        input = &std::cin;
    }

    double goalWeight;
    auto entries = parseLog(*input, goalWeight);

    if (calibrate) {
        params = calibrateParameters(entries, params);
    }

    auto estimates = kalman(entries, params);

    for (size_t i = 0; i < entries.size(); ++i) {
        size_t days = std::min(i, size_t(7));
        double dt = (entries[i].time - entries[i - days].time);
        double scale = 7.0 / dt;
        double diff =
            (estimates[i].weight - estimates[i - days].weight) * scale;
        double err =
            std::hypot(estimates[i].sdw, estimates[i - days].sdw) * scale;

        auto time = std::chrono::system_clock::to_time_t(entries[i].date);
        std::cout << std::put_time(std::localtime(&time), "%Y-%m-%d") << " "
                  << std::fixed << std::setprecision(2) << entries[i].weight
                  << " " << (entries[i].cals)
                  << " - TDEE: " << (estimates[i].tdee) << " ± "
                  << (estimates[i].sdtdee) << "  EstW: " << estimates[i].weight
                  << " ± " << estimates[i].sdw << "  ΔW7d: " << diff << " ± "
                  << err << "\n";
    }

    printAdvice(goalWeight, estimates.back().weight, estimates.back().tdee,
                maxTDEEDeltaPct, maxDailyChangePct, params.calPerFatKg);

    return 0;
}
