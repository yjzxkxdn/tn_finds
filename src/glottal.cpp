// SPDX-License-Identifier: GPL-3.0-only
#include "glottal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

namespace tn::detail {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double rd_min = .02, rd_step = .04, rd_max = 3.02;
constexpr int fit_orders = 80;

// Remapping coordinate at Rd = 0.02 + 0.04*i. Its interpolants reproduce
// the flag's Rd changes.
constexpr std::array<double,76> tension_axis = {
    .02, .17746, .43753, .66692, .84729, .98347, 1.08087, 1.17392,
    1.25649, 1.32908, 1.40173, 1.47488, 1.54063, 1.59777, 1.65342, 1.71218,
    1.77129, 1.82458, 1.86939, 1.90908, 1.94913, 1.99283, 2.03990, 2.08734,
    2.13110, 2.16819, 2.19832, 2.22371, 2.24745, 2.27216, 2.29938, 2.32964,
    2.36271, 2.39785, 2.43399, 2.46993, 2.50445, 2.53652, 2.56546, 2.59096,
    2.61312, 2.63236, 2.64923, 2.66434, 2.67822, 2.69131, 2.70397, 2.71641,
    2.72881, 2.74126, 2.75381, 2.76648, 2.77930, 2.79226, 2.80536, 2.81859,
    2.83193, 2.84541, 2.85898, 2.87267, 2.88645, 2.90032, 2.91425, 2.92820,
    2.94213, 2.95597, 2.96960, 2.98181, 2.98587, 2.99078, 2.99628, 3.00205,
    3.00774, 3.01291, 3.01715, 3.02
};

double shifted_rd(double rd, double mt) {
    const double at = (std::clamp(rd, rd_min, rd_max) - rd_min) / rd_step;
    const int j = std::min(static_cast<int>(tension_axis.size())-2, static_cast<int>(at));
    const double mapped = tension_axis[j] + (at-j)*(tension_axis[j+1]-tension_axis[j]) - .01*mt;
    const auto next = std::upper_bound(tension_axis.begin(), tension_axis.end(), mapped);
    if (next == tension_axis.begin()) return rd_min;
    if (next == tension_axis.end()) return rd_max;
    return rd_min + rd_step*((next-tension_axis.begin()-1) + (mapped-next[-1])/(*next-next[-1]));
}

// Liljencrants-Fant flow derivative, using time in units of one pitch period.
// Rd conversion uses the equations referenced by the author's open ciglet.
class LfPulse {
    double te_, ta_, wg_, sine_, cosine_, epsilon_, alpha_, end_, open_decay_;

    template<class Positive>
    static double bisect(double lo, double hi, Positive positive) {
        for (int j = 0; j < 48; ++j) {
            const double mid = (lo+hi)*.5;
            if (positive(mid)) lo = mid;
            else hi = mid;
        }
        return (lo+hi)*.5;
    }

public:
    explicit LfPulse(double rd) {
        double rk, rg;
        if (rd < 2.7) {
            ta_ = rd < .21 ? 1.e-6 : (-1.+4.8*rd)*.01;
            rk = (.224+.118*rd);
            rg = .25*rk/(.11*rd/(.5+1.2*rk)-ta_);
        } else {
            ta_ = .323/rd;
            const double oq = 1.-1./(2.17*rd);
            rg = .0093552+5.96/(7.96-2.*oq);
            rk = 2.*rg*oq-1.0428;
        }
        te_ = (1.+rk)/(2.*rg);
        wg_ = 2.*pi*rg;
        sine_ = std::sin(wg_*te_);
        cosine_ = std::cos(wg_*te_);
        const double tail = 1.-te_;
        epsilon_ = bisect(0., 1./ta_, [&](double e) { return -std::expm1(-e*tail) > e*ta_; });
        end_ = std::exp(-epsilon_*tail);
        const double return_area = (1.-tail*end_/ta_)/epsilon_;
        alpha_ = bisect(0., 2048., [&](double a) {
            const double open_area = -(a*sine_-wg_*cosine_+wg_*std::exp(-a*te_)) /
                                     (sine_*(a*a+wg_*wg_));
            return open_area > return_area;
        });
        open_decay_ = wg_*std::exp(-alpha_*te_);
    }

    std::complex<double> coefficient(int order) const {
        using C = std::complex<double>;
        const double w = 2.*pi*order;
        const C q(0., w), z = std::polar(1., -w*te_);
        const C u = alpha_-q;
        // Fourier integrals of the damped sinusoidal opening and exponential
        // return. At an integer harmonic, exp(-i*w*(1-Te)) = conj(z).
        const C open = -((u*sine_-wg_*cosine_)*z+open_decay_) /
                        (sine_*(u*u+wg_*wg_));
        const C close = -((z-end_)/(epsilon_+q)-end_*(z-1.)/q)/(epsilon_*ta_);
        return open+close;
    }

    double amplitude(int order) const {
        return std::max(1.e-14, std::abs(coefficient(order)));
    }
};

using Candidates = std::array<double, tension_axis.size()>;

struct FitShapes {
    // Keep candidates contiguous so the fit can vectorize independent sums
    // without changing the addition order within each candidate.
    std::array<Candidates, fit_orders> inverse;
    std::array<std::array<double, fit_orders>, tension_axis.size()> scale, offset;
    FitShapes() {
        for (std::size_t j = 0; j < tension_axis.size(); ++j) {
            const LfPulse pulse(rd_min+rd_step*j);
            double sum = 0., sum_log = 0.;
            for (int k = 0; k < fit_orders; ++k) {
                const double a = pulse.amplitude(k+1)/(k+1);
                const double p = a*a;
                inverse[k][j] = 1./p;
                sum += p;
                sum_log += std::log(p);
                scale[j][k] = sum;
                offset[j][k] = sum_log-(k+1)*std::log(sum);
            }
        }
    }
};

const FitShapes& fit_shapes() {
    static const FitShapes shapes;
    return shapes;
}
} // namespace

double GlottalTrack::sample(double at) const {
    if (rd.empty()) return 0.;
    const double i = std::clamp(at/hop, 0., double(rd.size()-1));
    const auto j = static_cast<std::size_t>(i);
    const auto next = std::min(j+1, rd.size()-1);
    if (rd[j] <= 0. || rd[next] <= 0.) return rd[i-j < .5 ? j : next];
    return rd[j]+(i-j)*(rd[next]-rd[j]);
}

double estimate_rd(const std::vector<Complex>& spectrum, int harmonics, double f0) {
    const int n = std::min({harmonics, fit_orders, static_cast<int>(8000./f0)});
    if (n < 4) return 0.;
    // Model matching uses flow power, after removing the lip radiation curve.
    constexpr double rr = 128./(9.*pi*pi), lr = 8.*.015/(3.*pi*340.);
    const double lip = 2.*pi*f0*lr/rr;
    std::array<double,fit_orders> power;
    double sum = 0., higher = 0.;
    for (int j = 0; j < n; ++j) {
        const double k = j+1.;
        power[j] = std::norm(spectrum[j+1])*(1./(k*k)+lip*lip);
        sum += power[j];
        if (j > 0) higher += power[j];
    }
    if (sum < 1.e-12 || higher < sum*1.e-5) return 0.;
    for (int j = 0; j < n; ++j) power[j] = std::max(1.e-10, power[j]/sum);
    const auto& shapes = fit_shapes();
    Candidates ratios{};
    for (int k = 0; k < n; ++k) {
        const auto& inverse = shapes.inverse[k];
        for (std::size_t j = 0; j < ratios.size(); ++j)
            ratios[j] += power[k]*inverse[j];
    }
    Candidates distance;
    int best = 0;
    for (int j = 0; j < static_cast<int>(distance.size()); ++j) {
        // Itakura-Saito divergence, omitting terms common to all candidates.
        distance[j] = ratios[j]*shapes.scale[j][n-1]+shapes.offset[j][n-1];
        if (distance[j] < distance[best]) best = j;
    }
    double fine = 0.;
    if (best > 0 && best+1 < static_cast<int>(distance.size())) {
        const double l = distance[best-1]-distance[best], r = distance[best+1]-distance[best];
        if (l+r > 1.e-12) fine = .5*(l-r)/(l+r);
    }
    return rd_min+rd_step*(best+fine);
}

void lf_spectrum(double rd, int harmonics, std::vector<std::complex<double>>& out) {
    const LfPulse pulse(rd);
    out.resize(harmonics+1);
    out[0] = {};
    for (int k = 1; k <= harmonics; ++k) out[k] = pulse.coefficient(k);
}

void tension_gain(double rd, double mt, int harmonics, std::vector<double>& log_gain) {
    log_gain.assign(harmonics+1, 0.);
    if (mt == 0. || rd <= 0. || harmonics < 2) return;
    const double edited = shifted_rd(rd, mt);
    if (edited == rd) return;
    const LfPulse original(rd), changed(edited);
    const double reference = std::log(original.amplitude(1)/changed.amplitude(1));
    for (int k = 2; k <= harmonics; ++k)
        log_gain[k] = std::log(changed.amplitude(k)/original.amplitude(k))+reference;
}
} // namespace tn::detail
