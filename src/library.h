#ifndef VNA_MATH_LIBRARY_H
#define VNA_MATH_LIBRARY_H

#include <math.h>
#include "Fast4ier.h"
#include "string"

enum LoadType {
    LOAD_R,
    LOAD_L,
    LOAD_C
};

std::string loadTypeStr[] = {"R", "L", "C"};

//enum {
//    LOAD_R = 0,
//    LOAD_L,
//    LOAD_C
//} LOAD_TYPE;

template<typename T>
T *blackman(const int M) {
    T *result = new T[M];
    memset(result, 0, sizeof(T) * M);
    if (M < 1) {
        return result;
    }
    if (M == 1) {
        result[0] = 1.0;
        return result;
    }
    for (int i = 0; i < M; i++) {
        T n = 1 - M + 2 * i;
        result[i] = 0.42 + 0.5 * cos(M_PI * n / (M - 1)) + 0.08 * cos(2.0 * M_PI * n / (M - 1));
    }
    return result;
}


template<typename T>
T *calc_td(complex *windowed_s11, const int n, const int fft_points) {
//    auto *ifft_result = new complex[fft_points];
    auto *input = new complex[fft_points];
    for (int i = 0; i < fft_points; ++i) {
        input[i] = complex(0., 0.);
    }
    for (int i = 0; i < n; ++i) {
        input[i] = windowed_s11[i];
    }

    Fast4::IFFT(input, fft_points);
    auto *result = new T[fft_points];
    for (int i = 0; i < fft_points; ++i) {
//        result[i] = ifft_result[i].abs();
        result[i] = input[i].abs();
    }
    return result;
}


//template<typename T>
//complex *calc_step_response(T *td, const int fft_points) {
////    auto *ifft_result = new complex[fft_points];
//    auto *step = new T[fft_points];
//    memset(step, 1.0, sizeof(T) * fft_points);
//
////    convolve(td, step)
////    sp1 = fft(in1, fshape, axes=axes)
////    sp2 = fft(in2, fshape, axes=axes)
////
////    ret = ifft(sp1 * sp2, fshape, axes=axes)
//    const auto fshape = 2 * fft_points - 1;
//    auto *result = new complex[fshape];
//
//    auto *sp_td = new complex[fshape];
//    for (int i = 0; i < fft_points; ++i) {
//        sp_td[i] = complex(td[i], 0.);
//    }
//    auto *sp_step = new complex[fshape];
//    for (int i = 0; i < fft_points; ++i) {
//        sp_step[i] = complex(step[i], 0.);
//    }
//    Fast4::FFT(sp_td, fshape);
//    Fast4::FFT(sp_step, fshape);
//    for (int i = 0; i < fshape; ++i) {
//        result[i] = sp_td[i] * sp_step[i];
//    }
//    Fast4::IFFT(result, fshape);
//
//
//    return result;
//}

template<typename T,typename FreqT>
T calc_len(T *td, const FreqT step_size, const T v, const int fft_points) {
#define LIGHT_SPEED 299792458
//    auto *ifft_result = new complex[fft_points];
    auto *time_axis = new T[fft_points];
    memset(time_axis, 0.0, sizeof(T) * fft_points);
    for (int i = 0; i < fft_points; ++i) {
        time_axis[i] = LIGHT_SPEED / step_size / fft_points * i * v;
    }
    int td_max_idx = 0;
    T td_max = 0;
    for (int i = 0; i < fft_points; ++i) {
        if (td[i] > td_max) {
            td_max = td[i];
            td_max_idx = i;
        }
    }

    T cable_len = time_axis[td_max_idx] / 2;
//    time_axis = np.linspace(0, 1 / step_size, FFT_POINTS)
//    distance_axis = time_axis * v * speed_of_light
//    index_peak = np.argmax(td)
//    cable_len = round(distance_axis[index_peak] / 2, 3)


    return cable_len;
}


template<typename T>
complex gamma_to_impedance(const complex gamma, T ref_impedance = 50.0) {
    return ((-gamma - 1) / (gamma - 1)) * ref_impedance;
}


//template<typename T>
complex serial_to_parallel(const complex z) {
    auto z_sq_sum = z.re() * z.re() + z.im() * z.im();
    if (z.re() == 0 and z.im() == 0) {

        return {INFINITY, INFINITY};
    }
    if (z.im() == 0) {
        return {z_sq_sum / z.re(), copysign(INFINITY, z_sq_sum)};
    }
    if (z.re() == 0) {
        return {copysign(INFINITY, z_sq_sum), z_sq_sum / z.im()};
    }
    return {z_sq_sum / z.re(), z_sq_sum / z.im()};
}

//def impedance_to_capacitance(z: complex, freq: float) -> float:
//"""Calculate capacitive equivalent for reactance"""
//if freq == 0:
//return -math.inf
//return math.inf if z.imag == 0 else -(1 / (freq * 2 * math.pi * z.imag))
template<typename T>
T impedance_to_capacitance(const complex z, const T freq) {
//    Calculate capacitive equivalent for reactance
    if (freq == 0) {
        return -INFINITY;
    }
    if (z.im() == 0) {
        return INFINITY;
    } else {
        return -(1 / (freq * 2 * M_PI * z.im()));
    }
}


//def impedance_to_inductance(z: complex, freq: float) -> float:
//"""Calculate inductive equivalent for reactance"""
//return 0 if freq == 0 else z.imag * 1 / (freq * 2 * math.pi)
template<typename T>
T impedance_to_inductance(const complex z, const T freq) {
//    Calculate inductive equivalent for reactance
    if (freq == 0) {
        return 0;
    } else {
        return z.im() * 1 / (freq * 2 * M_PI);
    }
}

//x_p_str = cap_p_str if imp_p.imag < 0 else ind_p_str
template<typename T>
T calc_parlc(const complex imp_p, const T freq) {
    if (imp_p.im() < 0) {
        return impedance_to_capacitance(imp_p, freq);
    } else {
        return impedance_to_inductance(imp_p, freq);
    }
}

template<typename FreqT>
LoadType l_or_c(const complex imp_p, const FreqT freq) {
    if (imp_p.im() < 0) {
        return LOAD_C;
    } else {
        return LOAD_L;
    }
}


template<typename T>
T calc_parr(const complex imp_p) {
    return imp_p.re();
}

template<typename T, typename FreqT>
class DataPoint {
public:

    //   Constructors
    DataPoint() : m_z(0., 0.), m_freq(0) {}

    DataPoint(T re, T im, FreqT freq) : m_z(re, im), m_freq(freq) {}

    DataPoint(complex z, FreqT freq) : m_z(z.re(), z.im()), m_freq(freq) {}

    //   Methods
    complex z() const {
        return m_z;
    }

    FreqT freq() const {
        return m_freq;
    }

    complex impedance(T ref_impedance = 50) const {
        return gamma_to_impedance<T>(m_z);
    }

    complex impedance_parallel() const {
        return serial_to_parallel(this->impedance());
    }

    T capacitiveEquivalent(T ref_impedance = 50) const {
        return impedance_to_capacitance<T>(this->impedance(ref_impedance), m_freq);
    }

    T inductiveEquivalent(T ref_impedance = 50) const {
        return impedance_to_inductance<T>(this->impedance(ref_impedance), m_freq);
    }

    LoadType loadType() const {
        return l_or_c<FreqT>(this->impedance_parallel(), m_freq);
    }

    friend std::ostream &operator<<(std::ostream &out, const DataPoint &dp) {
        out << dp.z() << " @ " << dp.freq();
        return out;
    }

//    operator complex<T>*() {
//        // 创建新的complex数组
//        std::complex<double>* complexArray = new std::complex<double>[1];
//        complexArray[0] = value;
//        return complexArray;
//    }


private:
    complex m_z;
    FreqT m_freq;
};

template<typename T, typename FreqT>
T calc_len(DataPoint<T, FreqT> *dps, const int n, const T v, const int fft_points) {

//#define NUM 101
//#define FFT_POINTS 16384
    const FreqT step_size = dps[1].freq() - dps[0].freq();

    complex s11[n];
    for (int i = 0; i < n; ++i) {
        s11[i] = dps[i].z();
    }
    auto window = blackman<T>(n);
    complex windowed_s11[n];
    for (int i = 0; i < n; ++i) {
        windowed_s11[i] = s11[i] * window[i];
    }

    auto td = calc_td<T>(windowed_s11, n, fft_points);
    return calc_len<T>(td, step_size, v, fft_points);

}

#endif //VNA_MATH_LIBRARY_H
