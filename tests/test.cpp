#include "gtest/gtest.h"
#include "library.h"
#include "iostream"
#include "test_data.h"

namespace {
    TEST(blackmanTest, 0) {
        int n = 0;
        auto res = blackman<double>(n);
        for (int i = 0; i < n; ++i) {
            std::cout << res[i] << " ";
        }

    }

    TEST(blackmanTest, 1) {
        int n = 1;
        auto res = blackman<double>(n);
        for (int i = 0; i < n; ++i) {
            std::cout << res[i] << " ";
        }
    }

    TEST(blackmanTest, 5) {
        int n = 5;
        auto res = blackman<double>(n);
        for (int i = 0; i < n; ++i) {
            std::cout << res[i] << " ";
        }
    }

    TEST(blackmanTest, 10) {
        int n = 10;
        auto res = blackman<double>(n);
        for (int i = 0; i < n; ++i) {
            std::cout << res[i] << " ";
        }
    }

    TEST(clac_tdTest, 391) {
#define NUM 101
#define FFT_POINTS 16384
        complex s11[NUM];
        for (int i = 0; i < NUM; ++i) {
            s11[i] = complex(tdr_test_data_391[2 * i], tdr_test_data_391[2 * i + 1]);
        }
        auto window = blackman<FLT>(NUM);
        complex windowed_s11[NUM];
        for (int i = 0; i < NUM; ++i) {
            windowed_s11[i] = s11[i] * window[i];
        }

        auto res = calc_td<FLT>(windowed_s11, NUM, FFT_POINTS);
        for (int i = 0; i < FFT_POINTS; ++i) {
            std::cout << res[i] << " ";
        }
    }


    TEST(clac_lenTest, 391) {
#define NUM 101
#define FFT_POINTS 16384
        complex s11[NUM];
        for (int i = 0; i < NUM; ++i) {
            s11[i] = complex(tdr_test_data_366[2 * i], tdr_test_data_366[2 * i + 1]);
        }
        auto window = blackman<FLT>(NUM);
        complex windowed_s11[NUM];
        for (int i = 0; i < NUM; ++i) {
            windowed_s11[i] = s11[i] * window[i];
        }

        auto td = calc_td<FLT>(windowed_s11, NUM, FFT_POINTS);
//        for (int i = 0; i < FFT_POINTS; ++i) {
//            std::cout << td[i] << " ";
//        }

        auto len = calc_len<FLT>(td, 14850990.0990099, 0.695, FFT_POINTS);
//        for (int i = 0; i < FFT_POINTS; ++i) {
//            std::cout << sr[i] << " ";
//        }
        std::cout << len << std::endl;
    }

    TEST(gamma_to_impTest, 0) {
        std::cout << gamma_to_impedance(complex(0.9, -0.1), 50) << std::endl;
    }

    TEST(serial_to_parallelTest, 0) {
        std::cout << serial_to_parallel(complex(0.9, -0.1)) << std::endl;
    }

    TEST(impedance_to_capacitanceTest, 0) {
        std::cout << impedance_to_capacitance<FLT>(complex(0.9, -0.1), 123500000.0) << std::endl;
    }

    TEST(impedance_to_inductanceTest, 0) {
        std::cout << impedance_to_inductance<FLT>(complex(0.9, -0.1), 123500000.0) << std::endl;
    }

    TEST(calc_parlcTest, 0) {
        std::cout << calc_parlc<FLT>(complex(0.9, -0.1), 123500000.0) << std::endl;
    }

    TEST(data_pointTest, 0) {
        DataPoint<int> dp(0.9, -0.1, 123500000);
        std::cout << dp << std::endl;
        std::cout << "cap: " << dp.capacitiveEquivalent() << std::endl;
        std::cout << "ind: " << dp.inductiveEquivalent() << std::endl;
        std::cout << "typ: " << loadTypeStr[dp.loadType()] << std::endl;
    }


}