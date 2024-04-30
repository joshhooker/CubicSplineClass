#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <vector>

#define ASSERT_WITH_MESSAGE(condition, message)                                 \
  do {                                                                         \
    if (!(condition)) {                                                        \
      printf((message));                                                       \
    }                                                                          \
    assert((condition));                                                       \
  } while (false)

class CubicSpline {
public:
    CubicSpline();

    template<typename T>
    CubicSpline(const std::vector<T> &x, const std::vector<T> &y, bool monotonic = false);

    template<typename T, int N, int M>
    CubicSpline(const T (&x)[N], const T (&y)[M], bool monotonic = false);

    template<typename T, std::size_t N, std::size_t M>
    CubicSpline(const std::array<T, N> &x, const std::array<T, M> &y, bool monotonic = false);

    template<typename T>
    void SetPoints(const std::vector<T> &x, const std::vector<T> &y, bool monotonic = false);

    template<typename T, int N, int M>
    void SetPoints(const T (&x)[N], const T (&y)[M], bool monotonic = false);

    template<typename T, std::size_t N, std::size_t M>
    void SetPoints(const std::array<T, N> &x, const std::array<T, M> &y, bool monotonic = false);

    void SetMonotonic() { monotonic_ = true; }

    template<typename T>
    double operator()(T x) const;

    ~CubicSpline();

private:
    size_t size_{0};
    std::vector<double> x_vec_, y_vec_;
    std::vector<double> b_vec_, c_vec_, d_vec_;

    bool monotonic_{false};

    void SetSpline();

    void SetSplineCubic();

    void SetSplineQuadratic();

    void SetSplineLinear();
};


inline CubicSpline::CubicSpline() = default;

template<typename T>
inline CubicSpline::CubicSpline(const std::vector<T> &x, const std::vector<T> &y, const bool monotonic) {
    ASSERT_WITH_MESSAGE(x.size() == y.size(), "In CubicSpline initialization, x vector size != y vector size\n");
    ASSERT_WITH_MESSAGE(x.size() > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = x.size();
    x_vec_ = x;
    y_vec_ = y;
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    monotonic_ = monotonic;

    SetSpline();
}

template<typename T, int N, int M>
inline CubicSpline::CubicSpline(const T (&x)[N], const T (&y)[M], const bool monotonic) {
    ASSERT_WITH_MESSAGE(N == M, "In CubicSpline initialization, x array size != y array size\n");
    ASSERT_WITH_MESSAGE(N > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = N;
    x_vec_.assign(x, x + N);
    y_vec_.assign(y, y + M);
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    monotonic_ = monotonic;

    SetSpline();
}

template<typename T, std::size_t N, std::size_t M>
inline CubicSpline::CubicSpline(const std::array<T, N> &x, const std::array<T, M> &y, const bool monotonic) {
    ASSERT_WITH_MESSAGE(N == M, "In CubicSpline initialization, x array size != y array size\n");
    ASSERT_WITH_MESSAGE(N > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = N;
    x_vec_.resize(size_);
    y_vec_.resize(size_);
    std::copy(x.begin(), x.begin() + size_, x_vec_.begin());
    std::copy(y.begin(), y.begin() + size_, y_vec_.begin());
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    monotonic_ = monotonic;

    SetSpline();
}

template<typename T>
inline void CubicSpline::SetPoints(const std::vector<T> &x, const std::vector<T> &y, const bool monotonic) {
    ASSERT_WITH_MESSAGE(x.size() == y.size(), "In CubicSpline SetPoints, x vector size != y vector size\n");
    ASSERT_WITH_MESSAGE(x.size() > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = x.size();
    x_vec_ = x;
    y_vec_ = y;
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    monotonic_ = monotonic;

    SetSpline();
}

template<typename T, int N, int M>
inline void CubicSpline::SetPoints(const T (&x)[N], const T (&y)[M], const bool monotonic) {
    ASSERT_WITH_MESSAGE(N == M, "In CubicSpline SetPoints, x array size != y array size\n");
    ASSERT_WITH_MESSAGE(N > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = N;
    x_vec_.assign(x, x + N);
    y_vec_.assign(y, y + M);
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    monotonic_ = monotonic;

    SetSpline();
}

template<typename T, std::size_t N, std::size_t M>
inline void CubicSpline::SetPoints(const std::array<T, N> &x, const std::array<T, M> &y, const bool monotonic) {
    ASSERT_WITH_MESSAGE(N == M, "In CubicSpline SetPoints, x array size != y array size\n");
    ASSERT_WITH_MESSAGE(N > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = N;
    x_vec_.resize(size_);
    y_vec_.resize(size_);
    std::copy(x.begin(), x.begin() + size_, x_vec_.begin());
    std::copy(y.begin(), y.begin() + size_, y_vec_.begin());
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    monotonic_ = monotonic;

    SetSpline();
}

void inline CubicSpline::SetSpline() {
    SetSplineCubic();
}

inline void CubicSpline::SetSplineCubic() {
    std::vector<double> h(size_), alpha(size_), l(size_), z(size_), u(size_);

    h[0] = x_vec_[1] - x_vec_[0];
    l[0] = 1.;
    u[0] = 0.;
    z[0] = 0.;
    l[size_ - 1] = 1.;
    u[size_ - 1] = 0.;
    c_vec_[size_ - 1] = 0.;

    for (unsigned int i = 1; i < size_ - 1; i++) {
        ASSERT_WITH_MESSAGE(x_vec_[i + 1] > x_vec_[i],
                            "In CubicSpline SetSpline, x array is not sorted from smallest to largest\n");
        h[i] = x_vec_[i + 1] - x_vec_[i];

        alpha[i] = 3. / h[i] * (y_vec_[i + 1] - y_vec_[i]) - 3. / h[i - 1] * (y_vec_[i] - y_vec_[i - 1]);
        l[i] = 2. * (x_vec_[i + 1] - x_vec_[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    for (int i = size_ - 2; i > -1; i--) {
        c_vec_[i] = z[i] - u[i] * c_vec_[i + 1];
        b_vec_[i] = (y_vec_[i + 1] - y_vec_[i]) / h[i] - h[i] * (c_vec_[i + 1] + 2. * c_vec_[i]) / 3.;
        d_vec_[i] = (c_vec_[i + 1] - c_vec_[i]) / (3. * h[i]);
    }

    // Monotonic correction
    if (monotonic_) {
        bool changed = false;
        for (unsigned int i = 0; i < size_; i++) {
            int i_low = std::max(static_cast<int>(i - 1), 0);
            int i_high = std::min(static_cast<int>(i + 1), static_cast<int>(size_) - 1);

            if (((y_vec_[i_low] <= y_vec_[i]) && (y_vec_[i] <= y_vec_[i_high]) && b_vec_[i] < 0.0) ||
                ((y_vec_[i_low] >= y_vec_[i]) && (y_vec_[i] >= y_vec_[i_high]) && b_vec_[i] > 0.0)) {
                b_vec_[i] = 0.0;
                changed = true;
            }

            double slope = (y_vec_[i_high] - y_vec_[i]) / h[i];
            if (slope == 0.0 && (b_vec_[i] != 0.0 || b_vec_[i + 1] != 0.0)) {
                b_vec_[i] = 0.0;
                b_vec_[i + 1] = 0.0;
                changed = true;
            } else {
                double r = sqrt(b_vec_[i] * b_vec_[i] + b_vec_[i + 1] * b_vec_[i + 1]) / fabs(slope);
                if (r > 3.0) {
                    b_vec_[i] *= 3.0 / r;
                    b_vec_[i + 1] *= 3.0 / r;
                    changed = true;
                }
            }
        }

        if (changed) {
            for (unsigned int i = 0; i < size_; i++) {
                c_vec_[i] =
                    (3.0 * (y_vec_[i + 1] - y_vec_[i]) / h[i] - 2.0 * b_vec_[i] - b_vec_[i + 1]) /
                    h[i];
                d_vec_[i] = ((b_vec_[i + 1] - b_vec_[i]) / (3.0 * h[i]) - 2.0 * c_vec_[i] / 3.0) / h[i];
            }
        }
    }

    d_vec_[0] = 0.;
    d_vec_[size_ - 1] = 0.;
}

template<typename T>
inline double CubicSpline::operator()(T x) const {
    const auto xs = static_cast<double>(x);

    int l = 0;
    int h = size_;
    while (l < h) {
        int mid = (l + h) / 2;
        if (xs <= x_vec_[mid]) {
            h = mid;
        } else {
            l = mid + 1;
        }
    }

    size_t idx = l == 0 ? 0 : l - 1;

    double xi = xs - x_vec_[idx];

    return y_vec_[idx] + b_vec_[idx] * xi + c_vec_[idx] * xi * xi + d_vec_[idx] * xi * xi * xi;
}

inline CubicSpline::~CubicSpline() = default;
