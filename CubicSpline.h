#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <vector>

#define ASSERT_WITH_MESSAGE(condition, message)                                \
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
    CubicSpline(const std::vector<T> &x, const std::vector<T> &y, bool quadratic = false, bool linear = false);

    template<typename T, int N, int M>
    CubicSpline(const T (&x)[N], const T (&y)[M], bool quadratic = false, bool linear = false);

    template<typename T, std::size_t N, std::size_t M>
    CubicSpline(const std::array<T, N> &x, const std::array<T, M> &y, bool quadratic = false, bool linear = false);

    template<typename T>
    void SetPoints(const std::vector<T> &x, const std::vector<T> &y, bool quadratic = false, bool linear = false);

    template<typename T, int N, int M>
    void SetPoints(const T (&x)[N], const T (&y)[M], bool quadratic = false, bool linear = false);

    template<typename T, std::size_t N, std::size_t M>
    void SetPoints(const std::array<T, N> &x, const std::array<T, M> &y, bool quadratic = false, bool linear = false);

    template<typename T>
    double operator()(T x) const;

    ~CubicSpline();

private:
    size_t size_{0};
    std::vector<double> x_vec_, y_vec_;
    std::vector<double> b_vec_, c_vec_, d_vec_;

    void SetSpline(bool quadratic = false, bool linear = false);
    void SetSplineCubic();
    void SetSplineQuadratic();
    void SetSplineLinear();
};


inline CubicSpline::CubicSpline() = default;

template<typename T>
inline CubicSpline::CubicSpline(const std::vector<T> &x, const std::vector<T> &y, const bool quadratic, const bool linear) {
    ASSERT_WITH_MESSAGE(x.size() == y.size(), "In CubicSpline initialization, x vector size != y vector size\n");
    ASSERT_WITH_MESSAGE(x.size() > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = x.size();
    x_vec_ = x;
    y_vec_ = y;
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    SetSpline(quadratic, linear);
}

template<typename T, int N, int M>
inline CubicSpline::CubicSpline(const T (&x)[N], const T (&y)[M], const bool quadratic, const bool linear) {
    ASSERT_WITH_MESSAGE(N == M, "In CubicSpline initialization, x array size != y array size\n");
    ASSERT_WITH_MESSAGE(N > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = N;
    x_vec_.assign(x, x + N);
    y_vec_.assign(y, y + M);
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    SetSpline(quadratic, linear);
}

template<typename T, std::size_t N, std::size_t M>
inline CubicSpline::CubicSpline(const std::array<T, N> &x, const std::array<T, M> &y, const bool quadratic, const bool linear) {
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

    SetSpline(quadratic, linear);
}

template<typename T>
inline void CubicSpline::SetPoints(const std::vector<T> &x, const std::vector<T> &y, const bool quadratic, const bool linear) {
    ASSERT_WITH_MESSAGE(x.size() == y.size(), "In CubicSpline SetPoints, x vector size != y vector size\n");
    ASSERT_WITH_MESSAGE(x.size() > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = x.size();
    x_vec_ = x;
    y_vec_ = y;
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    SetSpline(quadratic, linear);
}

template<typename T, int N, int M>
inline void CubicSpline::SetPoints(const T (&x)[N], const T (&y)[M], const bool quadratic, const bool linear) {
    ASSERT_WITH_MESSAGE(N == M, "In CubicSpline SetPoints, x array size != y array size\n");
    ASSERT_WITH_MESSAGE(N > 1, "In CubicSpline initialization, array size must be larger than 1\n");
    size_ = N;
    x_vec_.assign(x, x + N);
    y_vec_.assign(y, y + M);
    b_vec_.resize(size_);
    c_vec_.resize(size_);
    d_vec_.resize(size_);

    SetSpline(quadratic, linear);
}

template<typename T, std::size_t N, std::size_t M>
inline void CubicSpline::SetPoints(const std::array<T, N> &x, const std::array<T, M> &y, const bool quadratic, const bool linear) {
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

    SetSpline(quadratic, linear);
}

void inline CubicSpline::SetSpline(bool quadratic, bool linear) {
    if (quadratic) {
        SetSplineQuadratic();
    } else if (linear) {
        SetSplineLinear();
    } else {
        SetSplineCubic();
    }
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
        ASSERT_WITH_MESSAGE(x_vec_[i + 1] > x_vec_[i], "In CubicSpline SetSpline, x array is not sorted from smallest to largest\n");
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

    bool changed = false;
    for (unsigned int i = 0; i < size_ - 1; i++) {
        double slope = (y_vec_[i + 1] - y_vec_[i]) / h[i];
        if (slope == 0.0) {
            b_vec_[i] = 0.0;
            b_vec_[i + 1] = 0.0;
            changed = true;
        } else if ((b_vec_[i] >= 0.0 && b_vec_[i + 1] >= 0.0 && slope > 0.0) || (b_vec_[i] <= 0.0 && b_vec_[i + 1] <= 0.0 && slope < 0.0)) {
            double r = sqrt(b_vec_[i] * b_vec_[i] + b_vec_[i + 1] * b_vec_[i + 1]) / fabs(slope);
            if (r > 3.0) {
                b_vec_[i] *= 3.0 / r;
                b_vec_[i + 1] *= 3.0 / r;
                changed = true;
            }
        }
    }
    if (changed) {
        for (int i = size_ - 2; i > -1; i--) {
            c_vec_[i] = z[i] - u[i] * c_vec_[i + 1];
            d_vec_[i] = (c_vec_[i + 1] - c_vec_[i]) / (3. * h[i]);
        }
    }

    d_vec_[0] = 0.;
    d_vec_[size_ - 1] = 0.;
}

inline void CubicSpline::SetSplineQuadratic() {
    std::vector<double> h(size_), alpha(size_), l(size_), z(size_), u(size_);

    h[0] = x_vec_[1] - x_vec_[0];
    l[0] = 1.;
    u[0] = 0.;
    z[0] = 0.;
    l[size_ - 1] = 1.;
    u[size_ - 1] = 0.;
    c_vec_[size_ - 1] = 0.;

    for (unsigned int i = 1; i < size_ - 1; i++) {
        ASSERT_WITH_MESSAGE(x_vec_[i + 1] > x_vec_[i], "In CubicSpline SetSpline, x array is not sorted from smallest to largest\n");
        h[i] = x_vec_[i + 1] - x_vec_[i];

        alpha[i] = 3. / h[i] * (y_vec_[i + 1] - y_vec_[i]) - 3. / h[i - 1] * (y_vec_[i] - y_vec_[i - 1]);
        l[i] = 2. * (x_vec_[i + 1] - x_vec_[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    b_vec_[0] = (y_vec_[1] - y_vec_[0]) / h[0];
    c_vec_[0] = 0.;
    d_vec_[0] = 0.;
    for (int i = 1; i < size_; i++) {
        c_vec_[i] = (y_vec_[i + 1] - y_vec_[i] - b_vec_[i - 1] - 2. * c_vec_[i - 1] * h[i - 1]) / h[i] * h[i];
        b_vec_[i] = (y_vec_[i + 1] - y_vec_[i]) / h[i] - c_vec_[i] * h[i];
        d_vec_[i] = 0.;
    }

    bool changed = false;
    for (unsigned int i = 0; i < size_ - 1; i++) {
        double slope = (y_vec_[i + 1] - y_vec_[i]) / h[i];
        if (slope == 0.0) {
            b_vec_[i] = 0.0;
            b_vec_[i + 1] = 0.0;
            changed = true;
        } else if ((b_vec_[i] >= 0.0 && b_vec_[i + 1] >= 0.0 && slope > 0.0) || (b_vec_[i] <= 0.0 && b_vec_[i + 1] <= 0.0 && slope < 0.0)) {
            double r = sqrt(b_vec_[i] * b_vec_[i] + b_vec_[i + 1] * b_vec_[i + 1]) / fabs(slope);
            if (r > 3.0) {
                b_vec_[i] *= 3.0 / r;
                b_vec_[i + 1] *= 3.0 / r;
                changed = true;
            }
        }
    }
    if (changed) {
        for (int i = size_ - 2; i > -1; i--) {
            c_vec_[i] = (y_vec_[i + 1] - y_vec_[i] - b_vec_[i - 1] - 2. * c_vec_[i - 1] * h[i - 1]) / h[i] * h[i];
        }
    }
}

inline void CubicSpline::SetSplineLinear() {
    std::vector<double> h(size_), alpha(size_), l(size_), z(size_), u(size_);

    h[0] = x_vec_[1] - x_vec_[0];
    l[0] = 1.;
    u[0] = 0.;
    z[0] = 0.;
    l[size_ - 1] = 1.;
    u[size_ - 1] = 0.;
    c_vec_[size_ - 1] = 0.;

    for (unsigned int i = 1; i < size_ - 1; i++) {
        ASSERT_WITH_MESSAGE(x_vec_[i + 1] > x_vec_[i], "In CubicSpline SetSpline, x array is not sorted from smallest to largest\n");
        h[i] = x_vec_[i + 1] - x_vec_[i];

        alpha[i] = 3. / h[i] * (y_vec_[i + 1] - y_vec_[i]) - 3. / h[i - 1] * (y_vec_[i] - y_vec_[i - 1]);
        l[i] = 2. * (x_vec_[i + 1] - x_vec_[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    for (int i = 0; i < size_; i++) {
        c_vec_[i] = 0.;
        d_vec_[i] = 0.;
        b_vec_[i] = (y_vec_[i + 1] - y_vec_[i]) / h[i];
    }
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
