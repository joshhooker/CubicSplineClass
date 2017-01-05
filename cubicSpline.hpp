/*
cubicSpline.hpp

cubic spline interpolation. supports std::vector, std::array and standard arrays.

-----------------------------------------------------------------------------
MIT License

Copyright (c) 2017 Joshua Hooker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
-----------------------------------------------------------------------------
*/

#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include <algorithm>
#include <array>
#include <assert.h>
#include <type_traits>
#include <vector>

class CubicSpline {
public:
  CubicSpline(const std::vector<double> &x, const std::vector<double> &y);
  template <typename T, int N, int M> CubicSpline(T (&x) [N], T (&y) [M]);
  template <std::size_t N, std::size_t M> CubicSpline(std::array<double, N>& x, std::array<double, M>& y);
  template <typename T> double operator()(T x) const;
  ~CubicSpline();

private:
  size_t size;
  std::vector<double> xVec, yVec;
  std::vector<double> bVec, cVec, dVec;
};

CubicSpline::CubicSpline(const std::vector<double> &x, const std::vector<double> &y) {
  assert(x.size() == y.size());
  size = x.size();
  xVec = x; yVec = y;
  bVec.resize(size); cVec.resize(size); dVec.resize(size);
  std::vector<double> h(size), alpha(size), l(size), z(size), u(size);

  l[0] = 1.;
  u[0] = 0.;
  z[0] = 0.;
  l[size-1] = 1.;
  u[size-1] = 0.;
  cVec[size-1] = 0.;
  for(int i=0; i<size-1; i++) {
    assert(xVec[i+1]>xVec[i]);
    h[i] = xVec[i+1]-xVec[i];
    if(i>0) {
      alpha[i] = (3./h[i])*(yVec[i+1]-yVec[i]) - (3./h[i-1])*(yVec[i]-yVec[i-1]);
      l[i] = 2.*(xVec[i+1]-xVec[i-1]) - h[i-1]*u[i-1];
      u[i] = h[i]/l[i];
      z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }
  }
  for(int i=size-2; i>-1; i--) {
    cVec[i] = z[i]-u[i]*cVec[i+1];
    bVec[i] = (yVec[i+1]-yVec[i])/h[i] - h[i]*(cVec[i+1]+2.*cVec[i])/3.;
    dVec[i] = (cVec[i+1]-cVec[i])/(3.*h[i]);
  }
}

template <typename T, int N, int M> CubicSpline::CubicSpline(T (&x) [N], T (&y) [M]) {
  assert(N == M);
  size = N;
  xVec.assign(x, x+N); yVec.assign(y, y+M);
  bVec.resize(size); cVec.resize(size); dVec.resize(size);
  std::vector<double> h(size), alpha(size), l(size), z(size), u(size);

  l[0] = 1.;
  u[0] = 0.;
  z[0] = 0.;
  l[size-1] = 1.;
  u[size-1] = 0.;
  cVec[size-1] = 0.;
  for(int i=0; i<size-1; i++) {
    assert(xVec[i+1]>xVec[i]);
    h[i] = xVec[i+1]-xVec[i];
    if(i>0) {
      alpha[i] = (3./h[i])*(yVec[i+1]-yVec[i]) - (3./h[i-1])*(yVec[i]-yVec[i-1]);
      l[i] = 2.*(xVec[i+1]-xVec[i-1]) - h[i-1]*u[i-1];
      u[i] = h[i]/l[i];
      z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }
  }
  for(int i=size-2; i>-1; i--) {
    cVec[i] = z[i]-u[i]*cVec[i+1];
    bVec[i] = (yVec[i+1]-yVec[i])/h[i] - h[i]*(cVec[i+1]+2.*cVec[i])/3.;
    dVec[i] = (cVec[i+1]-cVec[i])/(3.*h[i]);
  }
}

template <std::size_t N, std::size_t M> CubicSpline::CubicSpline(std::array<double, N>& x, std::array<double, M>& y) {
  assert(N == M);
  size = N;
  xVec.resize(size); yVec.resize(size);
  std::copy(x.begin(), x.begin()+size, xVec.begin());
  std::copy(y.begin(), y.begin()+size, yVec.begin());
  bVec.resize(size); cVec.resize(size); dVec.resize(size);
  std::vector<double> h(size), alpha(size), l(size), z(size), u(size);

  l[0] = 1.;
  u[0] = 0.;
  z[0] = 0.;
  l[size-1] = 1.;
  u[size-1] = 0.;
  cVec[size-1] = 0.;
  for(int i=0; i<size-1; i++) {
    assert(xVec[i+1]>xVec[i]);
    h[i] = xVec[i+1]-xVec[i];
    if(i>0) {
      alpha[i] = (3./h[i])*(yVec[i+1]-yVec[i]) - (3./h[i-1])*(yVec[i]-yVec[i-1]);
      l[i] = 2.*(xVec[i+1]-xVec[i-1]) - h[i-1]*u[i-1];
      u[i] = h[i]/l[i];
      z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }
  }
  for(int i=size-2; i>-1; i--) {
    cVec[i] = z[i]-u[i]*cVec[i+1];
    bVec[i] = (yVec[i+1]-yVec[i])/h[i] - h[i]*(cVec[i+1]+2.*cVec[i])/3.;
    dVec[i] = (cVec[i+1]-cVec[i])/(3.*h[i]);
  }
}

template <typename T> double CubicSpline::operator()(T x) const{
  assert(std::is_arithmetic<T>::value);
  double xs = static_cast<double>(x);

  auto id = std::lower_bound(xVec.cbegin(), xVec.cend(), xs);
  auto idx = std::max(int(id-xVec.cbegin())-1, 0);

  double xi = xs-xVec[idx];
  double result;
  if(xs<xVec[0]) {
    result = yVec[0]+bVec[0]*xi+cVec[0]*xi*xi;
  }
  else if(xs>xVec[size-1]) {
    result = yVec[size-1]+bVec[size-1]*xi+cVec[size-1]*xi*xi;
  }
  else {
    result = yVec[idx]+bVec[idx]*xi+cVec[idx]*xi*xi+dVec[idx]*xi*xi*xi;
  }
  return result;
}

CubicSpline::~CubicSpline() {
}

#endif //SPLINE_HPP