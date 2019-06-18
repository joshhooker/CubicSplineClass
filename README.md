# Cubic Spline Class
Cubic spline class that is not dependent on other software. Uses standard C++ libraries and C++11.

## Usage
CubicSpline accepts std::vector, std::array (C++11) or standard arrays (i.e. double array[10]).

There are two ways to initiate the CubicSpline class:

1. You can begin with an empty constructor such as:
```C++
CubicSpline f;
```
and set the two arrays via:
```C++
f.SetPoints(x, y);
```

2. Can input the two arrays in the constructor:
```C++
CubicSpline f(x, y);
```

You are required to have the x array/vector in order. In a future release, it can sort the vectors if they are not sorted.

To get the cubic spline interpolation at position xi, call:
```C++
double result = f(xi)
```
xi can be an integer, float or double but will always return a double and f is the CubicSpline class.

The spline does output values for values outside of your initial range but it cannot be well trusted and is not recommended.

main.cpp included has an example on how to run the CubicSpline class with a small benchmark.

## License
CubicSpline
Copyright (c) 2017 Josh Hooker

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