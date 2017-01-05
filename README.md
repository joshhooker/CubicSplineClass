# Cubic Spline Class
Cubic spline class that is not dependent on other software.

## Usage
CubicSpline accepts std::vector, std::array or standard arrays (i.e. double array[10]).

To call, you must of two of the same type and the same length:
```C++
CubicSpline f(x, y);
```

You are required to have the x array/vector in order.

To get the cubic spline interpolation at position xi, call:
```C++
double result = f(xi)
```
xi can be an integer, float or double but will always return a double.

The spline does output values for positions outside of your initial range, it cannot be well trusted and not recommended.

main.cpp included has an example on how to run the Cubic spline class with a small benchmark.

## License
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