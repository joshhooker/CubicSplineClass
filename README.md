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
See the [LICENSE](https://github.com/joshhooker/CubicSplineClass/blob/master/LICENSE.md) file for license rights and limitations (MIT).