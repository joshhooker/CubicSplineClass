#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>

#include "CubicSpline.h"

int main() {

    std::vector<double> x, y;
    std::array<double, 110> xArray, yArray;
    double xArr[110], yArr[110];

    int j = 0;
    for(double i = -2.5; i <= 8.5; i += 0.1) {
        x.push_back(i);
        y.push_back(sin(2.*cos(3.*i)));
        if(j < 110) {
            xArray[j] = i;
            yArray[j] = sin(2.*cos(3.*i));
            xArr[j] = i;
            yArr[j] = sin(2.*cos(3.*i));
            j++;
        }
    }

    int iter = 5000;

    std::cout << "Using " << iter << " iterations" << std::endl << std::endl;

    std::cout << "New Spline Class:" << std::endl;
    std::cout << " Vectors:" << std::endl;
    double totalSetup = 0;
    double totalEval = 0;
    double diff;
    int num;
    for(int i = 0; i < iter; i++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        CubicSpline f(x, y);
        auto t2 = std::chrono::high_resolution_clock::now();
        diff = 0.0;
        num = 0;
        for(double i = -1.; i <= 5.0; i += 0.0034) {
            diff += fabs(sin(2.*cos(3.*i)) - f(i));
            num++;
        }
        auto t3 = std::chrono::high_resolution_clock::now();
        totalSetup += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        totalEval += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    }
    std::cout << "   Avg. Setup Time: " << totalSetup/iter << " us. Avg. Eval Time for "<< num << " calls: " << totalEval/iter << " us. Avg. diff: " << diff/num << std::endl;

    std::cout << " std::Array:" << std::endl;
    totalSetup = 0;
    totalEval = 0;
    for(int i = 0; i < iter; i++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        CubicSpline f(xArray, yArray);
        auto t2 = std::chrono::high_resolution_clock::now();
        diff = 0.0;
        num = 0;
        for(double i = -1.; i <= 5.0; i += 0.0034) {
            diff += fabs(sin(2.*cos(3.*i)) - f(i));
            num++;
        }
        auto t3 = std::chrono::high_resolution_clock::now();
        totalSetup += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        totalEval += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    }
    std::cout << "   Avg. Setup Time: " << totalSetup/iter << " us. Avg. Eval Time for "<< num << " calls: " << totalEval/iter << " us. Avg. diff: " << diff/num << std::endl;

    std::cout << " Array:" << std::endl;
    totalSetup = 0;
    totalEval = 0;
    for(int i = 0; i < iter; i++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        CubicSpline f(xArr, yArr);
        auto t2 = std::chrono::high_resolution_clock::now();
        diff = 0.0;
        num = 0;
        for(double i = -1.; i <= 5.0; i += 0.0034) {
            diff += fabs(sin(2.*cos(3.*i)) - f(i));
            num++;
        }
        auto t3 = std::chrono::high_resolution_clock::now();
        totalSetup += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        totalEval += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    }
    std::cout << "   Avg. Setup Time: " << totalSetup/iter << " us. Avg. Eval Time for "<< num << " calls: " << totalEval/iter << " us. Avg. diff: " << diff/num << std::endl;

    return 0;
}
