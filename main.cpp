#include <iostream>
#include <chrono>
#include <pthread.h>
#include "unistd.h"

using namespace std;

// 步长设定
// 1. x方向（宽边方向）5mm
// 2. y方向（窄边方向）5mm
// 3. z方向（拉坯方向），板坯切片厚度方向5mm / 10mm

// 参数解释
// 1. 从节点[i, j] 到 [x, y] 实际等效导热系数 lambda
// 2. 每个节点的密度，density
// 3. 边界节点热流密度，Q
// 4. 边界节点综合换热系数，heff

// 全局变量
// 步长 和 铸坯的尺寸，单位mm
const int XStep = 5;
const int YStep = 5;
const int ZStep = 10;

const int Length = 2700 / 2;
const int Width = 420 / 2;

const int edgeWidth = 20;
int step = 1; // 当edgeWidth > 0, step = 2;

const int ZLength = 40000;

// 数组长度
const int arrayLength = 1600 / 5;

// 计算需要的参数
double Density[arrayLength]; // 密度
double Enthalpy[arrayLength]; // 焓
double Lambda[arrayLength]; // 导热系数
double HEff[arrayLength]; // 综合换热系数, 注意：简单处理了！
double Q[arrayLength]; // 热流密度
double C[arrayLength]; // 比热容

float ThermalField[ZLength / ZStep][Width / YStep][Length / XStep];
float ThermalField1[ZLength / ZStep][Width / YStep][Length / XStep];

// 标记目前温度场温度存储在那个三维数组中
bool alternating = false; // 没计算一个 ▲t 进行一次异或运算

//v int 拉速

// 获取等效步长
int getEx(int x) {
    if (x == 0 || x == Length / XStep - 1) {
        return 2 * XStep;
    }
    return XStep;
}

int getEy(int y) {
    if (y == 0 || y == Width / YStep - 1) {
        return 2 * YStep;
    }
    return YStep;
}

// 计算实际传热系数
double GetLambda(int index1, int index2, int x1, int y1, int x2, int y2) {
    float K = 0.9; // 修正系数K
    // 等效空间步长
    int ex1 = x1 == 0 || x1 == Length / XStep - 1 ? 2 * XStep : XStep;
    int ex2 = x2 == 0 || x2 == Length / XStep - 1 ? 2 * XStep : XStep;
    int ey1 = y1 == 0 || y1 == Width / YStep - 1 ? 2 * YStep : YStep;
    int ey2 = y2 == 0 || y2 == Width / YStep - 1 ? 2 * YStep : YStep;
    if (x1 != x2) {
        return (K * Lambda[index1] * Lambda[index2] * (ex1 + ex2) /
                (Lambda[index1] * ex2 + Lambda[index2] * ex1));
    }
    if (y1 != y2) {
        return (K * Lambda[index1] * Lambda[index2] * (ey1 + ey2) /
                (Lambda[index1] * ey2 + Lambda[index2] * ey1));
    }
    return 1.0; // input error
}

// 计算时间步长
double GetDeltaT(int z, int x, int y) {
    float t = ThermalField[z][y][x];
    int index = int(t) / 5 - 1;
    int index1, index2, index3, index4;

    double denominator = 1.0;
    if (x == 0 && y == 0) { // case 1
        index1 = int(ThermalField[z][y][x + 1]) / 5 - 1;
        index2 = int(ThermalField[z][y + 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x + 1, y) / (XStep * (getEx(x) + getEx(x + 1))) +
                      2 * GetLambda(index, index2, x, y, x, y + 1) / (YStep * (getEy(y) + getEy(y + 1)));
    } else if (x > 0 && x < Length / XStep - 1 && y == 0) { // case 2
        index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
        index2 = int(ThermalField[z][y][x + 1]) / 5 - 1;
        index3 = int(ThermalField[z][y + 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x - 1, y) / (XStep * (getEx(x) + getEx(x - 1))) +
                      2 * GetLambda(index, index2, x, y, x + 1, y) / (XStep * (getEx(x) + getEx(x + 1))) +
                      2 * GetLambda(index, index3, x, y, x, y + 1) / (YStep * (getEy(y) + getEy(y + 1)));
    } else if (x == Length / XStep - 1 && y == 0) { // case 3
        index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
        index2 = int(ThermalField[z][y + 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x - 1, y) / (XStep * (getEx(x) + getEx(x - 1))) +
                      2 * GetLambda(index, index2, x, y, x, y + 1) / (YStep * (getEy(y) + getEy(y + 1))) +
                      HEff[index] / (XStep);
    } else if (x == 0 && y > 0 && y < Width / YStep - 1) { // case 4
        index1 = int(ThermalField[z][y][x + 1]) / 5 - 1;
        index2 = int(ThermalField[z][y + 1][x]) / 5 - 1;
        index3 = int(ThermalField[z][y - 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x + 1, y) / (XStep * (getEx(x) + getEx(x + 1))) +
                      2 * GetLambda(index, index2, x, y, x, y + 1) / (YStep * (getEy(y) + getEy(y + 1))) +
                      2 * GetLambda(index, index3, x, y, x, y - 1) / (YStep * (getEy(y) + getEy(y - 1)));
    } else if (x > 0 && x < Length / XStep - 1 && y > 0 && y < Width / YStep - 1) { // case 5
        index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
        index2 = int(ThermalField[z][y][x + 1]) / 5 - 1;
        index3 = int(ThermalField[z][y + 1][x]) / 5 - 1;
        index4 = int(ThermalField[z][y - 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x - 1, y) / (XStep * (getEx(x) + getEx(x - 1))) +
                      2 * GetLambda(index, index2, x, y, x + 1, y) / (XStep * (getEx(x) + getEx(x + 1))) +
                      2 * GetLambda(index, index3, x, y, x, y + 1) / (YStep * (getEy(y) + getEy(y + 1))) +
                      2 * GetLambda(index, index4, x, y, x, y - 1) / (YStep * (getEy(y) + getEy(y - 1)));
    } else if (x == Length / XStep - 1 && y > 0 && y < Width / YStep - 1) { // case6
        index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
        index2 = int(ThermalField[z][y + 1][x]) / 5 - 1;
        index3 = int(ThermalField[z][y - 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x - 1, y) / (XStep * (getEx(x) + getEx(x - 1))) +
                      2 * GetLambda(index, index2, x, y, x, y + 1) / (YStep * (getEy(y) + getEy(y + 1))) +
                      2 * GetLambda(index, index3, x, y, x, y - 1) / (YStep * (getEy(y) + getEy(y - 1))) +
                      HEff[index] / (XStep);
    } else if (x == 0 && y == Width / YStep - 1) { // case7
        index1 = int(ThermalField[z][y][x + 1]) / 5 - 1;
        index2 = int(ThermalField[z][y - 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x + 1, y) / (XStep * (getEx(x) + getEx(x + 1))) +
                      2 * GetLambda(index, index2, x, y, x, y - 1) / (YStep * (getEy(y) + getEy(y - 1))) +
                      HEff[index] / (YStep);
    } else if (x > 0 && x < Length / XStep - 1 && y == Width / YStep - 1) { // case 8
        index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
        index2 = int(ThermalField[z][y][x + 1]) / 5 - 1;
        index3 = int(ThermalField[z][y - 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x - 1, y) / (XStep * (getEx(x) + getEx(x - 1))) +
                      2 * GetLambda(index, index2, x, y, x + 1, y) / (XStep * (getEx(x) + getEx(x + 1))) +
                      2 * GetLambda(index, index3, x, y, x, y - 1) / (YStep * (getEy(y) + getEy(y - 1))) +
                      HEff[index] / (YStep);
    } else if (x == Length / XStep - 1 && y == Width / YStep - 1) { // case 9
        index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
        index2 = int(ThermalField[z][y - 1][x]) / 5 - 1;
        denominator = 2 * GetLambda(index, index1, x, y, x - 1, y) / (XStep * (getEx(x) + getEx(x - 1))) +
                      2 * GetLambda(index, index2, x, y, x, y - 1) / (YStep * (getEy(y) + getEy(y - 1))) +
                      HEff[index] / (XStep) +
                      HEff[index] / (YStep);
    }

    return (Density[index] * Enthalpy[index]) / (t * denominator);
}

// 计算时间步长
double calculateTimeStep() {
    // 计算时间步长 - start
    double deltaTArr[5];
    deltaTArr[0] = GetDeltaT(0, 0, Width / YStep - 1);
    deltaTArr[1] = GetDeltaT(0, 1, Width / YStep - 1);
    deltaTArr[2] = GetDeltaT(0, Length / XStep - 1, Width / YStep - 1);
    deltaTArr[3] = GetDeltaT(0, Length / XStep - 1, 1);
    deltaTArr[4] = GetDeltaT(0, Length / XStep - 1, 0);
    double min = 1000.0; // 模拟一个很大的数
    for (double i : deltaTArr) {
        if (min > i) {
            min = i;
        }
    }
    return min;
    // 计算时间步长 - end
}

// 计算一个left top点的温度变化
void calculatePointLT(double deltaT, int z) {
    int index = int(ThermalField[z][Width / YStep - 1][0]) / 5 - 1;
    int index1 = int(ThermalField[z][Width / YStep - 1][1]) / 5 - 1;
    int index2 = int(ThermalField[z][Width / YStep - 2][0]) / 5 - 1;
    double deltaHlt = GetLambda(index, index1, 0, Width / YStep - 1, 1, Width / YStep - 1) *
                      int(ThermalField[z][Width / YStep - 1][1] - ThermalField[z][Width / YStep - 1][0]) /
                      (XStep * (getEx(1) + getEx(0))) +
                      GetLambda(index, index2, 0, Width / YStep - 1, 0, Width / YStep - 2) *
                      int(ThermalField[z][Width / YStep - 2][0] - ThermalField[z][Width / YStep - 1][0]) /
                      (YStep * (getEy(Width / YStep - 1) + getEx(Width / YStep - 1))) +
                      Q[index] / (2 * YStep);

    deltaHlt = deltaHlt * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][Width / YStep - 1][0] = ThermalField[z][Width / YStep - 1][0] - float(deltaHlt / C[index])
                : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][Width / YStep - 1][0] = ThermalField1[z][Width / YStep - 1][0] - float(deltaHlt / C[index]);
}

// 计算上表面点温度变化
void calculatePointTA(double deltaT, int x, int z) {
    int index = int(ThermalField[z][Width / YStep - 1][x]) / 5 - 1;
    int index1 = int(ThermalField[z][Width / YStep - 1][x - 1]) / 5 - 1;
    int index2 = int(ThermalField[z][Width / YStep - 1][x + 1]) / 5 - 1;
    int index3 = int(ThermalField[z][Width / YStep - 2][x]) / 5 - 1;
    double deltaHta = GetLambda(index, index1, x, Width / YStep - 1, x - 1, Width / YStep - 1) *
                      int(ThermalField[z][Width / YStep - 1][x - 1] - ThermalField[z][Width / YStep - 1][x]) /
                      (XStep * (getEx(x - 1) + getEx(x))) +
                      GetLambda(index, index2, x, Width / YStep - 1, x + 1, Width / YStep - 1) *
                      int(ThermalField[z][Width / YStep - 1][x + 1] - ThermalField[z][Width / YStep - 1][x]) /
                      (XStep * (getEx(x) + getEx(x + 1))) +
                      GetLambda(index, index3, x, Width / YStep - 1, x, Width / YStep - 2) *
                      int(ThermalField[z][Width / YStep - 2][x] - ThermalField[z][Width / YStep - 1][x]) /
                      (YStep * (getEy(Width / YStep - 2) + getEy(Width / YStep - 1))) +
                      Q[index] / (2 * YStep);

    deltaHta = deltaHta * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][Width / YStep - 1][x] = ThermalField[z][Width / YStep - 1][x] - float(deltaHta / C[index])
                : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][Width / YStep - 1][x] = ThermalField1[z][Width / YStep - 1][x] - float(deltaHta / C[index]);
}

// 计算right top点的温度变化
void calculatePointRT(double deltaT, int z) {
    int index = int(ThermalField[z][Width / YStep - 1][Length / XStep - 1]) / 5 - 1;
    int index1 = int(ThermalField[z][Width / YStep - 1][Length / XStep - 2]) / 5 - 1;
    int index2 = int(ThermalField[z][Width / YStep - 2][Length / XStep - 1]) / 5 - 1;
    double deltaHrt = GetLambda(index, index1, Length / XStep - 1, Width / YStep - 1, Length / XStep - 2, Width / YStep - 1) *
                      int(ThermalField[z][Width / YStep - 1][Length / XStep - 2] - ThermalField[z][Width / YStep - 1][Length / XStep - 1]) /
                      (XStep * (getEx(Length / XStep - 2) + getEx(Length / XStep - 1))) +
                      GetLambda(index, index2, Length / XStep - 1, Width / YStep - 1, Length / XStep - 1, Width / YStep - 2) *
                      int(ThermalField[z][Width / YStep - 2][Length / XStep - 1] - ThermalField[z][Width / YStep - 1][Length / XStep - 1]) /
                      (YStep * (getEy(Width / YStep - 2) + getEy(Width / YStep - 1))) +
                      Q[index] / (2 * YStep) +
                      Q[index] / (2 * XStep);

    deltaHrt = deltaHrt * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][Width / YStep - 1][Width / YStep - 1] =
                    ThermalField[z][Width / YStep - 1][Width / YStep - 1] - float(deltaHrt / C[index])
                : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][Width / YStep - 1][Width / YStep - 1] =
                    ThermalField1[z][Width / YStep - 1][Width / YStep - 1] - float(deltaHrt / C[index]);
}

// 计算右表面点的温度变化
void calculatePointRA(double deltaT, int y, int z) {
    int index = int(ThermalField[z][y][Length / XStep - 1]) / 5 - 1;
    int index1 = int(ThermalField[z][y][Length / XStep - 2]) / 5 - 1;
    int index2 = int(ThermalField[z][y - 1][Length / XStep - 1]) / 5 - 1;
    int index3 = int(ThermalField[z][y + 1][Length / XStep - 1]) / 5 - 1;
    double deltaHra = GetLambda(index, index1, Length / XStep - 1, y, Length / XStep - 2, y) *
                      int(ThermalField[z][y][Length / XStep - 2] - ThermalField[z][y][Length / XStep - 1]) /
                      (XStep * (getEx(Length / XStep - 2) + getEx(Length / XStep - 1))) +
                      GetLambda(index, index2, Length / XStep - 1, y, Length / XStep - 1, y - 1) *
                      int(ThermalField[z][y - 1][Length / XStep - 1] - ThermalField[z][y][Length / XStep - 1]) /
                      (YStep * (getEy(y - 1) + getEy(y))) +
                      GetLambda(index, index3, Length / XStep - 1, y, Length / XStep - 1, y + 1) *
                      int(ThermalField[z][y + 1][Length / XStep - 1] - ThermalField[z][y][Length / XStep - 1]) /
                      (YStep * (getEy(y + 1) + getEy(y))) +
                      Q[index] / (2 * XStep);

    deltaHra = deltaHra * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][y][Width / YStep - 1] = ThermalField[z][y][Width / YStep - 1] - float(deltaHra / C[index])
                : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][y][Width / YStep - 1] = ThermalField1[z][y][Width / YStep - 1] - float(deltaHra / C[index]);
}

// 计算right bottom点的温度变化
void calculatePointRB(double deltaT, int z) {
    int index = int(ThermalField[z][0][Length / XStep - 1]) / 5 - 1;
    int index1 = int(ThermalField[z][0][Length / XStep - 2]) / 5 - 1;
    int index2 = int(ThermalField[z][1][Length / XStep - 1]) / 5 - 1;
    double deltaHrb = GetLambda(index, index1, Length / XStep - 1, 0, Length / XStep - 2, 0) *
                      int(ThermalField[z][0][Length / XStep - 2] - ThermalField[z][0][Length / XStep - 1]) /
                      (XStep * (getEx(Length / XStep - 2) + getEx(Length / XStep - 1))) +
                      GetLambda(index, index2, Length / XStep - 1, 0, Length / XStep - 1, 1) *
                      int(ThermalField[z][1][Length / XStep - 1] - ThermalField[z][0][Length / XStep - 1]) /
                      (YStep * (getEy(1) + getEy(0))) +
                      Q[index] / (2 * XStep);

    deltaHrb = deltaHrb * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][0][Width / YStep - 1] = ThermalField[z][0][Width / YStep - 1] - float(deltaHrb / C[index])
                : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][0][Width / YStep - 1] = ThermalField1[z][0][Width / YStep - 1] - float(deltaHrb / C[index]);
}

// 计算下表面点的温度变化
void calculatePointBA(double deltaT, int x, int z) {
    int index = int(ThermalField[z][0][x]) / 5 - 1;
    int index1 = int(ThermalField[z][0][x - 1]) / 5 - 1;
    int index2 = int(ThermalField[z][0][x + 1]) / 5 - 1;
    int index3 = int(ThermalField[z][1][x]) / 5 - 1;

    double deltaHba = GetLambda(index, index1, x, 0, x - 1, 0) *
                      int(ThermalField[z][0][x - 1] - ThermalField[z][0][x]) /
                      (XStep * (getEx(x - 1) + getEx(x))) +
                      GetLambda(index, index2, x, 0, x + 1, 0) *
                      int(ThermalField[z][0][x + 1] - ThermalField[z][0][x]) /
                      (XStep * (getEx(x + 1) + getEx(x))) +
                      GetLambda(index, index3, x, 0, x, 1) *
                      int(ThermalField[z][1][x] - ThermalField[z][0][x]) /
                      (YStep * (getEy(1) + getEy(0)));

    deltaHba = deltaHba * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][0][x] = ThermalField[z][0][x] - float(deltaHba / C[index]) : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][0][x] = ThermalField1[z][0][x] - float(deltaHba / C[index]);
}

// 计算left bottom点的温度变化
void calculatePointLB(double deltaT, int z) {
    int index = int(ThermalField[z][0][0]) / 5 - 1;
    int index1 = int(ThermalField[z][0][1]) / 5 - 1;
    int index2 = int(ThermalField[z][1][0]) / 5 - 1;
    double deltaHlb = GetLambda(index, index1, 1, 0, 0, 0) *
                      int(ThermalField[z][0][1] - ThermalField[z][0][0]) /
                      (XStep * (getEx(0) + getEx(1))) +
                      GetLambda(index, index2, 0, 1, 0, 0) *
                      int(ThermalField[z][1][0] - ThermalField[z][0][0]) /
                      (YStep * (getEy(1) + getEy(0)));

    deltaHlb = deltaHlb * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][0][0] = ThermalField[z][0][0] - float(deltaHlb / C[index]) : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][0][0] = ThermalField1[z][0][0] - float(deltaHlb / C[index]);
}

// 计算左表面点温度的变化
void calculatePointLA(double deltaT, int y, int z) {
    int index = int(ThermalField[z][y][0]) / 5 - 1;
    int index1 = int(ThermalField[z][y][1]) / 5 - 1;
    int index2 = int(ThermalField[z][y - 1][0]) / 5 - 1;
    int index3 = int(ThermalField[z][y + 1][0]) / 5 - 1;
    double deltaHla = GetLambda(index, index1, 1, y, 0, y) *
                      int(ThermalField[z][y][1] - ThermalField[z][y][0]) /
                      (XStep * (getEx(0) + getEx(1))) +
                      GetLambda(index, index2, 0, y - 1, 0, y) *
                      int(ThermalField[z][y - 1][0] - ThermalField[z][y][0]) /
                      (YStep * (getEy(y) + getEy(y - 1))) +
                      GetLambda(index, index3, 0, y + 1, 0, y) *
                      int(ThermalField[z][y + 1][0] - ThermalField[z][y][0]) /
                      (YStep * (getEy(y) + getEy(y + 1)));

    deltaHla = deltaHla * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][y][0] = ThermalField[z][y][0] - float(deltaHla / C[index]) : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][y][0] = ThermalField1[z][y][0] - float(deltaHla / C[index]);
}

// 计算内部点的温度变化
void calculatePointIN(double deltaT, int x, int y, int z) {
    int index = int(ThermalField[z][y][x]) / 5 - 1;
    int index1 = int(ThermalField[z][y][x - 1]) / 5 - 1;
    int index2 = int(ThermalField[z][y][x + 1]) / 5 - 1;
    int index3 = int(ThermalField[z][y - 1][x]) / 5 - 1;
    int index4 = int(ThermalField[z][y + 1][x]) / 5 - 1;
    double deltaHin = GetLambda(index, index1, x - 1, y, x, y) *
                      int(ThermalField[z][y][x - 1] - ThermalField[z][y][x]) /
                      (XStep * (getEx(x) + getEx(x - 1))) +
                      GetLambda(index, index2, x + 1, y, x, y) *
                      int(ThermalField[z][y][x + 1] - ThermalField[z][y][x]) /
                      (XStep * (getEx(x) + getEx(x + 1))) +
                      GetLambda(index, index3, x, y - 1, x, y) *
                      int(ThermalField[z][y - 1][x] - ThermalField[z][y][x]) /
                      (YStep * (getEy(y) + getEy(y - 1))) +
                      GetLambda(index, index4, x, y + 1, x, y) *
                      int(ThermalField[z][y + 1][x] - ThermalField[z][y][x]) /
                      (YStep * (getEy(y) + getEy(y + 1)));

    deltaHin = deltaHin * (2 * deltaT / Density[index]);

    alternating ?
            ThermalField1[z][y][x] = ThermalField[z][y][x] - float(deltaHin / C[index]) : // 需要修改焓的变化到温度变化的映射关系
            ThermalField[z][y][x] = ThermalField1[z][y][x] - float(deltaHin / C[index]);
}

void calculateSerially() {
    auto start = chrono::system_clock::now();
    for (int count = 0; count < 4; count++) {
        double deltaT = calculateTimeStep();
        for (int k = 0; k < ZLength / ZStep; k++) {
            // 先计算点，再计算外表面，再计算里面的点
            calculatePointRT(deltaT, k);
            for (int i = Length / XStep / 2; i < Length / XStep - 1; i++) {
                calculatePointTA(deltaT, i, k);
            }
            for (int j = Width / YStep / 2; j < Width / YStep - 1; j++) {
                calculatePointRA(deltaT, j, k);
            }
            for (int j = Width / YStep - 1 - edgeWidth; j < Width / YStep - 1; j++) {
                for (int i = Length / XStep - 1 - edgeWidth; i < Length / XStep - 1; i++) {
                    calculatePointIN(deltaT, i, j, k);
                }
            }
            for (int j = Width / YStep / 2; j < Width / YStep - 1 - edgeWidth; j++) {
                for (int i = Length / XStep - 1 - edgeWidth; i < Length / XStep - 1; i++) {
                    calculatePointIN(deltaT, i, j, k);
                }
            }
            for (int j = Width / YStep / 2; j < Width / YStep - 1 - edgeWidth; j = j + 2) {
                for (int i = Length / XStep / 2; i < Length / XStep - 1 - edgeWidth; i = i + 2) {
                    calculatePointIN(deltaT, i, j, k);
                }
            }
            alternating ^= 1;
        }
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "串行计算时间: " << elapsed_seconds.count() << endl;
}

void *calculateCase1(void *arg) {
    auto start = chrono::system_clock::now();
    double deltaT = calculateTimeStep();
    int count = 0;
    for (int k = 0; k < ZLength / ZStep; k++) {
        // 先计算点，再计算外表面，再计算里面的点
        calculatePointLT(deltaT, k);
        count++;
        for (int i = 1; i < Length / XStep / 2; i++) {
            calculatePointTA(deltaT, i, k);
            count++;
        }
        for (int j = Width / YStep / 2; j < Width / YStep - 1; j++) {
            calculatePointLA(deltaT, j, k);
            count++;
        }
        for (int j = Width / YStep - 1 - edgeWidth; j < Width / YStep - 1; j++) {
            for (int i = 1; i < 1 + edgeWidth; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = Width / YStep / 2; j < Width / YStep - 1 - edgeWidth; j++) {
            for (int i = 1; i < 1 + edgeWidth; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = Width / YStep - 1 - edgeWidth; j < Width / YStep - 1; j++) {
            for (int i = 1 + edgeWidth; i < Length / XStep / 2; i = i + 1) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        } 
        for (int j = Width / YStep / 2; j < Width / YStep - 1 - edgeWidth; j = j + step) {
            for (int i = 1 + edgeWidth; i < Length / XStep / 2; i = i + step) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        alternating ^= 1;
    }
    cout << "case1:" << count << endl;
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "任务1执行时间: " << elapsed_seconds.count() << endl;
    return nullptr;
}

void *calculateCase2(void *arg) {
    auto start = chrono::system_clock::now();
    double deltaT = calculateTimeStep();
    int count = 0;
    for (int k = 0; k < ZLength / ZStep; k++) {
        // 先计算点，再计算外表面，再计算里面的点
        calculatePointRT(deltaT, k);
        count++;
        for (int i = Length / XStep / 2; i < Length / XStep - 1; i++) {
            calculatePointTA(deltaT, i, k);
            count++;
        }
        for (int j = Width / YStep / 2; j < Width / YStep - 1; j++) {
            calculatePointRA(deltaT, j, k);
            count++;
        }
        for (int j = Width / YStep - 1 - edgeWidth; j < Width / YStep - 1; j++) {
            for (int i = Length / XStep - 1 - edgeWidth; i < Length / XStep - 1; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = Width / YStep / 2; j < Width / YStep - 1 - edgeWidth; j++) {
            for (int i = Length / XStep - 1 - edgeWidth; i < Length / XStep - 1; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = Width / YStep - 1 - edgeWidth; j < Width / YStep - 1; j++) {
            for (int i = Length / XStep / 2; i < Length / XStep - 1 - edgeWidth; i = i + 1) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = Width / YStep / 2; j < Width / YStep - 1 - edgeWidth; j = j + step) {
            for (int i = Length / XStep / 2; i < Length / XStep - 1 - edgeWidth; i = i + step) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        alternating ^= 1;
    }
    cout << "case2:" << count << endl;
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "任务2执行时间: " << elapsed_seconds.count() << endl;
    return nullptr;
}

void *calculateCase3(void *arg) {
    auto start = chrono::system_clock::now();
    double deltaT = calculateTimeStep();
    int count = 0;
    for (int k = 0; k < ZLength / ZStep; k++) {
        // 先计算点，再计算外表面，再计算里面的点
        calculatePointRB(deltaT, k);
        count++;
        for (int i = Length / XStep / 2; i < Length / XStep - 1; i++) {
            calculatePointBA(deltaT, i, k);
            count++;
        }
        for (int j = 1; j < Width / YStep / 2; j++) {
            calculatePointRA(deltaT, j, k);
            count++;
        }
        for (int j = 1; j < 1 + edgeWidth; j++) {
            for (int i = Length / XStep - 1 - edgeWidth; i < Length / XStep - 1; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = 1 + edgeWidth; j < Width / YStep / 2; j++) {
            for (int i = Length / XStep - 1 - edgeWidth; i < Length / XStep - 1; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = 1; j < 1 + edgeWidth; j++) {
            for (int i = Length / XStep / 2; i < Length / XStep - 1 - edgeWidth; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = 1 + edgeWidth; j < Width / YStep / 2; j = j + step) {
            for (int i = Length / XStep / 2; i < Length / XStep - 1 - edgeWidth; i = i + step) {
                calculatePointIN(deltaT, i, j, k);
                count++;   
            }
        }
        alternating ^= 1;
    }
    cout << "case3:" << count << endl;
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "任务3执行时间: " << elapsed_seconds.count() << endl;
    return nullptr;
}

void *calculateCase4(void *arg) {
    auto start = chrono::system_clock::now();
    double deltaT = calculateTimeStep();
    int count = 0;
    for (int k = 0; k < ZLength / ZStep; k++) {
        // 先计算点，再计算外表面，再计算里面的点
        calculatePointLB(deltaT, k);
        count++;
        for (int i = 1; i < Length / XStep / 2; i++) {
            calculatePointBA(deltaT, i, k);
            count++;
        }
        for (int j = 1; j < Width / YStep / 2; j++) {
            calculatePointLA(deltaT, j, k);
            count++;
        }
        for (int j = 1; j < 1 + edgeWidth; j++) {
            for (int i = 1; i < 1 + edgeWidth; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = 1 + edgeWidth; j < Width / YStep / 2; j++) {
            for (int i = 1; i < 1 + edgeWidth; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = 1; j < 1 + edgeWidth; j++) {
            for (int i = 1 + edgeWidth; i < Length / XStep / 2; i++) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        for (int j = 1 + edgeWidth; j < Width / YStep / 2; j = j + step) {
            for (int i = 1 + edgeWidth; i < Length / XStep / 2; i = i + step) {
                calculatePointIN(deltaT, i, j, k);
                count++;
            }
        }
        alternating ^= 1;
    }
    cout << "case4:" << count << endl;
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "任务4执行时间: " << elapsed_seconds.count() << endl;
    return nullptr;
}

void calculateConcurrently() {
    int cpucorenum = sysconf(_SC_NPROCESSORS_CONF);  /*获取核数*/
    printf("system has %i processor(s). \n", cpucorenum);

    auto start = chrono::system_clock::now();
    // Create cpu sets for threads
    cpu_set_t cpu_set_1;
    cpu_set_t cpu_set_2;
    cpu_set_t cpu_set_3;
    cpu_set_t cpu_set_4;
    // Clears set, so that it contains no CPUs.
    CPU_ZERO(&cpu_set_1);
    CPU_ZERO(&cpu_set_2);
    CPU_ZERO(&cpu_set_3);
    CPU_ZERO(&cpu_set_4);
    // Add the CPU cores we want to pin the threads to
    // Init thread attr for 2 cpu sets
    pthread_attr_t attr1, attr2, attr3, attr4;
    pthread_attr_init(&attr1);
    pthread_attr_init(&attr2);
    pthread_attr_init(&attr3);
    pthread_attr_init(&attr4);
    // Add the CPU cores we want to pin the threads to
    CPU_SET(0, &cpu_set_1);
    CPU_SET(1, &cpu_set_2);
    CPU_SET(2, &cpu_set_3);
    CPU_SET(3, &cpu_set_4);
    // 计算温度场 - start
    pthread_t t1;
    pthread_t t2;
    pthread_t t3;
    pthread_t t4;

    pthread_create(&t1, &attr1, calculateCase1, NULL);
    pthread_create(&t2, &attr2, calculateCase2, NULL);
    pthread_create(&t3, &attr3, calculateCase3, NULL);
    pthread_create(&t4, &attr4, calculateCase4, NULL);

    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    pthread_join(t3, NULL);
    pthread_join(t4, NULL);
    // end
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "并行计算时间: " << elapsed_seconds.count() << endl;
}

int main() {
    // 方程初始条件为 T = Tw，Tw为钢水刚到弯月面处的温度。

    // 对于1/4模型，如果不考虑沿着拉坯方向的传热，则每个切片是首切片、中间切片和尾切片均相同，
    // 仅需要图中的四个角部短点、四个边界节点和内部节点的不同，给出9种不同位置的差分方程。

    // 初始化
    // 1. 初始化网格划分的各个节点的初始温度
    auto start = chrono::system_clock::now();

    // 初始化温度场
    for (int z = 0; z < ZLength / ZStep; z++) {
        for (int y = 0; y < Width / YStep; y++) {
            for (int x = 0; x < Length / XStep; x++) {
                ThermalField[z][y][x] = 1600.0;
                ThermalField1[z][y][x] = 1600.0;
            }
        }
    }

    // 2. 导热系数，200℃ 到 1600℃，随温度的上升而下降
    double LambdaStart = 50.0;
    double LambdaIter = (50.0 - 45.0) / int(1600 / 5);
    for (int i = 0; i < arrayLength; i++) {
        Lambda[i] = LambdaStart - i * LambdaIter;
    }

    // 3. 密度
    double DensityStart = 8.0;
    double DensityIter = (8.0 - 7.0) / int(1600 / 5);
    for (int i = 0; i < arrayLength; i++) {
        Density[i] = DensityStart - i * DensityIter;
    }

    // 4. 焓值
    double EnthalpyStart = 1000.0;
    double EnthalpyStep = (10000.0 - 1000.0) / int(1600 / 5);
    for (int i = 0; i < arrayLength; i++) {
        Enthalpy[i] = EnthalpyStart + i * EnthalpyStep;
    }

    // 5. 综合换热系数
    double HEffStart = 5.0;
    double HEffStep = (20.0 - 15.0) / int(1600 / 5);
    for (int i = 0; i < arrayLength; i++) {
        HEff[i] = HEffStart + i * HEffStep;
    }

    // 6. 热流密度
    double QStart = 12.0;
    double QStep = (25.0 - 20.0) / int(1600 / 5);
    for (int i = 0; i < arrayLength; i++) {
        Q[i] = QStart + i * QStep;
    }

    // 7. 比热容
    double CStart = 46.0;
    double CStep = 754.0 / int(1600 / 5);
    for (int i = 0; i < arrayLength; i++) {
        C[i] = CStart + i * CStep;
    }

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "初始化时间: " << elapsed_seconds.count() << endl;

    if (edgeWidth > 0) {
        step = 2;
    }

    // 四个核心一起计算
    calculateConcurrently();
    // 一个核心计算
    calculateSerially();

    return 0;
}
