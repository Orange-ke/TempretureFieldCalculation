#include <iostream>
#include <chrono>

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
// 步长 和 铸坯的尺寸
const int XStep = 5;
const int YStep = 5;
const int ZStep = 5;

const int Length = 2700;
const int Width  = 420;

const int ZLength = 40000;

// 数组长度
const int arrayLength = 1350 / 5; 

// 计算需要的参数
double Density[arrayLength]; // 密度
double Enthalpy[arrayLength]; // 焓
double Lambda[arrayLength]; // 导热系数
double HEff[arrayLength]; // 综合换热系数, 注意：简单处理了！
double Q[arrayLength]; // 热流密度
double C[arrayLength]; // 比热容

//v int 拉速

int ThermalField[ZLength / ZStep][Width / YStep][Length / XStep]; // 温度场

void calculateConcurrently() {
	// 计算温度场 - start
    // end
}

// 获取等效步长
int getEx(int x) {
	if (x == 0 || x == Length/XStep-1) {
		return 2 * XStep;
	}
	return XStep;
}

int getEy(int y) {
	if (y == 0 || y == Width/YStep-1) {
		return 2 * YStep;
	}
	return YStep;
}

// 计算实际传热系数
double GetLambda(int z, int x1, int y1, int x2, int y2) {
	double K = 0.9; // 修正系数K
	int t1 = ThermalField[z][y1][x1];
	int t2 = ThermalField[z][y2][x2];
	int index1 = (t1-200)/5;
    int index2 = (t2-200)/5;
	// 等效空间步长
	int ex1 = getEx(x1);
	int ex2 = getEx(x2);
	int ey1 = getEy(y1);
	int ey2 = getEy(y2);
	if (x1 != x2) {
		return (K * Lambda[index1] * Lambda[index2] * (ex1+ex2) /
			(Lambda[index1]*(ex2) + Lambda[index2]*(ex1)));
	}	if (y1 != y2) {
		return (K * Lambda[index1] * Lambda[index2] * (ey1+ey2) /
			(Lambda[index1]*(ey2) + Lambda[index2]*(ey1)));
	}
    return 0.0; // input error
}

// 计算时间步长
double GetDeltaT(int z, int x, int y) {
	int t = ThermalField[z][y][x];
	int index = (t - 200) / 5;

	double denominator = 1.0;
	if (x == 0 && y == 0) { // case 1
		denominator = 2*GetLambda(z, x, y, x+1, y)/(XStep*(getEx(x)+getEx(x+1))) +
			2*GetLambda(z, x, y, x, y+1)/(YStep*(getEy(y)+getEy(y+1)));
	} else if (x > 0 && x < Length/XStep-1 && y == 0) { // case 2
		denominator = 2*GetLambda(z, x, y, x-1, y)/(XStep*(getEx(x)+getEx(x-1))) +
			2*GetLambda(z, x, y, x+1, y)/(XStep*(getEx(x)+getEx(x+1))) +
			2*GetLambda(z, x, y, x, y+1)/(YStep*(getEy(y)+getEy(y+1)));
	} else if (x == Length/XStep-1 && y == 0) { // case 3
		denominator = 2*GetLambda(z, x, y, x-1, y)/(XStep*(getEx(x)+getEx(x-1))) +
			2*GetLambda(z, x, y, x, y+1)/(YStep*(getEy(y)+getEy(y+1))) +
			HEff[index]/(XStep);
	} else if (x == 0 && y > 0 && y < Width/YStep-1) { // case 4
		denominator = 2*GetLambda(z, x, y, x+1, y)/(XStep*(getEx(x)+getEx(x+1))) +
			2*GetLambda(z, x, y, x, y+1)/(YStep*(getEy(y)+getEy(y+1))) +
			2*GetLambda(z, x, y, x, y-1)/(YStep*(getEy(y)+getEy(y-1)));
	} else if (x > 0 && x < Length/XStep-1 && y > 0 && y < Width/YStep-1) { // case 5
	    denominator = 2*GetLambda(z, x, y, x-1, y)/(XStep*(getEx(x)+getEx(x-1))) +
			2*GetLambda(z, x, y, x+1, y)/(XStep*(getEx(x)+getEx(x+1))) +
			2*GetLambda(z, x, y, x, y+1)/(YStep*(getEy(y)+getEy(y+1))) +
			2*GetLambda(z, x, y, x, y-1)/(YStep*(getEy(y)+getEy(y-1)));
	} else if (x == Length/XStep-1 && y > 0 && y < Width/YStep-1) { // case6
		denominator = 2*GetLambda(z, x, y, x-1, y)/(XStep*(getEx(x)+getEx(x-1))) +
			2*GetLambda(z, x, y, x, y+1)/(YStep*(getEy(y)+getEy(y+1))) +
			2*GetLambda(z, x, y, x, y-1)/(YStep*(getEy(y)+getEy(y-1))) +
			HEff[index]/(XStep);
	} else if (x == 0 && y == Width/YStep-1) { // case7
		denominator = 2*GetLambda(z, x, y, x+1, y)/(XStep*(getEx(x)+getEx(x+1))) +
			2*GetLambda(z, x, y, x, y-1)/(YStep*(getEy(y)+getEy(y-1))) +
			HEff[index]/(YStep);
	} else if (x > 0 && x < Length/XStep-1 && y == Width/YStep-1) { // case 8
		denominator = 2*GetLambda(z, x, y, x-1, y)/(XStep*(getEx(x)+getEx(x-1))) +
			2*GetLambda(z, x, y, x+1, y)/(XStep*(getEx(x)+getEx(x+1))) +
			2*GetLambda(z, x, y, x, y-1)/(YStep*(getEy(y)+getEy(y-1))) +
			HEff[index]/(YStep);
	} else if (x == Length/XStep-1 && y == Width/YStep-1) {
		denominator = 2*GetLambda(z, x, y, x-1, y)/(XStep*(getEx(x)+getEx(x-1))) +
			2*GetLambda(z, x, y, x, y-1)/(YStep*(getEy(y)+getEy(y-1))) +
			HEff[index]/(XStep) +
			HEff[index]/(YStep);
	}

	return (Density[index] * Enthalpy[index]) / (t * denominator);
}

// 计算时间步长
double calculateTimeStep() {
	// 计算时间步长 - start
	double deltaTArr[5]; 
    deltaTArr[0] = GetDeltaT(0, 0, Width/YStep-1);
    deltaTArr[1] = GetDeltaT(0, 1, Width/YStep-1);
    deltaTArr[2] = GetDeltaT(0, Length/XStep-1, Width/YStep-1);
    deltaTArr[3] = GetDeltaT(0, Length/XStep-1, 1);
    deltaTArr[4] = GetDeltaT(0, Length/XStep-1, 0);
	double min = 1000.0;
	for (int i = 0; i < 5; i++) {
		if (min > deltaTArr[i]) {
            min = deltaTArr[i];
        }
    }
	return min;
	// 计算时间步长 - end
}

// 计算一个点的温度变化
void calculateOnePoint(double deltaT) {
	int indexlt = (ThermalField[0][Width/YStep-1][0] - 200) / 5;
	double deltaHlt = GetLambda(0, 0, Width/YStep-1, 1, Width/YStep-1)*
		(ThermalField[0][Width/YStep-1][1]-ThermalField[0][Width/YStep-1][0])/
		(XStep*(getEx(1)+getEx(0))) +
		GetLambda(0, 0, Width/YStep-1, 0, Width/YStep-1)*
			(ThermalField[0][Width/YStep-2][0]-ThermalField[0][Width/YStep-1][0])/
			(YStep*(getEy(Width/YStep-1)+getEx(Width/YStep-1))) +
		Q[indexlt]/(2*YStep);

	deltaHlt = deltaHlt * (2 * deltaT / Density[indexlt]);

	int T = (ThermalField[0][Width/YStep-1][0]) - deltaHlt/C[indexlt];
	ThermalField[0][Width/YStep-1][0] = T;
	ThermalField[0][Width/YStep-1][0] = 1550;
}

void calculateSerially() {
	auto start = chrono::steady_clock::now();
	for (int t = 0; t < 4000; t++) {
        double deltaT = calculateTimeStep();
		for (int i = 0; i < 84/2*540/2; i++) {
			calculateOnePoint(deltaT);
		}
	}
    auto end =  chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    cout << "串行计算时间: " << elapsed_seconds.count() << endl;
}

int main(int argc, char *argv[]) {
	// 方程初始条件为 T = Tw，Tw为钢水刚到弯月面处的温度。

	// 对于1/4模型，如果不考虑沿着拉坯方向的传热，则每个切片是首切片、中间切片和尾切片均相同，
	// 仅需要图中的四个角部短点、四个边界节点和内部节点的不同，给出9种不同位置的差分方程。

	// 初始化
	// 1. 初始化网格划分的各个节点的初始温度
	auto start = chrono::steady_clock::now();
	cout << "切片数 = " << sizeof(ThermalField) / sizeof(ThermalField[0])  << " 窄边方向网格数 " << sizeof(ThermalField[0]) / sizeof(ThermalField[0][0])  << " 宽边方向网格数 = " << sizeof(ThermalField[0][0]) / sizeof(ThermalField[0][0][0]) << endl;
	for (int z = 0; z < ZLength/ZStep; z++) {
		for (int y = 0; y < Width/YStep; y++) {
			for (int x = 0; x < Length/XStep; x++) {
				ThermalField[z][y][x] = 1550;
			}
		}
	}
	// 2. 导热系数，200℃ 到 1600℃，随温度的上升而下降
	double LambdaStart = 50.0;
	double LambdaIter = (50.0 - 20.0) / (1400 / 5);
	for (int i = 0; i < (1600-200)/5; i++) {
		Lambda[i] = LambdaStart - i*LambdaIter;
	}

	// 3. 密度
	double DensityStart = 8.0;
	double DensityIter = (8.0 - 5.0) / (1400 / 5);
	for (int i = 0; i < (1600-200)/5; i++) {
		Density[i] = DensityStart - i*DensityIter;
	}

	// 4. 焓值
	double EnthalpyStart = 1000.0;
	double EnthalpyStep = (10000.0 - 1000.0) / (1400 / 5);
	for (int i = 0; i < (1600-200)/5; i++) {
		Enthalpy[i] = EnthalpyStart + i*EnthalpyStep;
	}

	// 5. 综合换热系数
	double HEffStart = 5.0;
	double HEffStep = (20.0 - 5.0) / (1400 / 5);
	for (int i = 0; i < (1600-200)/5; i++) {
		HEff[i] = HEffStart + i*HEffStep;
	}

	// 6. 热流密度
	double QStart = 12.0;
	double QStep = (25.0 - 12.0) / (1400 / 5);
	for (int i = 0; i < (1600-200)/5; i++) {
		Q[i] = QStart + i*QStep;
	}

	// 7. 比热容
	double CStart = 0.46;
	double CStep = (0.5 - 0.46) / (1400 / 5);
	for (int i = 0; i < (1600-200)/5; i++) {
		C[i] = CStart + i*CStep;
	}

    auto end =  chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    cout << "初始化时间: " << elapsed_seconds.count() << endl;
    
	// 四个核心一起计算
	calculateConcurrently();
	// 一个核心计算
	calculateSerially();
}
