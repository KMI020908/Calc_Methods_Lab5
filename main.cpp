#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>
#include"testFuncs.h"
#include<algorithm>

// Процедура проверки решений уравнения методами бисекции и НьютонаЫ 
template<typename Type>
void checkTestEquations(Type (*f)(Type x), Type firstX, Type lastX, Type x0,
const std::string &B_EQ_FILE_PATH, const std::string &N_EQ_FILE_PATH, Type accuracy = 1e-6, Type h = 1e-4, std::size_t stopIt = 10000){
    // Метод биссекции
    Type solBisection = getEquationSolutionBissection(f, firstX, lastX, accuracy, stopIt);
    std::size_t numOfItersB = getIterationsBissection(f, firstX, lastX, accuracy, stopIt);
    writeEqResBissection(solBisection, firstX, lastX, accuracy, B_EQ_FILE_PATH);
    writeIters(numOfItersB, B_EQ_FILE_PATH);

    // Метод Ньютона
    Type solNewthon = getEquationSolutionNewthon(f, x0, accuracy, h, stopIt);
    std::size_t numOfItersN = getIterationsNewthon(f, x0, accuracy, h, stopIt);
    writeEqResNewthon(solNewthon, x0, accuracy, N_EQ_FILE_PATH);
    writeIters(numOfItersN, N_EQ_FILE_PATH);
}

// Процедура проверки систем уравнений методом Ньютона
template<typename Type>
void checkTestSystem(Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), Type x0, Type y0,
const std::string &N_SYS_FILE_PATH, Type accuracy = 1e-6){
    std::vector<Type> solution;
    std::size_t numOfIterations = getSystemSolutionNewthon(solution, f1, f2, x0, y0, accuracy);
    writeSysResNewthon(solution, x0, y0, accuracy, N_SYS_FILE_PATH);
}

template<typename Type>
void temp_main(){
    Type accuracy = 1e-6;
    std::size_t stopIt = 10000;
    Type firstX, lastX, x0, y0, h;

    // Уравнения
    firstX = 0.0;
    lastX = 1.0;
    x0 = 0.0;
    h = 1e-4;
    checkTestEquations(func1, firstX, lastX, x0, B_EQ_FILE_PATH_1, N_EQ_FILE_PATH_1, accuracy, h, stopIt);

    firstX = -1.0;
    lastX = 10.0;
    x0 = 1.0;
    h = 1e-6;
    checkTestEquations(func2, firstX, lastX, x0, B_EQ_FILE_PATH_2, N_EQ_FILE_PATH_2, accuracy, h, stopIt);

    firstX = 0.0;
    lastX = 1.0;
    x0 = 0.75;
    h = 1e-8;
    checkTestEquations(func3, firstX, lastX, x0, B_EQ_FILE_PATH_3, N_EQ_FILE_PATH_3, accuracy, h, stopIt);

    // Cистемы
    x0 = -2.0;
    y0 = 7.0;
    checkTestSystem(func41, func42, x0, y0, N_SYS_FILE_PATH_1, accuracy);

    x0 = -2.0;
    y0 = 7.0;
    checkTestSystem(func51, func52, x0, y0, N_SYS_FILE_PATH_2, accuracy);
}

int main(){
    temp_main<double>();

    double accuracy = 1e-6;
    std::vector<double> solution;
    std::vector<double> xGrid;
    std::vector<double> yGrid;
    getRectangularGrid(xGrid, yGrid, 10.0, 10.0, 4);
    writeVectorFile(xGrid, "D:\\Calc_Methods\\Lab5\\xGrid.txt");
    writeVectorFile(yGrid, "D:\\Calc_Methods\\Lab5\\yGrid.txt");
    std::vector<double> iterationVec(3);
    std::vector<std::vector<double>> iterationMatrix;
    for (std::size_t i = 0; i < xGrid.size(); i++){
        for (std::size_t j = 0; j < yGrid.size(); j++){
            iterationVec[0] = xGrid[i];
            iterationVec[1] = yGrid[j];
            iterationVec[2] = static_cast<double>(getSystemSolutionNewthon(solution, func41, func42, xGrid[i], yGrid[j], accuracy));
            iterationMatrix.push_back(iterationVec);
        }
    }
    writeMatrixFile(iterationMatrix, "D:\\Calc_Methods\\Lab5\\iterationData.txt");
    return 0;
}