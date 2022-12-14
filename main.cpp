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
const std::string &B_EQ_FILE_PATH, const std::string &N_EQ_FILE_PATH, Type accuracy = 1e-6, Type h = 1e-4, std::size_t stopIt = 10000, std::size_t n = 51){
    // Локализация корней
    std::vector<std::vector<Type>> segmentMatrix;
    std::size_t numOfRoots = locoliseRoots(f, firstX, lastX, n, segmentMatrix);

    // Метод биссекции
    if (numOfRoots != 0){
        Type solBisection = getEquationSolutionBissection(f, segmentMatrix[0][0], segmentMatrix[0][1], accuracy, stopIt);
        std::size_t numOfItersB = getIterationsBissection(f, segmentMatrix[0][0], segmentMatrix[0][1], accuracy, stopIt);
        writeEqResBissection(solBisection, segmentMatrix[0][0], segmentMatrix[0][1], accuracy, B_EQ_FILE_PATH);
        writeIters(numOfItersB, B_EQ_FILE_PATH);
        for (std::size_t i = 1; i < numOfRoots; i++){
            solBisection = getEquationSolutionBissection(f, segmentMatrix[i][0], segmentMatrix[i][1], accuracy, stopIt);
            numOfItersB = getIterationsBissection(f, segmentMatrix[i][0], segmentMatrix[i][1], accuracy, stopIt);
            writeEqResBissection(solBisection, segmentMatrix[i][0], segmentMatrix[i][1], accuracy, B_EQ_FILE_PATH, true);
            writeIters(numOfItersB, B_EQ_FILE_PATH);
        }    
    }

    // Метод Ньютона
    if (numOfRoots != 0){
        Type x0 = (segmentMatrix[0][0] + segmentMatrix[0][1]) / 2.0;
        Type solNewthon = getEquationSolutionNewthon(f, x0, h, accuracy, stopIt);
        std::size_t numOfItersN = getIterationsNewthon(f, x0, h, accuracy, stopIt);
        writeEqResNewthon(solNewthon, x0, accuracy, N_EQ_FILE_PATH);
        writeIters(numOfItersN, N_EQ_FILE_PATH);
        solNewthon = getEquationSolutionNewthonModified(f, x0, h, accuracy, stopIt);
        numOfItersN = getIterationsNewthonModified(f, x0, h, accuracy, stopIt);
        writeEqResNewthon(solNewthon, x0, accuracy, N_EQ_FILE_PATH, true);
        writeIters(numOfItersN, N_EQ_FILE_PATH);
        for (std::size_t i = 1; i < numOfRoots; i++){
            x0 = (segmentMatrix[i][0] + segmentMatrix[i][1]) / 2.0;
            solNewthon = getEquationSolutionNewthon(f, x0, h, accuracy, stopIt);
            numOfItersN = getIterationsNewthon(f, x0, h, accuracy, stopIt);
            writeEqResNewthon(solNewthon, x0, accuracy, N_EQ_FILE_PATH, true);
            writeIters(numOfItersN, N_EQ_FILE_PATH);
            solNewthon = getEquationSolutionNewthonModified(f, x0, h, accuracy, stopIt);
            numOfItersN = getIterationsNewthonModified(f, x0, h, accuracy, stopIt);
            writeEqResNewthon(solNewthon, x0, accuracy, N_EQ_FILE_PATH, true);
            writeIters(numOfItersN, N_EQ_FILE_PATH);
        }
    }
    
}

// Процедура проверки систем уравнений методом Ньютона
template<typename Type>
void checkTestSystem(Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), std::vector<Type> (*getJacobiMatrixElems)(Type x, Type y), 
Type x0, Type y0, Type L1, Type L2, std::size_t N, const std::string &N_SYS_FILE_PATH, const std::string &A_N_SYS_FILE_PATH, 
const std::string &IT_SYS_FILE_PATH, const std::string &A_IT_SYS_FILE_PATH, Type h = 1e-4, Type accuracy = 1e-6, std::size_t stopIt = 10000){
    // Решение системы
    std::vector<Type> solution;
    std::size_t numOfIterations = getSystemSolutionNewthon(solution, f1, f2, x0, y0, accuracy, h, 2.0, stopIt);
    writeSysResNewthon(solution, x0, y0, accuracy, N_SYS_FILE_PATH);
    writeIters(numOfIterations, N_SYS_FILE_PATH);
    numOfIterations = getSystemSolutionNewthonAnalytic(solution, getJacobiMatrixElems, f1, f2, x0, y0, accuracy, 2.0, stopIt);
    writeSysResNewthon(solution, x0, y0, accuracy, A_N_SYS_FILE_PATH);
    writeIters(numOfIterations, A_N_SYS_FILE_PATH);
    
    // Анализ сходимости
    std::vector<Type> xGrid;
    std::vector<Type> yGrid;
    Type h1 = 2.0 * L1 / N;
    Type h2 = 2.0 * L2 / N;
    for (std::size_t i = 0; i < N + 1; i++){
        xGrid.push_back(-L1 + i * h1);
        yGrid.push_back(-L2 + i * h2);
    }
    std::vector<Type> iterationVec(3);
    std::vector<std::vector<Type>> iterationMatrix, iterationMatrixAnalytic;
    for (std::size_t i = 0; i < xGrid.size(); i++){
        for (std::size_t j = 0; j < yGrid.size(); j++){
            iterationVec[0] = xGrid[i];
            iterationVec[1] = yGrid[j];
            iterationVec[2] = static_cast<Type>(getSystemSolutionNewthon(solution, f1, f2, xGrid[i], yGrid[j], h, accuracy, 2.0, stopIt));
            iterationMatrix.push_back(iterationVec);
        }
    }
    writeMatrixFile(iterationMatrix, IT_SYS_FILE_PATH);
    for (std::size_t i = 0; i < xGrid.size(); i++){
        for (std::size_t j = 0; j < yGrid.size(); j++){
            iterationVec[0] = xGrid[i];
            iterationVec[1] = yGrid[j];
            iterationVec[2] = static_cast<Type>(getSystemSolutionNewthonAnalytic(solution, getJacobiMatrixElems, f1, f2, xGrid[i], yGrid[j], accuracy, 2.0, stopIt));
            iterationMatrixAnalytic.push_back(iterationVec);
        }
    }
    writeMatrixFile(iterationMatrixAnalytic, A_IT_SYS_FILE_PATH);
}

template<typename Type>
void temp_main(){
    Type accuracy = 1e-6;
    std::size_t stopIt = 10000;
    Type firstX, lastX, x0, y0, h, L1, L2, N;
    std::size_t n;

    // Уравнения
    firstX = 0.0;
    lastX = 1.0;
    x0 = 0.0;
    h = 1e-4;
    n = 51;
    checkTestEquations(func1, firstX, lastX, x0, B_EQ_FILE_PATH_1, N_EQ_FILE_PATH_1, accuracy, h, stopIt, n);

    firstX = -1.0;
    lastX = 10.0;
    x0 = 1.0;
    h = 1e-6;
    n = 51;
    checkTestEquations(func2, firstX, lastX, x0, B_EQ_FILE_PATH_2, N_EQ_FILE_PATH_2, accuracy, h, stopIt, n);

    firstX = 0.0;
    lastX = 1.0;
    x0 = 0.75;
    h = 1e-8;
    n = 51;
    checkTestEquations(func3, firstX, lastX, x0, B_EQ_FILE_PATH_3, N_EQ_FILE_PATH_3, accuracy, h, stopIt, n);

    // Cистемы
    x0 = -2.0;
    y0 = 7.0;
    L1 = 10.0;
    L2 = 10.0;
    N = 50;
    h = 1e-4;
    checkTestSystem(func41, func42, getJacobiElems1, x0, y0, L1, L2, N, N_SYS_FILE_PATH_1, A_N_SYS_FILE_PATH_1, IT_SYS_FILE_PATH_1, A_IT_SYS_FILE_PATH_1, h, accuracy, stopIt);

    x0 = -2.0;
    y0 = 7.0;
    L1 = 10.0;
    L2 = 10.0;
    N = 50;
    h = 1e-4;
    checkTestSystem(func51, func52, getJacobiElems2, x0, y0, L1, L2, N, N_SYS_FILE_PATH_2, A_N_SYS_FILE_PATH_2, IT_SYS_FILE_PATH_2, A_IT_SYS_FILE_PATH_2, h, accuracy, stopIt);

    // Бассейн Ньютона
    Type R = 2;
    n = 200;
    h = 1e-4;
    //writeNewthonSwPool(func61, func62, R, n, h, accuracy, "D:\\Calc_Methods\\Lab5\\NewthonSwPool.txt", stopIt);
}

int main(){
    //temp_main<double>();
    std::cout << getEquationSolutionNewthon(func1, 1.0) << '\n';
    std::cout << getConvergEstimateNewthon(func1, 1.0, 0.75, 1e-4, 1e-6, 10000) << '\n'; // Функция из методы
    std::cout << getConvergEstimateNewthon(func4, 100.0, 1.0, 1e-4, 1e-6, 10000) << '\n'; // x^2
    std::cout << getConvergEstimateNewthon(func6, 1.0, 0.0, 1e-4, 1e-6, 10000) << '\n'; // sin
    std::cout << getConvergEstimateNewthon(func5, 2.0, 1.0, 1e-4, 1e-6, 10000) << '\n'; // корень


    return 0;
}