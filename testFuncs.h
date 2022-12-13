// Файл, содержащий функции для тестов
#ifndef TEST_FUNC_H
#define TEST_FUNC_H

#include<cmath>

// Функции для тестов уравнений
template<typename Type>
Type func1(Type x){
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

template<typename Type>
Type func2(Type x){
    return std::sqrt(x + 1.0) - 1.0;
}

template<typename Type>
Type func3(Type x){
    return 35.0 * std::pow(x, 3) - 67.0 * std::pow(x, 2) - 3.0 * x + 3.0;
}

// Тесты для системы

template<typename Type>
Type func41(Type x, Type y){
    return x * x - y * y - 15.0;
}

template<typename Type>
Type func42(Type x, Type y){
    return x * y + 4.0;
}

template<typename Type>
Type func51(Type x, Type y){
    return x * x + y * y + x + y - 8.0;
}

template<typename Type>
Type func52(Type x, Type y){
    return x * x + y * y + x * y - 7.0;
}

template<typename Type>
Type func61(Type x, Type y){
    return -1.0 + std::pow(x, 4) - 6 * std::pow(x, 2) * std::pow(y, 2) + std::pow(y, 4);
}

template<typename Type>
Type func62(Type x, Type y){
    return 4.0 * (std::pow(x, 3) * y - std::pow(y, 3) * x);
}

// Аналитически вычисленные элементы матрицы Якоби
template<typename Type>
std::vector<Type> getJacobiElems1(Type x, Type y){
    return std::vector<Type>{
        2.0 * x, -2.0 * y,
        y, x
    };
}

template<typename Type>
std::vector<Type> getJacobiElems2(Type x, Type y){
    return std::vector<Type>{
        2.0 * x + 1.0, 2.0 * y + 1.0,
        2.0 * x + y, 2.0 * y + x
    };
}
#endif