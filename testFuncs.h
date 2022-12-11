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

#endif