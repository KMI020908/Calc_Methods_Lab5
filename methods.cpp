#include "methods.h"

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    for (std::size_t k = 0; k < rows; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        for (std::size_t i = k + 1; i < rows; i++){   // Частичный выбор главного элемента
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainRow != k){ //Замена строк       
            Type temp;
            for (std::size_t j = 0; j < cols; j++){
                temp = lCoefs[k][j];
                lCoefs[k][j] = lCoefs[mainRow][j];
                lCoefs[mainRow][j] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        for (std::size_t i = k + 1; i < rows; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < cols; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (std::abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG gaussMethodFull(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    std::vector<std::size_t> switchCols; // Вектор, хранящий индексы перемещенных столбцов исходной матрицы
    for (std::size_t k = 0; k < rows; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        std::size_t mainСol = k; // Столбец главного элемента
        // Полный выбор
        // Поиск главного элмента в k-ом столбце  
        for (std::size_t i = k + 1; i < rows; i++){ 
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainRow != k){ //Замена строк
            Type temp;
            for (std::size_t j = 0; j < cols; j++){
                temp = lCoefs[k][j];
                lCoefs[k][j] = lCoefs[mainRow][j];
                lCoefs[mainRow][j] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        // Поиск главного элмента в k-ой строке 
        for (std::size_t j = k + 1; j < cols; j++){ 
            if (std::abs(lCoefs[k][j]) > std::abs(mainValue)){
                mainValue = lCoefs[k][j]; 
                mainСol = j;
            }
        }
        //Замена столбцов
        if (mainСol != k){ 
            Type temp;
            for (std::size_t i = 0; i < rows; i++){
                temp = lCoefs[i][k];
                lCoefs[i][k] = lCoefs[i][mainСol];
                lCoefs[i][mainСol] = temp;
            }
            switchCols.push_back(k);
            switchCols.push_back(mainСol);
        }

        // Прямой ход Гаусса 
        for (std::size_t i = k + 1; i < rows; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < cols; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (std::abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    // Обратная перестановка
    for (int i = switchCols.size() - 2; i >= 0; i -= 2){
        Type temp = solution[switchCols[i]];
        solution[switchCols[i]] = solution[switchCols[i + 1]];
        solution[switchCols[i + 1]] = temp;

    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG tridiagonalAlgoritm(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, std::vector<Type> &solution){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    solution.resize(rows, 0.0);
    std::vector<Type> alpha, beta;

    alpha.push_back(-lCoefs[0][1] / lCoefs[0][0]);
    beta.push_back(rCoefs[0] / lCoefs[0][0]);

    for (std::size_t i = 1; i < rows - 1; i++){
        Type coef = lCoefs[i][i - 1] * alpha[i - 1] + lCoefs[i][i];
        alpha.push_back(-lCoefs[i][i + 1] / coef);
        beta.push_back((rCoefs[i] - lCoefs[i][i - 1] * beta[i - 1]) / coef);
    }
    beta.push_back(
        (rCoefs[rows - 1] - lCoefs[rows - 1][rows - 2] * beta[rows - 2]) / (lCoefs[rows - 1][rows - 2] * alpha[rows - 2] + lCoefs[rows - 1][rows - 1])
    );
    solution[rows - 1] = beta[rows - 1];
    for (int i = rows - 2; i >= 0; i--){
        solution[i] = alpha[i] * solution[i + 1] + beta[i];
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG tridiagonalAlgoritm(const std::vector<Type> &diag, std::vector<Type> &lDiag, std::vector<Type> &uDiag, const std::vector<Type> &rVec, std::vector<Type> &solution){
    std::size_t dim = diag.size();
    if (lDiag.size() != dim - 1 || uDiag.size() != dim - 1 || rVec.size() != dim){
        return NO_SOLUTION;
    }
    lDiag.insert(lDiag.begin(), 1, 0.0);
    uDiag.push_back(0.0);
    solution.resize(dim);
    std::vector<Type> alpha = {-uDiag[0] / diag[0]};
    std::vector<Type> beta = {rVec[0] / diag[0]};
    for (std::size_t i = 1; i < dim - 1; i++){
        Type c = diag[i] + lDiag[i] * alpha[i - 1];
        alpha.push_back(-uDiag[i] / c); 
        beta.push_back((rVec[i] - lDiag[i] * beta[i - 1]) / c);
    }
    beta.push_back((rVec[dim - 1] - lDiag[dim - 1] * beta[dim - 2]) / (diag[dim - 1] + lDiag[dim - 1] * alpha[dim - 2]));
    solution[dim - 1] = beta[dim - 1];
    for (int i = dim - 2; i >= 0; i--){
        solution[i] = alpha[i] * solution[i + 1] + beta[i];
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    if (rows != cols){
        solution.resize(0);
        return NO_SOLUTION;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/std::sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                Type s = lCoefs[i][k]/std::sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                for (std::size_t j = k; j < cols; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c * lCoefs[k][j] + s * lCoefs[i][j];
                    lCoefs[i][j] = -s * temp + c * lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < accuracy)
                        lCoefs[i][j] = 0;
                }
                Type temp = rCoefs[k];
                rCoefs[k] = c*rCoefs[k] + s*rCoefs[i];
                rCoefs[i] = -s*temp + c*rCoefs[i];
            }
        }
    }
    if (std::abs(lCoefs[rows - 1][rows - 1]) < accuracy){  // detA = 0
        return NO_SOLUTION;
    }
    
     // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }   
    return HAS_SOLUTION;
}

template<typename Type>
QUADRATIC_FLAG findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t rows = lCoefs.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/std::sqrt(lCoefs[k][k] * lCoefs[k][k] + lCoefs[i][k] * lCoefs[i][k]);
                Type s = lCoefs[i][k]/std::sqrt(lCoefs[k][k] * lCoefs[k][k] + lCoefs[i][k] * lCoefs[i][k]);
                for (std::size_t j = 0; j < cols; j++){
                    Type temp = Q[k][j];
                    Q[k][j] = c * Q[k][j] + s * Q[i][j];
                    Q[i][j] = -s * temp + c * Q[i][j];
                    if (std::abs(Q[i][j]) < accuracy)
                        Q[i][j] = 0.0;
                }
                for (std::size_t j = k; j < cols; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c * lCoefs[k][j] + s * lCoefs[i][j];
                    lCoefs[i][j] = -s * temp + c * lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < accuracy)
                        lCoefs[i][j] = 0.0;
                }
            }
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findQBlock(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, std::size_t rows, Type accuracy){
    std::size_t cols = 0;
    if (rows != 0)
        cols = rows;
    else
        return NOT_QUADRATIC;
    if (lCoefs.size() < rows){
        return NOT_QUADRATIC;
    }
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/std::sqrt(lCoefs[k][k] * lCoefs[k][k] + lCoefs[i][k] * lCoefs[i][k]);
                Type s = lCoefs[i][k]/std::sqrt(lCoefs[k][k] * lCoefs[k][k] + lCoefs[i][k] * lCoefs[i][k]);
                for (std::size_t j = 0; j < cols; j++){
                    Type temp = Q[k][j];
                    Q[k][j] = c * Q[k][j] + s * Q[i][j];
                    Q[i][j] = -s * temp + c * Q[i][j];
                    if (std::abs(Q[i][j]) < accuracy)
                        Q[i][j] = 0.0;
                }
                for (std::size_t j = k; j < cols; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c * lCoefs[k][j] + s * lCoefs[i][j];
                    lCoefs[i][j] = -s * temp + c * lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < accuracy)
                        lCoefs[i][j] = 0.0;
                }
            }
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findQMatrix3Diag(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    // k < rows - 2
    for (std::size_t k = 0; k < rows - 2; k++){
        if (std::abs(matrix[k + 1][k]) >= accuracy){
            Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            for (std::size_t j = 0; j < cols; j++){
                Type temp = Q[k][j];
                Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
                Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
                if (std::abs(Q[k + 1][j]) < accuracy)
                    Q[k + 1][j] = 0.0;
            }
            for (std::size_t j = k; j < k + 3; j++){
                Type temp = matrix[k][j];
                matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
                matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
                if (std::abs(matrix[k + 1][j]) < accuracy)
                    matrix[k + 1][j] = 0.0;
            }
        }
    }
    // k = rows - 2
    std::size_t k = rows - 2;
    if (std::abs(matrix[k + 1][k]) >= accuracy){
        Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
        Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
        for (std::size_t j = 0; j < cols; j++){
            Type temp = Q[k][j];
            Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
            Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
            if (std::abs(Q[k + 1][j]) < accuracy)
                Q[k + 1][j] = 0.0;
        }
        for (std::size_t j = k; j < cols; j++){
            Type temp = matrix[k][j];
            matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
            matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
            if (std::abs(matrix[k + 1][j]) < accuracy)
                matrix[k + 1][j] = 0.0;
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findQBlock3Diag(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, std::size_t rows, Type accuracy){
    std::size_t cols = 0;
    if (rows != 0)
        cols = rows;
    else
        return NOT_QUADRATIC;
    if (matrix.size() < rows){
        return NOT_QUADRATIC;
    }    
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    // k < rows - 2
    for (std::size_t k = 0; k < rows - 2; k++){
        if (std::abs(matrix[k + 1][k]) >= accuracy){
            Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            for (std::size_t j = 0; j < cols; j++){
                Type temp = Q[k][j];
                Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
                Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
                if (std::abs(Q[k + 1][j]) < accuracy)
                    Q[k + 1][j] = 0.0;
            }
            for (std::size_t j = k; j < k + 3; j++){
                Type temp = matrix[k][j];
                matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
                matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
                if (std::abs(matrix[k + 1][j]) < accuracy)
                    matrix[k + 1][j] = 0.0;
            }
        }
    }
    // k = rows - 2
    std::size_t k = rows - 2;
    if (std::abs(matrix[k + 1][k]) >= accuracy){
        Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
        Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
        for (std::size_t j = 0; j < cols; j++){
            Type temp = Q[k][j];
            Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
            Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
            if (std::abs(Q[k + 1][j]) < accuracy)
                Q[k + 1][j] = 0.0;
        }
        for (std::size_t j = k; j < cols; j++){
            Type temp = matrix[k][j];
            matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
            matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
            if (std::abs(matrix[k + 1][j]) < accuracy)
                matrix[k + 1][j] = 0.0;
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findQMatrixHess(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    for (std::size_t k = 0; k < rows - 1; k++){
        if (std::abs(matrix[k + 1][k]) >= accuracy){
            Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            for (std::size_t j = 0; j < cols; j++){
                Type temp = Q[k][j];
                Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
                Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
                if (std::abs(Q[k + 1][j]) < accuracy)
                    Q[k + 1][j] = 0.0;
            }
            for (std::size_t j = k; j < cols; j++){
                Type temp = matrix[k][j];
                matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
                matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
                if (std::abs(matrix[k + 1][j]) < accuracy)
                    matrix[k + 1][j] = 0.0;
            }
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findQBlockHess(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, std::size_t rows, Type accuracy){
    std::size_t cols = 0;
    if (rows != 0)
        cols = rows;
    else
        return NOT_QUADRATIC;
    if (matrix.size() < rows){
        return NOT_QUADRATIC;
    }
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    for (std::size_t k = 0; k < rows - 1; k++){
        if (std::abs(matrix[k + 1][k]) >= accuracy){
            Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            for (std::size_t j = 0; j < cols; j++){
                Type temp = Q[k][j];
                Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
                Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
                if (std::abs(Q[k + 1][j]) < accuracy)
                    Q[k + 1][j] = 0.0;
            }
            for (std::size_t j = k; j < cols; j++){
                Type temp = matrix[k][j];
                matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
                matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
                if (std::abs(matrix[k + 1][j]) < accuracy)
                    matrix[k + 1][j] = 0.0;
            }
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
std::size_t transposeMatrix(std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return 0;
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = i + 1; j < cols; j++){
            Type temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
    return rows;
}

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    std::vector<Type> b1(rows);  // Правая часть после подстановки полученного решения
        for (std::size_t i = 0; i < rows; i++){
                Type sum = 0.0;
                for (std::size_t k = 0; k < cols; k++){
                    sum += lCoefs[i][k] * solution[k];
                }
                b1[i] = sum;
        }
        std::vector<Type> discrepancyVector(rows); // Вектор невязки
        for (std::size_t i = 0; i < rows; i++){
            discrepancyVector[i] = rCoefs[i] - b1[i];
        }
        Type discrepancy = 0.0; // Невязка
        for (std::size_t i = 0; i < rows; i++){
            discrepancy += discrepancyVector[i] * discrepancyVector[i];     
        }
    return std::sqrt(discrepancy);
}

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    Type norm1OfMatrix = 0;
    for (std::size_t j = 0; j < cols; j++){
        Type sum = 0.0;
        for (std::size_t i = 0; i < rows; i++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > norm1OfMatrix)
            norm1OfMatrix = sum;
    }
    return norm1OfMatrix;
}

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<std::vector<Type>> invMatrix; //Обратная к A матрица
    INVERTIBLE_FLAG flag = invertMatrix(matrix, invMatrix);
    Type norm1OfMatrix = findMatrixNorm1(matrix);
    Type norm1OfInvMatrix = 0;
    if (flag == IS_INVERTIBLE)
        norm1OfInvMatrix = findMatrixNorm1(invMatrix);
    else
        norm1OfInvMatrix = INFINITY;
    Type cond = norm1OfMatrix * norm1OfInvMatrix;
    return cond;
}

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    Type normInfOfMatrix = 0;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > normInfOfMatrix)
            normInfOfMatrix = sum;
    }
    return normInfOfMatrix;
}

template<typename Type>
Type normOfMatrix(const std::vector<std::vector<Type>> &matrix, double p){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (p == 1.0){
        return findMatrixNorm1(matrix);
    }
    if (p == INFINITY){
        return findMatrixNormInf(matrix);
    }
    return NAN;
}

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<std::vector<Type>> invMatrix; //Обратная к A матрица
    INVERTIBLE_FLAG flag = invertMatrix(matrix, invMatrix);
    Type normInfOfMatrix = findMatrixNormInf(matrix);
    Type normInfOfInvMatrix = 0;
    if (flag == IS_INVERTIBLE)
        normInfOfInvMatrix = findMatrixNormInf(invMatrix);
    else
        normInfOfInvMatrix = INFINITY;
    Type cond = normInfOfMatrix * normInfOfInvMatrix;
    return cond;
}

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix,
    SOLUTION_FLAG (*method)(std::vector<std::vector<Type>> &, std::vector<Type>&, std::vector<Type>&, Type accuracy)){
    std::size_t rows = inputMatrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = inputMatrix[0].size();
    else
        return NOT_INVERTIBLE;
    if (rows != cols)
        return NOT_INVERTIBLE;
    resMatrix.resize(rows);
    std::vector<std::vector<Type>> E(rows);
    for (std::size_t i = 0; i < rows; i++){
        E[i].resize(cols, 0);
        resMatrix[i].resize(cols);
    }
    for (std::size_t i = 0; i < rows; i++){
        E[i][i] = 1;
    }
    std::vector<Type> solution(rows);
    std::vector<std::vector<Type>> tempMatrix(rows);
    SOLUTION_FLAG flag;
    for (std::size_t i = 0; i < rows; i++){
        tempMatrix = inputMatrix;
        flag = (*method)(tempMatrix, E[i], solution, 1e-7);
        if (flag == NO_SOLUTION){
            for (std::size_t i = 0; i < rows; i++)
                resMatrix[i].clear();
            resMatrix.clear();
            return NOT_INVERTIBLE;
        }
        for (std::size_t k = 0; k < rows; k++)
            resMatrix[k][i] = solution[k];
    }
    return IS_INVERTIBLE;
}

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type sum = 0;
    for (std::size_t i = 0; i < vector.size(); i++)
        sum += std::abs(vector[i]);
    return sum;
}

template<typename Type>
Type norm2OfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type sum = 0;
    for (std::size_t i = 0; i < vector.size(); i++){
        sum += std::pow(vector[i], 2);
    }
    return std::sqrt(sum);
}

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type max = std::abs(vector[0]);
    for (std::size_t i = 1; i < vector.size(); i++)
        if (std::abs(vector[i]) > max)
            max = std::abs(vector[i]);
    return max;
}

template<typename Type>
Type normOfVector(const std::vector<Type> &vector, double p){
    if (!vector.size()){
        return NAN;
    }
    if (p == 2.0){
        return norm2OfVector(vector);
    }
    if (p == 1.0){
        return norm1OfVector(vector);
    }
    if (p == INFINITY){
        return normInfOfVector(vector);
    }
    return NAN;
}

template<typename Type>
Type findLowerBoundOfcond1(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1, Type delta2, Type delta3){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<Type> rCoefs1(rows), rCoefs2(rows), rCoefs3(rows);
    for (std::size_t i = 0; i < rows; i++){
        rCoefs1[i] = rCoefs[i] + delta1;
        rCoefs2[i] = rCoefs[i] + delta2;
        rCoefs3[i] = rCoefs[i] + delta3;
    }
    std::vector<std::vector<Type>> Q;
    findQMatrix(lCoefs, Q);
    transposeMatrix(Q);

    std::vector<Type> solution(rows), perturbSolution(rows), deltaSolution(rows);
    std::vector<Type> tempRCoefs;

    multiplyMatrix(Q, rCoefs, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, solution);

    multiplyMatrix(Q, rCoefs1, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs1(rows, delta1);
    Type lowerBound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs1));

    multiplyMatrix(Q, rCoefs2, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs2(rows, delta2);
    Type bound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs2));
    if (bound > lowerBound)
        lowerBound = bound;

    multiplyMatrix(Q, rCoefs3, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs3(rows, delta3);
    bound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs3));
    if (bound > lowerBound)
        lowerBound = bound;
    
    return lowerBound;
}

template<typename Type>
Type findLowerBoundOfcondInf(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1, Type delta2, Type delta3){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<Type> rCoefs1(rows), rCoefs2(rows), rCoefs3(rows);
    for (std::size_t i = 0; i < rows; i++){
        rCoefs1[i] = rCoefs[i] + delta1;
        rCoefs2[i] = rCoefs[i] + delta2;
        rCoefs3[i] = rCoefs[i] + delta3;
    }
    std::vector<std::vector<Type>> Q;
    findQMatrix(lCoefs, Q);
    transposeMatrix(Q);

    std::vector<Type> solution(rows), perturbSolution(rows), deltaSolution(rows);
    std::vector<Type> tempRCoefs;

    multiplyMatrix(Q, rCoefs, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, solution);

    multiplyMatrix(Q, rCoefs1, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs1(rows, delta1);
    Type lowerBound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs1));

    multiplyMatrix(Q, rCoefs2, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs2(rows, delta2);
    Type bound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs2));
    if (bound > lowerBound)
        lowerBound = bound;

    multiplyMatrix(Q, rCoefs3, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs3(rows, delta3);
    bound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs3));
    if (bound > lowerBound)
        lowerBound = bound;
    
    return lowerBound;
}

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec, std::vector<Type> &result){
    std::size_t rows1 = matrix.size();
    std::size_t cols = 0;
    if (rows1 != 0)
        cols = matrix[0].size();
    else
        return NOT_MULTIPLIED;
    std::size_t rows2 = vec.size();
    if (cols != rows2)
        return NOT_MULTIPLIED;
    result.resize(rows1);
    std::size_t rows = result.size();
    for (size_t i = 0; i < rows; i++){
        Type sum = 0;
        for (size_t k = 0; k < cols; k++){
            sum += matrix[i][k] * vec[k];
        }
        result[i] = sum;
    }
    return IS_MULTIPLIED;
}

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else{
        os << 0;
        return os;
    }
    for (std::size_t i = 0; i < rows - 1; i++){
        for (std::size_t j = 0; j < cols; j++){
            os << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
    for (std::size_t j = 0; j < cols; j++){
            os << matrix[rows - 1][j] << ' ';
        }
    return os;
}

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<Type> &vector){
    std::size_t rows = vector.size();
    if (!rows){
        os << 0;
        return os;
    }  
    os << "{ ";
    for (std::size_t i = 0; i < rows - 1; i++)
        os << vector[i] << ", ";
    os << vector[rows - 1] << ' ';
    os << '}';
    return os;
}

template<typename Type>
std::vector<Type> operator+(const std::vector<Type>& vec1, const std::vector<Type>& vec2)
{
    auto result = vec1;
    const std::size_t size_max = std::max<std::size_t>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (std::size_t i = 0; i < result.size() && i < vec2.size(); i++){   
        result[i] += vec2[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator-(const std::vector<Type>& vec1, const std::vector<Type>& vec2)
{
    auto result = vec1;
    const std::size_t size_max = std::max<std::size_t>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (std::size_t i = 0; i < result.size() && i < vec2.size(); i++){   
        result[i] -= vec2[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec)
{
    std::size_t rows1 = matrix.size();
    std::size_t cols = 0;
    if (rows1 != 0)
        cols = matrix[0].size();
    else
        return std::vector<double>(1, NAN);
    std::size_t rows2 = vec.size();
    if (cols != rows2)
        return std::vector<double>(rows2, NAN);
    std::vector<Type> result(rows1);
    for (std::size_t i = 0; i < rows1; i++){
        Type sum = 0;
        for (std::size_t k = 0; k < cols; k++){
            sum += matrix[i][k] * vec[k];
        }
        result[i] = sum;
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(Type num, const std::vector<Type> &vec){
    std::size_t size = vec.size();
    if (!size)
        return std::vector<double>(1, NAN);
    std::vector<Type> res(size);
    for (std::size_t i = 0; i < size; i++){
        res[i] = num * vec[i];
    }
    return res;
}

template<typename Type>
std::vector<std::vector<Type>> operator*(Type num, const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (!rows)
        return std::vector<std::vector<Type>>(1, std::vector<Type>(1, NAN));
    else
        cols = matrix[0].size();
    std::vector<std::vector<Type>> res;
    std::vector<Type> tempVec(cols);
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++)
        {
            tempVec[j] = num * matrix[i][j];
        }
        res.push_back(tempVec);
    }
    return res;
}

template<typename Type>
Type dot(const std::vector<Type> &v1, const std::vector<Type> &v2){
    std::size_t size1 = v1.size();
    std::size_t size2 = v2.size();
    std::size_t minSize = 0;
    if (size1 < size2){
        minSize = size1;
    }
    else{
        minSize = size2;
    }
    Type sum = 0.0;
    for (std::size_t i = 0; i < minSize; i++){
        sum += v1[i] * v2[i];
    }
    return sum;
}

template<typename Type>
std::size_t simpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type tao, Type accuracy, double p, Type epsilon_0, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t JacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy, double p, Type epsilon_0, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t relaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy, Type omega, double p, Type epsilon_0, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt){
            break;
        }
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t relaxationMethodFor3Diag(const std::vector<Type> &a, const std::vector<Type> &b, const std::vector<Type> &c, const std::vector<Type> &d, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy, Type omega, double p, Type epsilon_0, std::size_t stopIt){
    if (!b.size() || b.size() != d.size() || a.size() != b.size() - 1 || c.size() != a.size())
        return NO_SOLUTION;
    std::size_t dim = b.size();
    solution.resize(dim); // Искомое решение
    std::vector<Type> prev_solution = firstVec;
    solution[0] = (1 - omega) * prev_solution[0] - (omega / b[0]) * (c[0] * prev_solution[1] - d[0]);
    for (std::size_t i = 1; i < dim - 1; i++){    
        solution[i] = (1 - omega) * prev_solution[i] - (omega / b[i]) * (a[i - 1] * solution[i - 1] + c[i] * prev_solution[i + 1] - d[i]);
    }
    solution[dim - 1] = (1 - omega) * prev_solution[dim - 1] - (omega / b[dim - 1]) * (a[dim - 2] * prev_solution[dim - 2] - d[dim - 1]);
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        solution[0] = (1 - omega) * prev_solution[0] - (omega / b[0]) * (c[0] * prev_solution[1] - d[0]);
        for (std::size_t i = 1; i < dim - 1; i++){    
            solution[i] = (1 - omega) * prev_solution[i] - (omega / b[i]) * (a[i - 1] * solution[i - 1] + c[i] * prev_solution[i + 1] - d[i]);
        }
        solution[dim - 1] = (1 - omega) * prev_solution[dim - 1] - (omega / b[dim - 1]) * (a[dim - 2] * prev_solution[dim - 2] - d[dim - 1]);
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}
 
template<typename Type>
QUADRATIC_FLAG findCanonicalFormSimpleIt(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type tao){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    } 
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                C[i][j] = -tao * lCoefs[i][j];
            }
            else{
                C[i][j] = -tao * lCoefs[i][j] + 1;
            }
        }
    }
    for (std::size_t i = 0; i < cols; i++){
        y[i] = tao * rCoefs[i];
    }
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findCanonicalFormJacobi(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    } 
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                C[i][j] = -lCoefs[i][j] / lCoefs[i][i];
            }
            else{
                C[i][j] = 0.0;
            }
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        y[i] = rCoefs[i] / lCoefs[i][i];
    }
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    }
    std::vector<Type> columnOfMatrix(rows, 0);
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * C[j][k];
            }
            Type sum2 = 0.0;
            if (k > i){
                sum2 = lCoefs[i][k];  
                C[i][k] = -(omega / lCoefs[i][i]) * (sum1 + sum2);
            }
            else{
                if (k < i){
                    C[i][k] = -(omega / lCoefs[i][i]) * sum1;
                }
                else{
                    C[i][k] = (1 - omega) - (omega / lCoefs[i][i]) * sum1;
                }
            }
        }
    }
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        Type sum = omega * rCoefs[i] / lCoefs[i][i];
        for (std::size_t k = 0; k < i; k++){
            sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * y[k]; 
        }
        y[i] = sum;
    }
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation2(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    }
    C[0][0] = 1 - omega; 
    for (std::size_t j = 1; j < cols; j++){
        C[0][j] = -omega * lCoefs[0][j] / lCoefs[0][0];
    }
    for (std::size_t i = 1; i < rows; i++){
        for (std::size_t j = 0; j < i; j++){
            Type sum = 0.0;
            for (std::size_t k = 0; k < i; k++){
                sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * C[k][j];
            }
            C[i][j] = sum;
        }
        for (std::size_t j = i; j < cols; j++){
            Type sum = 0.0;
            if (i != j){
                sum = -omega * lCoefs[i][j] / lCoefs[i][i];
            }
            else{
                sum = 1 - omega;
            }
            for (std::size_t k = 0; k < i; k++){
                sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * C[k][j];
            }
            C[i][j] = sum;
        }
    }
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        Type sum = omega * rCoefs[i] / lCoefs[i][i];
        for (std::size_t k = 0; k < i; k++){
            sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * y[k]; 
        }
        y[i] = sum;
    }
    return IS_QUADRATIC;
}

template<typename Type>
Type findLowerBoundOfIterations(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, const std::vector<Type> &firstVec, Type accuracy, ITERATION_METHOD_FLAG method, Type tao, Type omega, double p){
    std::size_t rows = lCoefs.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0.0;
    std::vector<Type> prev_solution = firstVec;
    std::vector<Type> solution(rows);
    std::vector<std::vector<Type>> C;
    std::vector<Type> y;
    switch (method){
        case SIMPLE_IT:
            for (std::size_t i = 0; i < rows; i++){
                Type sum = 0.0;
                for (std::size_t j = 0; j < cols; j++){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
                solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
            }
            findCanonicalFormSimpleIt(lCoefs, rCoefs, C, y, tao);
            break;

        case JACOBI:
            for (std::size_t i = 0; i < rows; i++){
                    Type sum = 0.0;
                    for (std::size_t j = 0; j < cols; j++){
                        if (i != j){
                        sum += lCoefs[i][j] * prev_solution[j];
                        }
                    }
                solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
            }
            findCanonicalFormJacobi(lCoefs, rCoefs, C, y);
            break;

        case RELAXATION:
            for (std::size_t i = 0; i < rows; i++){
                Type sum1 = 0.0;
                for (std::size_t j = 0; j < i; j++){
                    sum1 += lCoefs[i][j] * solution[j];
                }
                Type sum2 = 0.0;
                for (std::size_t j = i + 1; j < cols; j++){
                    sum2 += lCoefs[i][j] * solution[j];
                }
                solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
            }
            findCanonicalFormRelaxation(lCoefs, rCoefs, C, y, omega); 
            break;
            default:
                break;
    }
    Type distance = normOfVector(solution - prev_solution, p);
    Type q = normOfMatrix(C, p);
    return std::ceil(std::log(accuracy * (1 - q) / distance) / std::log(q));
}

template<typename Type>
std::size_t findExactItersSimpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type tao, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (normOfVector(solution - rightSolution, p) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t findExactItersJacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (normOfVector(solution - rightSolution, p) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t findExactRelaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, Type omega, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (normOfVector(solution - rightSolution, p) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt){
            break;
        }
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
Type findNormOfErrAfterEstIt_SIT(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type tao, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (numOfIt < bound + 1){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return normOfVector(solution - rightSolution, p);
}

template<typename Type>
Type findNormOfErrAfterEstIt_JAC(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return INFINITY;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (numOfIt < bound + 1){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return normOfVector(solution - rightSolution, p);
}

template<typename Type>
Type findNormOfErrAfterEstIt_REL(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type omega, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (numOfIt < bound + 1){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt){
            break;
        }
        numOfIt++;
    }
    return normOfVector(solution - rightSolution, p);
}

template<typename Type>
std::size_t findEigenNumsQRMethod(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy, bool hasShift, bool is3Diag){
    std::size_t numOfIters = 0; // Количество итераций
    Type shift = 0.0; // Сдвиг
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return numOfIters;
    eigenList.resize(rows, 0.0);
    std::size_t eigenRow = rows - 1; // Строка в которой по итогу вращений должно получиться собственное значение
    std::vector<std::vector<Type>> Q;
    std::vector<std::vector<Type>> R;
    R.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        R[i].resize(cols, 0.0);
    }
    while (eigenRow != 0){
        numOfIters++;
        // Сдвиг
        if (hasShift){
            shift = matrix[eigenRow][eigenRow];
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] -= shift; 
        }
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                R[i][j] = matrix[i][j];
            }
        }
        if (!is3Diag){
            findQBlock(R, Q, rows, accuracy);
        }
        else{
            findQBlock3Diag(R, Q, rows, accuracy);
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                Type sum = 0.0;
                for (std::size_t k = i; k < rows; k++){
                    sum += R[i][k] * Q[k][j];
                }
                matrix[i][j] = sum;
            }
        }
        // Обратный сдвиг
        if (hasShift){
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] += shift; 
            }
        }
        Type sumOfEigenRow = 0.0;
        for (std::size_t j = 0; j < eigenRow; j++){
            sumOfEigenRow += std::abs(matrix[eigenRow][j]);
        }
        if (sumOfEigenRow < accuracy){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow]; 
            eigenRow--;
            rows--;
            cols--;
        }
        if (eigenRow == 0){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow];
            break;
        }
    }
    return numOfIters;
}

template<typename Type>
QUADRATIC_FLAG getHessenbergMatrix(std::vector<std::vector<Type>> &matrix, Type accuracy, bool isSymmetric){
    std::size_t rows = matrix.size(); 
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    if (!isSymmetric){
        for (std::size_t k = 1; k < rows; k++){
            for (std::size_t i = k + 1; i < rows; i++){
                if (std::abs(matrix[i][k - 1]) >= accuracy){
                    Type c = matrix[k][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    Type s = matrix[i][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    for (std::size_t j = k - 1; j < cols; j++){
                        Type temp = matrix[k][j];
                        matrix[k][j] = c * matrix[k][j] + s * matrix[i][j];
                        matrix[i][j] = -s * temp + c * matrix[i][j];
                    }
                    for (std::size_t j = k - 1; j < rows; j++){
                        Type temp = matrix[j][k];
                        matrix[j][k] = c * matrix[j][k] + s * matrix[j][i];
                        matrix[j][i] = -s * temp + c * matrix[j][i];
                    }
                }
            }
        }    
    }
    else{
       for (std::size_t k = 1; k < rows; k++){
            for (std::size_t i = k + 1; i < rows; i++){
                if (std::abs(matrix[i][k - 1]) >= accuracy){
                    Type c = matrix[k][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    Type s = matrix[i][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    for (std::size_t j = k - 1; j < cols; j++){
                        Type temp = matrix[k][j];
                        matrix[k][j] = c * matrix[k][j] + s * matrix[i][j];
                        matrix[j][k] = matrix[k][j];
                        matrix[i][j] = -s * temp + c * matrix[i][j];
                        matrix[j][i] = matrix[i][j];
                    }
                }
            }
        } 
    }
    return IS_QUADRATIC;
}

template<typename Type>
std::size_t findEigenNumsQRMethodHessenberg(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy, bool hasShift, bool isSymmetric){
    std::size_t numOfIters = 0; // Количество итераций
    Type shift = 0.0; // Сдвиг
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return numOfIters; 
    eigenList.resize(rows, 0.0);
    std::size_t eigenRow = rows - 1; // Строка в которой по итогу вращений должно получиться собственное значение
    std::vector<std::vector<Type>> Q;
    std::vector<std::vector<Type>> R;
    R.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        R[i].resize(cols, 0.0);
    }
    // Приводим к матрице Хессенберга
    getHessenbergMatrix(matrix, accuracy, isSymmetric);
    while (eigenRow != 0){
        numOfIters++;
        // Сдвиг
        if (hasShift){
            shift = matrix[eigenRow][eigenRow];
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] -= shift; 
        }
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                R[i][j] = matrix[i][j];
            }
        }
        if (isSymmetric){
            findQBlock3Diag(R, Q, rows, accuracy);
        }
        else{
            findQBlockHess(R, Q, rows, accuracy);;
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                Type sum = 0.0;
                for (std::size_t k = i; k < rows; k++){
                    sum += R[i][k] * Q[k][j];
                }
                matrix[i][j] = sum;
            }
        }
        // Обратный сдвиг
        if (hasShift){
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] += shift; 
            }
        }
        if (std::abs(matrix[eigenRow][eigenRow - 1]) < accuracy){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow];
            eigenRow--;
            rows--;
            cols--;
        }
        if (eigenRow == 0){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow];
            break;
        }
    }
    return numOfIters;
}

template<typename Type>
std::size_t invertItersMethod(const std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &eigenMatrix, const std::vector<Type> &startEigenList,
Type accuracy, bool is3Diag){
    std::size_t numOfIters = 0; // Количество итераций
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return numOfIters;
    eigenMatrix.clear();
    std::vector<Type> prevEigenVec(rows);
    std::vector<Type> eigenVec(rows);
    std::vector<std::vector<Type>> lambdaMatrix(rows);
    for (std::size_t i = 0; i < cols; i++){
        lambdaMatrix[i].resize(cols, 0.0);
        for (std::size_t j = 0; j < cols; j++){
            lambdaMatrix[i][j] = matrix[i][j];
        }  
    }
    for (std::size_t i = 0; i < rows; i++){
        eigenVec[0] = 1.0;
        for (std::size_t k = 1; k < rows; k++){
            eigenVec[k] = 0.0;
        }
        //eigenVec = { 0.175032, -0.154008, -0.0218063, -0.972198 };
        for (std::size_t k = 0; k < rows; k++){
            lambdaMatrix[k][k] = matrix[k][k] - startEigenList[i];
        }
        while (normOfVector(eigenVec - prevEigenVec) > accuracy && normOfVector(eigenVec + prevEigenVec) > accuracy){
            for (std::size_t k = 0; k < rows; k++){
                prevEigenVec[k] = eigenVec[k];
            }
            if (!is3Diag){
                qrMethod(lambdaMatrix, prevEigenVec, eigenVec, accuracy);
            }
            else{
                tridiagonalAlgoritm(lambdaMatrix, prevEigenVec, eigenVec);
            }
            Type normOfEigVec = normOfVector(eigenVec);
            for (std::size_t k = 0; k < rows; k++){
                eigenVec[k] /= normOfEigVec;
            }
            numOfIters++;
        }
        eigenMatrix.push_back(eigenVec);
    }
    return numOfIters;
}

template<typename Type>
Type invertItersMethodRayleigh(const std::vector<std::vector<Type>> &matrix, std::vector<Type> &startVec, 
std::vector<Type> &eigenVec, Type accuracy, bool is3Diag){
    std::size_t numOfIters = 0; // Количество итераций
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (startVec.size() != rows){
        startVec.resize(rows, 1.0);
    }
    std::vector<std::vector<Type>> lambdaMatrix(rows);
    for (std::size_t i = 0; i < cols; i++){
        lambdaMatrix[i].resize(cols, 0.0);
        for (std::size_t j = 0; j < cols; j++){
            lambdaMatrix[i][j] = matrix[i][j];
        }  
    }
    eigenVec.resize(rows, 0.0);
    for (std::size_t k = 0; k < rows; k++){
        eigenVec[k] = startVec[k];
    }
    Type normOfEigVec = normOfVector(eigenVec);
    if (std::abs(normOfEigVec - 1.0) > accuracy){
        for (std::size_t k = 0; k < rows; k++){
            eigenVec[k] /= normOfEigVec;
        }
    }
    std::vector<Type> prevEigenVec(rows, 0.0);
    Type eigenApprox = dot(matrix * eigenVec, eigenVec);
    for (std::size_t k = 0; k < rows; k++){
        lambdaMatrix[k][k] = matrix[k][k] - eigenApprox;
    }
    while (normOfVector(eigenVec - prevEigenVec) > accuracy && normOfVector(eigenVec + prevEigenVec) > accuracy){
        for (std::size_t k = 0; k < rows; k++){
            prevEigenVec[k] = eigenVec[k];
        }
        if (!is3Diag){
            qrMethod(lambdaMatrix, prevEigenVec, eigenVec, accuracy);
        }
        else{
            tridiagonalAlgoritm(lambdaMatrix, prevEigenVec, eigenVec);
        }
        normOfEigVec = normOfVector(eigenVec);
        for (std::size_t k = 0; k < rows; k++){
            eigenVec[k] /= normOfEigVec;
        }
        eigenApprox = dot(matrix * prevEigenVec, prevEigenVec);
        for (std::size_t k = 0; k < rows; k++){
            lambdaMatrix[k][k] = matrix[k][k] - eigenApprox;
        }
        numOfIters++;
    }
    return eigenApprox;
}

template<typename Type>
FILE_FLAG writeRayleighSwPool(const std::vector<std::vector<Type>> &matrix, Type step, const std::string& OUT_FILE_PATH, 
Type accuracy, bool is3Diag){
    std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (matrix.size() != 3){
        return IS_CLOSED;
    }
    Type teta = 0.05;
    std::vector<Type> startVec(3, 0.0);
    std::vector<Type> eigVec;
    while (teta <= M_PI){
        Type phi = 0.0;
        while (phi <= 2 * M_PI){
            startVec[0] = sin(teta) * cos(phi);
            startVec[1] = sin(teta) * sin(phi);
            startVec[2] = cos(teta);
            Type lambda = invertItersMethodRayleigh(matrix, startVec, eigVec, accuracy, is3Diag);
            file << startVec[0] << '\t' << startVec[1] << '\t' << startVec[2] << '\t' << lambda << '\n';
            phi += step;
        }
        teta += step;
    }
    file.close();
    return IS_CLOSED;
}

// Лаба 3

template<typename Type>
std::size_t getUniformGrid(std::vector<Type> &xGrid, std::vector<Type> &fGrid, 
Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfFinitElems){
    if (lastX <= firstX){
        return 0;
    }
    xGrid.clear();
    fGrid.clear();
    Type h = (lastX - firstX) / numOfFinitElems;
    for (std::size_t i = 0; i < numOfFinitElems + 1; i++){
        Type tempX = firstX + h * i;
        xGrid.push_back(tempX);
        fGrid.push_back(f(tempX));
    }
    return numOfFinitElems + 1;
}

template<typename Type>
std::size_t getChebyshevGrid(std::vector<Type> &xGrid, std::vector<Type> &fGrid, 
Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfFinitElems){
    if (lastX <= firstX){
        return 0;
    }
    xGrid.clear();
    fGrid.clear();
    for (std::size_t i = 0; i < numOfFinitElems + 1; i++){
        Type tempX = (firstX + lastX) / 2.0 + (firstX - lastX) / 2.0 * std::cos(((2.0 * i + 1.0) * M_PI) / (2.0 * (numOfFinitElems + 1.0))); 
        xGrid.push_back(tempX);
        fGrid.push_back(f(tempX));
    }
    return numOfFinitElems + 1;  
}

// Коэффициенты Лагранжа 
template<typename Type>
Type c(Type x, const std::vector<Type> &xGrid, std::size_t k){
    std::size_t numOfFinitElems = xGrid.size();
    if (k > numOfFinitElems + 1){
        return NAN;
    }
    Type mult = 1.0;
    for (std::size_t i = 0; i < k; i++){
        mult *= (x - xGrid[i]) / (xGrid[k] - xGrid[i]);
    }
    for (std::size_t i = k + 1; i < numOfFinitElems; i++){
        mult *= (x - xGrid[i]) / (xGrid[k] - xGrid[i]);
    }
    return mult;
}

// Полином Лагранжа
template<typename Type>
Type LagrangePolynom(Type x, const std::vector<Type> &xGrid, const std::vector<Type> &fGrid){
    std::size_t numOfFinitElems = xGrid.size();
    if (numOfFinitElems != fGrid.size()){
        return NAN;
    }
    Type sum = 0.0;
    for (std::size_t k = 0; k < numOfFinitElems; k++){
        sum += c(x, xGrid, k) * fGrid[k]; 
    }
    return sum;
}

// Интерполяция полиномом Лагранжа
template<typename Type>
SOLUTION_FLAG getLagrangeInterpolation(const std::vector<Type> &xVec, std::vector<Type> &fVec, const std::vector<Type> &xGrid, const std::vector<Type> &fGrid){
    std::size_t numOfFinitElems = xGrid.size();
    if (numOfFinitElems != fGrid.size()){
        return NO_SOLUTION;
    }
    fVec.clear();
    std::size_t numOfIntervalVal = xVec.size();
    for (std::size_t i = 0; i < numOfIntervalVal; i++){
        Type sum = 0.0;
        for (std::size_t k = 0; k < numOfFinitElems; k++){
            Type mult = 1.0;
            for (std::size_t j = 0; j < k; j++){
                mult *= (xVec[i] - xGrid[j]) / (xGrid[k] - xGrid[j]);
            }
            for (std::size_t j = k + 1; j < numOfFinitElems; j++){
                mult *= (xVec[i] - xGrid[j]) / (xGrid[k] - xGrid[j]);
            }
            sum += mult * fGrid[k]; 
        }
        fVec.push_back(sum);
    }
    return HAS_SOLUTION;
}

// Сплайн интерполяция

template<typename Type>
SOLUTION_FLAG findSplineCoefs(const std::vector<Type> &xGrid, const std::vector<Type> &fGrid, 
std::vector<Type> &a, std::vector<Type> &b, std::vector<Type> &c, std::vector<Type> &d){
    std::size_t numOfUnits = xGrid.size();
    if (numOfUnits != fGrid.size()){
        return NO_SOLUTION;
    }
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    for (std::size_t i = 1; i < numOfUnits; i++){
        a.push_back(fGrid[i - 1]);
    }
    std::vector<Type> h;
    for (std::size_t i = 1; i < numOfUnits; i++){
        h.push_back(xGrid[i] - xGrid[i - 1]);
    }
    std::vector<Type> g;
    for (std::size_t i = 0; i < numOfUnits - 1; i++){
        g.push_back((fGrid[i + 1] - fGrid[i]) / h[i]);
    }

    // Заполнение матрицы с учетом условия S''(x_0) = S''(x_n) = 0, т.е. c_1 = c_(n+1) = 0 
    std::vector<Type> diag;
    for (std::size_t i = 1; i < numOfUnits - 1; i++){
        diag.push_back(2.0 * (h[i - 1] + h[i]));
    }

    std::vector<Type> lDiag, uDiag;
    for (std::size_t i = 1; i < numOfUnits - 2; i++){
        lDiag.push_back(h[i]);
        uDiag.push_back(h[i]);
    }

    std::vector<Type> rVec;
    for (std::size_t i = 1; i < numOfUnits - 1; i++){
        rVec.push_back(3.0 * (g[i] - g[i - 1]));
    }

    tridiagonalAlgoritm(diag, lDiag, uDiag, rVec, c);
    c.insert(c.begin(), 1, 0.0);

    for (std::size_t i = 0; i < numOfUnits - 1; i++){
        b.push_back(g[i] - (c[i + 1] + 2.0 * c[i]) * h[i] / 3.0);
        d.push_back((c[i + 1] - c[i]) / (3.0 * h[i]));
    }

    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG getCubeSplineInterpolation(const std::vector<Type> &xVec, std::vector<Type> &fVec, const std::vector<Type> &xGrid, const std::vector<Type> &fGrid, 
Type accuracy){
    std::size_t numOfUnits = xGrid.size();
    std::vector<Type> a, b, c, d;
    findSplineCoefs(xGrid, fGrid, a, b, c, d);
    std::size_t i = 0;
    fVec.clear();
    for (std::size_t k = 1; k < numOfUnits; k++){
        while (xGrid[k - 1] <= xVec[i] && xVec[i] <= xGrid[k] + accuracy){
            Type val = xVec[i] - xGrid[k - 1];
            fVec.push_back(a[k - 1] + b[k - 1] * val + c[k - 1] * std::pow(val, 2) + d[k - 1] * std::pow(val, 3));
            i++;
        }
    }
    return HAS_SOLUTION;
}

template<typename Type>
std::size_t getInterpolationErrorsLagrangeUniform(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &uniError, Type accuracy){
    uniError.clear();
    std::vector<Type> xGrid, fGrid;
    for (std::size_t n = 1; n <= numOfErr; n++){
        Type maxErr = 0.0; 
        getUniformGrid(xGrid, fGrid, f, firstX, lastX, n);
        Type tempX = xGrid[0];
        Type h = (xGrid[n] - xGrid[0]) / (2.0 * numOfErr);
        while (tempX <= xGrid[n] + accuracy){
            Type tempErr = std::abs(LagrangePolynom(tempX, xGrid, fGrid) - f(tempX));
            if (tempErr > maxErr){
                maxErr = tempErr;
            }
            tempX += h;
        }
        uniError.push_back(maxErr);
    }
    return uniError.size();
}

template<typename Type>
std::size_t getInterpolationErrorsLagrangeChebyshev(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &chebError, Type accuracy){
    chebError.clear();
    std::vector<Type> xGrid, fGrid;
    for (std::size_t n = 1; n <= numOfErr; n++){
        Type maxErr = 0.0; 
        getChebyshevGrid(xGrid, fGrid, f, firstX, lastX, n);
        Type tempX = xGrid[0];
        Type h = (xGrid[n] - xGrid[0]) / (2.0 * numOfErr);
        while (tempX <= xGrid[n] + accuracy){
            Type tempErr = std::abs(LagrangePolynom(tempX, xGrid, fGrid) - f(tempX));
            if (tempErr > maxErr){
                maxErr = tempErr;
            }
            tempX += h;
        }
        chebError.push_back(maxErr);
    }
    return chebError.size();
}

template<typename Type>
std::size_t getInterpolationErrorsSplineUniform(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &uniError, Type accuracy){
    uniError.clear();
    std::vector<Type> xGrid, fGrid;
    std::vector<Type> a, b, c, d;
    for (std::size_t n = 1; n <= numOfErr; n++){
        Type maxErr = 0.0; 
        getUniformGrid(xGrid, fGrid, f, firstX, lastX, n);
        findSplineCoefs(xGrid, fGrid, a, b, c, d);
        Type tempX = xGrid[0];
        Type h = (xGrid[n] - xGrid[0]) / (2.0 * numOfErr);
        for (std::size_t k = 1; k < n + 1; k++){
            while (xGrid[k - 1] <= tempX && tempX <= xGrid[k] + accuracy){
                Type val = tempX - xGrid[k - 1];
                Type tempErr = std::abs(a[k - 1] + b[k - 1] * val + c[k - 1] * std::pow(val, 2) + d[k - 1] * std::pow(val, 3) - f(tempX));
                if (tempErr > maxErr){
                    maxErr = tempErr;
                }
                tempX += h;
            }
        }
        uniError.push_back(maxErr);
    }
    return uniError.size();
}

template<typename Type>
std::size_t getInterpolationErrorsSplineChebyshev(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfErr,
std::vector<Type> &chebError, Type accuracy){
    chebError.clear();
    std::vector<Type> xGrid, fGrid;
    std::vector<Type> a, b, c, d;
    for (std::size_t n = 1; n <= numOfErr; n++){
        Type maxErr = 0.0; 
        getChebyshevGrid(xGrid, fGrid, f, firstX, lastX, n);
        findSplineCoefs(xGrid, fGrid, a, b, c, d);
        Type tempX = xGrid[0];
        Type h = (xGrid[n] - xGrid[0]) / (2.0 * numOfErr);
        for (std::size_t k = 1; k < n + 1; k++){
            while (xGrid[k - 1] <= tempX && tempX <= xGrid[k] + accuracy){
                Type val = tempX - xGrid[k - 1];
                Type tempErr = std::abs(a[k - 1] + b[k - 1] * val + c[k - 1] * std::pow(val, 2) + d[k - 1] * std::pow(val, 3) - f(tempX));
                if (tempErr > maxErr){
                    maxErr = tempErr;
                }
                tempX += h;
            }
        }
        chebError.push_back(maxErr);
    }
    return chebError.size();
}

template<typename Type>
std::size_t getSpeedEstimateInPoint(Type (*f)(Type x), Type firstX, Type lastX, Type xi, std::size_t numOfFinEl0, 
std::vector<Type> &stepVec, std::vector<Type> &errResult, std::vector<Type> &speedResult, std::size_t stopIt, Type accuracy){
    errResult.clear();
    speedResult.clear();
    std::vector<Type> xGrid, fGrid;
    std::vector<Type> a, b, c, d;
    std::size_t n = numOfFinEl0;
    for (std::size_t i = 0; i < stopIt; i++){
        getUniformGrid(xGrid, fGrid, f, firstX, lastX, n);
        findSplineCoefs(xGrid, fGrid, a, b, c, d);
        stepVec.push_back((xGrid[n] - xGrid[0]) / n);
        //stepVec.push_back(n);
        for (std::size_t k = 1; k < n + 1; k++){
            std::cout << k << '\n';
            if (xGrid[k - 1] <= xi && xi <= xGrid[k] + accuracy){
                Type val = xi - xGrid[k - 1];
                errResult.push_back(
                    std::abs(f(xi) - (a[k - 1] + b[k - 1] * val + c[k - 1] * std::pow(val, 2) + d[k - 1] * std::pow(val, 3)))
                );
                break;
            }
        }
        n = 2 * n;
    }
    for (std::size_t i = 1; i < errResult.size(); i++){
        speedResult.push_back(errResult[i - 1] / errResult[i]);
    }
    return stopIt;
}

template<typename Type>
std::size_t getSpeedEstimate(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfFinEl0, std::vector<Type> &xiVec,
std::vector<Type> &err1, std::vector<Type> &err2, std::vector<Type> &speedResult, Type accuracy){
    speedResult.clear();
    std::vector<Type> xGrid, fGrid;
    std::vector<Type> a, b, c, d;
    std::size_t n = numOfFinEl0;
    getUniformGrid(xGrid, fGrid, f, firstX, lastX, n);
    findSplineCoefs(xGrid, fGrid, a, b, c, d);
    xiVec.clear();
    Type h = (lastX - firstX) / n;
    for (std::size_t i = 0; i < n; i++){
        xiVec.push_back(xGrid[i] + h / 3.0);
        xiVec.push_back(xGrid[i] + 2.0 * h / 3.0);
    }
    std::size_t i = 0;
    err1.clear();
    for (std::size_t k = 1; k < n + 1; k++){
        while (xGrid[k - 1] <= xiVec[i] && xiVec[i] <= xGrid[k] + accuracy){
            Type val = xiVec[i] - xGrid[k - 1];
            err1.push_back(
                std::abs(f(xiVec[i]) - (a[k - 1] + b[k - 1] * val + c[k - 1] * std::pow(val, 2) + d[k - 1] * std::pow(val, 3)))
            );
            i++;
        }
    }
    n *= 2;
    getUniformGrid(xGrid, fGrid, f, firstX, lastX, n);
    findSplineCoefs(xGrid, fGrid, a, b, c, d);
    i = 0;
    err2.clear();
    for (std::size_t k = 1; k < n + 1; k++){
        while (xGrid[k - 1] <= xiVec[i] && xiVec[i] <= xGrid[k] + accuracy){
            Type val = xiVec[i] - xGrid[k - 1];
            err2.push_back(
                std::abs(f(xiVec[i]) - (a[k - 1] + b[k - 1] * val + c[k - 1] * std::pow(val, 2) + d[k - 1] * std::pow(val, 3)))
            );
            i++;
        }
    }
    for (std::size_t i = 0; i < err1.size(); i++){
        speedResult.push_back(err1[i] / err2[i]);
    }
    return xiVec.size();
}

// Лаба 5

template<typename Type>
Type getEquationSolutionBissection(Type (*f)(Type x), Type firstX, Type lastX, Type accuracy, std::size_t stopIteration){
    if (lastX <= firstX){
        return NAN;
    }
    if (std::abs(f(firstX)) < accuracy){
       return firstX;
    }
    if (std::abs(f(lastX)) < accuracy){
       return lastX;
    }
    std::size_t numOfIters = 0;
    Type leftX = firstX;
    Type rightX = lastX;
    Type tempX = (leftX + rightX) / 2.0;
    while ((rightX - leftX) / 2.0 > accuracy){
        if (f(tempX) * f(rightX) > 0.0){
            rightX = tempX;
        }
        if (f(tempX) * f(leftX) > 0.0){
            leftX = tempX;
        }
        tempX = (leftX + rightX) / 2.0;
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return tempX;
}

template<typename Type>
std::size_t getIterationsBissection(Type (*f)(Type x), Type firstX, Type lastX, Type accuracy, std::size_t stopIteration){
    if (lastX <= firstX){
        return 0;
    }
    if (std::abs(f(firstX)) < accuracy){
       return firstX;
    }
    if (std::abs(f(lastX)) < accuracy){
       return lastX;
    }
    std::size_t numOfIters = 0;
    Type leftX = firstX;
    Type rightX = lastX;
    Type tempX = (leftX + rightX) / 2.0;
    while ((rightX - leftX) / 2.0 > accuracy){
        if (f(tempX) * f(rightX) < 0.0){
            leftX = tempX;
        }
        if (f(tempX) * f(leftX) < 0.0){
            rightX = tempX;
        }
        tempX = (leftX + rightX) / 2.0;
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return numOfIters;
}

template<typename Type>
Type diff(Type (*f)(Type x), Type x, Type h){
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

template<typename Type>
Type partialDiff2D(Type (*f)(Type x, Type y), Type x, Type y, VARIABLE_FLAG flag, Type h){
    if (flag == X){
        return (f(x + h, y) - f(x - h, y)) / (2.0 * h);
    }
    if (flag == Y){
        return (f(x, y + h) - f(x, y - h)) / (2.0 * h);
    }
    return NAN;
}

template<typename Type>
void getJacobiMatrix2D(std::vector<std::vector<Type>> &matrix, Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), Type x, Type y, Type h){
    for (std::size_t i = 0; i < matrix.size(); i++){
        matrix[i].clear();
    }
    matrix.clear();
    std::vector<Type> tempRow;
    tempRow.push_back(partialDiff2D(f1, x, y, X, h));
    tempRow.push_back(partialDiff2D(f1, x, y, Y, h));
    matrix.push_back(tempRow);
    tempRow[0] = partialDiff2D(f2, x, y, X, h);
    tempRow[1] = partialDiff2D(f2, x, y, Y, h);
    matrix.push_back(tempRow);
}

template<typename Type>
Type getEquationSolutionNewthon(Type (*f)(Type x), Type firstX, Type h, Type accuracy, std::size_t stopIteration){
    Type tempX = firstX - f(firstX) / diff(f, firstX, h);
    Type prevX = firstX;
    std::size_t numOfIters = 0;
    while (std::abs(tempX - prevX) > accuracy){
        prevX = tempX;
        tempX = prevX - f(prevX) / diff(f, prevX, h);
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return tempX;
}

template<typename Type>
Type getEquationSolutionNewthonModified(Type (*f)(Type x), Type firstX, Type h, Type accuracy, std::size_t stopIteration){
    Type tempX = firstX - f(firstX) / diff(f, firstX, h);
    Type prevX = firstX;
    std::size_t numOfIters = 0;
    Type diff0 = diff(f, firstX, h);
    while (std::abs(tempX - prevX) > accuracy){
        prevX = tempX;
        tempX = prevX - f(prevX) / diff0;
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return tempX;
}

template<typename Type>
std::size_t getIterationsNewthon(Type (*f)(Type x), Type firstX, Type h, Type accuracy, std::size_t stopIteration){
    Type tempX = firstX - f(firstX) / diff(f, firstX, h);
    Type prevX = firstX;
    std::size_t numOfIters = 0;
    while (std::abs(tempX - prevX) > accuracy){
        prevX = tempX;
        tempX = prevX - f(prevX) / diff(f, prevX, h);
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return numOfIters;
}

template<typename Type>
std::size_t getIterationsNewthonModified(Type (*f)(Type x), Type firstX, Type h, Type accuracy, std::size_t stopIteration){
    Type tempX = firstX - f(firstX) / diff(f, firstX, accuracy);
    Type prevX = firstX;
    std::size_t numOfIters = 0;
    Type diff0 = diff(f, firstX, h);
    while (std::abs(tempX - prevX) > accuracy){
        prevX = tempX;
        tempX = prevX - f(prevX) / diff0;
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return numOfIters;
}

template<typename Type>
std::size_t getSystemSolutionNewthon(std::vector<Type> &solution, Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y),
Type firstX, Type firstY, Type h, Type accuracy, Type p, std::size_t stopIteration){
    std::size_t numOfIterations = 1;
    solution.resize(2);
    std::vector<Type> prevSolution = {firstX, firstY};
    std::vector<std::vector<Type>> JacobiMatrix;
    getJacobiMatrix2D(JacobiMatrix , f1, f2, prevSolution[0], prevSolution[1], h);
    std::vector<Type> deltaVec;
    std::vector<Type> rightVec = {-f1(prevSolution[0], prevSolution[1]), -f2(prevSolution[0], prevSolution[1])};
    gaussMethod(JacobiMatrix, rightVec, deltaVec, accuracy);
    for (std::size_t i = 0; i < 2; i++){
        solution[i] = deltaVec[i] + prevSolution[i];
    }
    while(normOfVector(solution - prevSolution, p) > accuracy){
        for (std::size_t i = 0; i < 2; i++){
            prevSolution[i] = solution[i];
        }
        rightVec[0] = -f1(prevSolution[0], prevSolution[1]);
        rightVec[1] = -f2(prevSolution[0], prevSolution[1]);
        getJacobiMatrix2D(JacobiMatrix , f1, f2, prevSolution[0], prevSolution[1], h);
        gaussMethod(JacobiMatrix, rightVec, deltaVec, accuracy);
        for (std::size_t i = 0; i < 2; i++){
            solution[i] = deltaVec[i] + prevSolution[i];
        }
        numOfIterations++;
        if (numOfIterations == stopIteration){
            break;
        }
    }
    return numOfIterations;
}

template<typename Type>
std::size_t getSystemSolutionNewthonAnalytic(std::vector<Type> &solution, std::vector<Type> (*getJacobiMatrixElems)(Type x, Type y), 
Type (*f1)(Type x, Type y), Type (*f2)(Type x, Type y), Type firstX, Type firstY, Type accuracy, Type p, std::size_t stopIteration){
    std::size_t numOfIterations = 1;
    solution.resize(2);
    std::vector<Type> prevSolution = {firstX, firstY};
    std::vector<std::vector<Type>> JacobiMatrix;
    std::vector<Type> JacobiElems = getJacobiMatrixElems(prevSolution[0], prevSolution[1]);
    JacobiMatrix.push_back(std::vector<Type>{JacobiElems[0], JacobiElems[1]});
    JacobiMatrix.push_back(std::vector<Type>{JacobiElems[2], JacobiElems[3]});
    std::vector<Type> deltaVec;
    std::vector<Type> rightVec = {-f1(prevSolution[0], prevSolution[1]), -f2(prevSolution[0], prevSolution[1])};
    gaussMethod(JacobiMatrix, rightVec, deltaVec, accuracy);
    for (std::size_t i = 0; i < 2; i++){
        solution[i] = deltaVec[i] + prevSolution[i];
    }
    while(normOfVector(solution - prevSolution, p) > accuracy){
        for (std::size_t i = 0; i < 2; i++){
            prevSolution[i] = solution[i];
        }
        rightVec[0] = -f1(prevSolution[0], prevSolution[1]);
        rightVec[1] = -f2(prevSolution[0], prevSolution[1]);
        JacobiElems = getJacobiMatrixElems(prevSolution[0], prevSolution[1]);
        JacobiMatrix[0][0] = JacobiElems[0];
        JacobiMatrix[0][1] = JacobiElems[1];
        JacobiMatrix[1][0] = JacobiElems[2];
        JacobiMatrix[1][1] = JacobiElems[3];
        gaussMethod(JacobiMatrix, rightVec, deltaVec, accuracy);
        for (std::size_t i = 0; i < 2; i++){
            solution[i] = deltaVec[i] + prevSolution[i];
        }
        numOfIterations++;
        if (numOfIterations == stopIteration){
            break;
        }
    }
    return numOfIterations;
}

template<typename Type>
FILE_FLAG writeNewthonSwPool(Type (*reF)(Type x, Type y), Type (*imF)(Type x, Type y), Type R, std::size_t n, 
Type h, Type accuracy, const std::string &OUT_FILE_PATH, std::size_t stopIteration){
    std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::vector<Type> rGrid;
    Type step = R / n;
    for (std::size_t i = 1; i < n + 1; i++){
        Type tempR = step * i;
        rGrid.push_back(tempR);
    }
    std::vector<Type> solution;
    for (std::size_t i = 0; i < rGrid.size(); i++){
        Type phi = 0.0;
        while (phi < 2 * M_PI){
            Type x = rGrid[i] * std::cos(phi);
            Type y = rGrid[i] * std::sin(phi);
            getSystemSolutionNewthon(solution, reF, imF, x, y, h, accuracy, 2.0, stopIteration);
            file << x << ' ' << y << ' ' << solution[0] << ' ' << solution[1] << '\n';
            phi += step;
        }
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
std::size_t locoliseRoots(Type (*f)(Type x), Type firstX, Type lastX, std::size_t numOfElems, 
std::vector<std::vector<Type>> &segmentMatrix){
    if (lastX <= firstX){
        return 0;
    }
    for (std::size_t i = 0; i < segmentMatrix.size(); i++){
        segmentMatrix[i].clear();
    }
    segmentMatrix.clear();
    std::vector<Type> xGrid, fGrid;
    getUniformGrid(xGrid, fGrid, f, firstX, lastX, numOfElems);
    std::vector<Type> tempSegment = {0.0, 0.0};
    for (std::size_t i = 0; i < numOfElems; i++){
        if (fGrid[i] * fGrid[i + 1] <= 0){
            tempSegment[0] = xGrid[i];
            tempSegment[1] = xGrid[i + 1];
            segmentMatrix.push_back(tempSegment);
        }
    }
    return segmentMatrix.size();
}               

template<typename Type>
Type getConvergEstimateNewthon(Type (*f)(Type x), Type x_0, Type realX, Type h, Type accuracy, std::size_t stopIteration){
    Type prevX = x_0;
    Type tempX = x_0 - f(x_0) / diff(f, x_0, h);
    Type nextX = tempX - f(tempX) / diff(f, tempX, h);
    Type s = std::log(std::abs((tempX - realX) / (nextX - realX))) / std::log(std::abs((prevX - realX) / (tempX - realX)));
    std::size_t numOfIters = 0;
    while (std::abs(tempX - prevX) > accuracy){
        prevX = tempX;
        tempX = prevX - f(prevX) / diff(f, prevX, h);
        nextX = tempX - f(tempX) / diff(f, tempX, h);
        if (std::abs(nextX - realX) > 2.0 * accuracy){
            s = std::log(std::abs((tempX - realX) / (nextX - realX))) / std::log(std::abs((prevX - realX) / (tempX - realX)));
        }
        if (numOfIters == stopIteration){
            return tempX;
        }
        numOfIters++;
    }
    return s;
} 
