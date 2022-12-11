#ifndef FLAGS_H
#define FLAGS_H

enum SOLUTION_FLAG{
    NO_SOLUTION, // 0
    HAS_SOLUTION // 1 
};

enum FILE_FLAG{
    NOT_OPEN, // 0
    IS_CLOSED // 1
};

enum INVERTIBLE_FLAG{
    NOT_INVERTIBLE, // 0
    IS_INVERTIBLE // 1
};

enum QUADRATIC_FLAG{
    NOT_QUADRATIC, // 0
    IS_QUADRATIC // 1
};

enum MULTIPLIED_FLAG{
    NOT_MULTIPLIED, // 0
    IS_MULTIPLIED // 1
};

enum CORRECT_INPUT_FLAG{
    ERROR, // 0
    CORRECT // 1
};

enum ITERATION_METHOD_FLAG{
    SIMPLE_IT, // 0
    JACOBI, // 1
    RELAXATION // 2
};

enum GRID_FLAG{
    UNIFORM, // 0
    CHEBYSHEV // 1
};

enum VARIABLE_FLAG{
    X, // 0
    Y  // 1  
};

#endif 