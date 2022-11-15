#ifndef _LINALG_H
#define _LINALG_H

#include <cstring>
#include <stdexcept>

template <typename T>
std::string to_string(T var);

template <typename T>
T gcf(T a, T b);

template <typename T>
T lcm(T a, T b);

template <typename T>
struct Vector;

template <typename T>
struct Row {
    int cols;
    T* elems;
    T divisor;

    Row(){}
    Row(const Row& row);
    Row(int cols, T* elems, T divisor);
    Row(int cols, T* elems);
    Row(int cols, T divisor);
    Row(int cols);
    ~Row() {
        free(this->elems);
    }

    void operator= (const Row& row);
    T* operator[] (int index);
    Row operator+ (Row r);
    Row operator- (Row r);
    Row operator* (T n);
    Row operator/ (T n);

    void set(int col, T elem);
};

template <typename T>
struct Matrix {
    int rows;
    int cols;
    T* elems;
    T divisor;

    Matrix(){}
    Matrix(const Matrix& matrix);
    Matrix(int rows, int cols, T* elems, T divisor);
    Matrix(int rows, int cols, T* elems);
    Matrix(int rows, int cols, T divisor);
    Matrix(int rows, int cols);
    Matrix(int side);
    ~Matrix() {
        free(this->elems);
    }

    Row<T> operator[] (int index);
    Matrix operator+ (Matrix m);
    Matrix operator- (Matrix m);
    Matrix operator* (Matrix m);
    Matrix operator* (T n);
    Matrix operator/ (T n);
    bool operator== (Matrix m);

    void setRow(int row, Row<T> r);
    void setElem(int row, int col, T elem);
    Matrix<T> transpose();
    Matrix<T> invert();
};

template <typename T>
struct Vector {
    int dimensions;
    T* elems;
    T divisor;

    Vector(){}
    Vector(const Vector& vector, T divisor);
    Vector(const Vector& vector);
    Vector(int dimensions, T* elems, T divisor);
    Vector(int dimensions, T* elems);
    Vector(int dimensions, T divisor);
    Vector(int dimensions);
    ~Vector() {
        free(this->elems);
    }

    T* operator[] (int index);
    Matrix<T> operator* (Row<T> row);
    Vector operator/ (T n);

    void set(int dimension, T elem);
    Row<T> transpose();
};

template <typename T>
struct Identity : public Matrix<T> {
    Identity(int side);
};

template <typename T>
Row<T>::Row(const Row& row) {
    this->cols = row.cols;
    this->elems = (T*)calloc(cols, sizeof(T));
    for (int i = 0; i < this->cols; i++)
        this->elems[i] = (T)row.elems[i];
    this->divisor = (T)1;
}
template <typename T>
Row<T>::Row(int cols, T* elems, T divisor) {
    this->cols = cols;
    this->elems = (T*)calloc(cols, sizeof(T));
    for (int i = 0; i < cols; i++)
        this->elems[i] = (T)elems[i];
    this->divisor = divisor;
}
template <typename T>
Row<T>::Row(int cols, T* elems) {
    this->cols = cols;
    this->elems = (T*)calloc(cols, sizeof(T));
    for (int i = 0; i < this->cols; i++)
        this->elems[i] = (T)elems[i];
    this->divisor = (T)1;
}

template <typename T>
Row<T>::Row(int cols, T divisor) {
    this->cols = cols;
    this->elems = (T*)calloc(cols, sizeof(T));
    for (int i = 0; i < cols; i++)
        this->elems[i] = T();
    this->divisor = divisor;
}

template <typename T>
Row<T>::Row(int cols) {
    this->cols = cols;
    this->elems = (T*)calloc(cols, sizeof(T));
    for (int i = 0; i < cols; i++)
        this->elems[i] = T();
    this->divisor = (T)1;
}

template <typename T>
Matrix<T>::Matrix(const Matrix& matrix) {
    this->rows = matrix.rows;
    this->cols = matrix.cols;
    this->elems = (T*)calloc(this->rows * this->cols, sizeof(T));
    for (int i = 0; i < this->cols * this->rows; i++)
        this->elems[i] = (T)matrix.elems[i];
    this->divisor = matrix.divisor;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, T* elems, T divisor) {
    this->rows = rows;
    this->cols = cols;
    this->elems = (T*)calloc(rows * cols, sizeof(T));
    for (int i = 0; i < this->cols * this->rows; i++)
        this->elems[i] = (T)elems[i];
    this->divisor = divisor;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, T* elems) {
    this->rows = rows;
    this->cols = cols;
    this->elems = (T*)calloc(rows * cols, sizeof(T));
    for (int i = 0; i < this->cols * this->rows; i++)
        this->elems[i] = (T)elems[i];
    this->divisor = (T)1;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, T divisor) {
    this->rows = rows;
    this->cols = cols;
    this->elems = (T*)calloc(rows * cols, sizeof(T));
    for (int i = 0; i < this->cols * this->rows; i++)
        this->elems[i] = T();
    this->divisor = divisor;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    this->elems = (T*)calloc(rows * cols, sizeof(T));
    for (int i = 0; i < this->cols * this->rows; i++)
        this->elems[i] = T();
    this->divisor = (T)1;
}

template <typename T>
Matrix<T>::Matrix(int side) {
    this->rows = side;
    this->cols = side;
    this->elems = (T*)calloc(side * side, sizeof(T));
    for (int i = 0; i < this->cols * this->rows; i++)
        this->elems[i] = T();
    this->divisor = (T)1;
}

template <typename T>
Vector<T>::Vector(const Vector& vector, T divisor) {
    this->dimensions = vector.dimensions;
    this->elems = (T*)calloc(dimensions, sizeof(T));
    for (int i = 0; i < vector.dimensions; i++)
        this->elems[i] = (T)vector.elems[i];
    this->divisor = divisor;
}

template <typename T>
Vector<T>::Vector(const Vector& vector) {
    this->dimensions = vector.dimensions;
    this->elems = (T*)calloc(dimensions, sizeof(T));
    for (int i = 0; i < vector.dimensions; i++)
        this->elems[i] = (T)vector.elems[i];
    this->divisor = vector.divisor;
}

template <typename T>
Vector<T>::Vector(int dimensions, T* elems, T divisor) {
    this->dimensions = dimensions;
    this->elems = (T*)calloc(dimensions, sizeof(T));
    for (int i = 0; i < dimensions; i++)
        this->elems[i] = (T)elems[i];
    this->divisor = divisor;
}

template <typename T>
Vector<T>::Vector(int dimensions, T* elems) {
    this->dimensions = dimensions;
    this->elems = (T*)calloc(dimensions, sizeof(T));
    memcpy(&this->elems, &elems, this->dimensions * sizeof(T));
    this->divisor = (T)1;
}

template <typename T>
Vector<T>::Vector(int dimensions, T divisor) {
    this->dimensions = dimensions;
    this->elems = (T*)calloc(dimensions, sizeof(T));
    for (int i = 0; i < this->dimensions; i++)
        this->elems[i] = T();
    this->divisor = divisor;
}

template <typename T>
Vector<T>::Vector(int dimensions) {
    this->dimensions = dimensions;
    this->elems = (T*)calloc(dimensions, sizeof(T));
    for (int i = 0; i < this->dimensions; i++)
        this->elems[i] = T();
    this->divisor = (T)1;
}

template <typename T>
Identity<T>::Identity(int side) {
    this->rows = side;
    this->cols = side;
    this->elems = (T*)calloc(side * side, sizeof(T));
    for (int i = 0; i < side; i++)
        for (int j = 0; j < side; j++)
            this->elems[i * side + j] = (T)((i == j) ? 1 : 0);
    this->divisor = T(1);
}

template <typename T>
void Row<T>::operator= (const Row& row){
    this->cols = row.cols;
    this->elems = (T*)calloc(this->cols, sizeof(T));
    for (int i = 0; i < this->cols; i++)
        this->elems[i] = row.elems[i];
    this->divisor = row.divisor;
}

template <typename T>
T* Row<T>::operator[] (int index){
    try {
        if (index < this->cols)
            return this->elems + index;
        else
            throw std::runtime_error("Trying to access column " + to_string(index) + " in matrix with " + to_string(this->cols) + " columns.");
    } catch (const std::exception& e) {
        printf("Index out of range: %s\n", e.what());
    }
    return new T();
}

template <typename T>
Row<T> Row<T>::operator+ (Row r){
    Row sum;
    try {
        if (this->cols == r.cols) {
            sum = Row<T>(this->cols, lcm(this->divisor, r.divisor));
            T left = sum.divisor / this->divisor;
            T right = sum.divisor / r.divisor;
            for (int i = 0; i < sum.cols; i++) {
                *sum.set(i, *(*this)[i] * left + *r[i] * right);
            }
        } else {
            throw std::runtime_error(
                "(" + to_string(this->cols) + " and "
                + to_string(r.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
    return sum;
}

template <typename T>
Row<T> Row<T>::operator- (Row r){
    Row difference;
    try {
        if (this->cols == r.cols) {
            difference = Row<T>(this->cols, lcm(this->divisor, r.divisor));
            T left = difference.divisor / this->divisor;
            T right = difference.divisor / r.divisor;
            for (int i = 0; i < difference.cols; i++) {
                difference.set(i, *(*this)[i] * left - *r[i] * right);
            }
        } else {
            throw std::runtime_error(
                "(" + to_string(this->cols) + " and "
                + to_string(r.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
    return difference;
}

template <typename T>
Row<T> Row<T>::operator* (T n){
    T factor = gcf(this->divisor, n);
    Row product = Row<T>(this->cols, this->divisor / factor);
    factor = n / factor;
    for (int i = 0; i < product.cols; i++)
        product.set(i, *(*this)[i] * factor);
    return product;
}

template <typename T>
Row<T> Row<T>::operator/ (T n){
    Row quotient;
    try {
        if (n == 1)
            return &this;
        int factor = n;
        for (int i = 0; i < this->cols; i++) {
            if (*(*this)[i] != 0)
                factor = gcf(factor, *(*this)[i]);
            if (factor == (T)1)
                return Row<T>(this->cols, this->elems, this->divisor * n);
        }
        Row quotient = Row<T>(this->cols, this->divisor * n / factor);
        for (int i = 0; i < quotient.cols; i++)
            quotient.set(i, *(*this)[i] / factor);
        return quotient;
    } catch (const std::exception& e) {
        printf("Floating point error: %s\n", e.what());
    }
    return quotient;
}

template <typename T>
Row<T> Matrix<T>::operator[] (int index){
    try {
        if (index < this->rows) {
            return Row<T>(this->cols, this->elems + index * this->cols);
        } else {
            throw std::runtime_error("Trying to access row " + to_string(index) + " in matrix with " + to_string(this->rows) + " rows.");
        }
    } catch (const std::exception& e) {
        printf("Index out of range: %s\n", e.what());
    }
    return Row<T>();
}

template <typename T>
Matrix<T> Matrix<T>::operator+ (Matrix m){
    Matrix sum;
    try {
        if (this->rows == m.rows && this->cols == m.cols) {
            sum = Matrix<T>(this->rows, this->cols, lcm(this->divisor, m.divisor));
            T left = sum.divisor / this->divisor;
            T right = sum.divisor / m.divisor;
            for (int i = 0; i < sum.rows; i++)
                for (int j = 0; j < sum.cols; j++)
                    sum.setElem(i, j, *(*this)[i][j] * left + *m[i][j] * right);
        } else {
            throw std::runtime_error(
                "(" + to_string(this->rows) + "x" + to_string(this->cols)+ " and "
                + to_string(m.rows) + "x" + to_string(m.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
    return sum;
}

template <typename T>
Matrix<T> Matrix<T>::operator- (Matrix m){
    Matrix difference;
    try {
        if (this->rows == m.rows && this->cols == m.cols) {
            difference = Matrix<T>(this->rows, this->cols, lcm(this->divisor, m.divisor));
            T left = difference.divisor / this->divisor;
            T right = difference.divisor / m.divisor;
            for (int i = 0; i < difference.rows; i++) {
                for (int j = 0; j < difference.cols; j++) {
                    difference.setElem(i, j, *(*this)[i][j] * left - *m[i][j] * right);
                }
            }
        } else {
            throw std::runtime_error(
                "(" + to_string(this->rows) + "x" + to_string(this->cols)+ " and "
                + to_string(m.rows) + "x" + to_string(m.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
    return difference;
}

template <typename T>
Matrix<T> Matrix<T>::operator* (Matrix m){
    Matrix product;
    try {
        if (this->rows == m.cols && this->cols == m.rows && this->rows >= this->cols) {
            product = Matrix<T>(this->rows);
            for (int i = 0; i < product.rows; i++) {
                for (int j = 0; j < product.cols; j++) {
                    for (int k = 0; k < this->cols; k++) {
                        product.setElem(i, j, *product[i][j] + *(*this)[i][k] * *m[k][j]);
                    }
                }
            }
            if (this->divisor != (T)1 || m.divisor != (T)1)
                product = product / (this->divisor * m.divisor);
        } else {
            throw std::runtime_error(
                "(" + to_string(this->rows) + "x" + to_string(this->cols)+ " and "
                + to_string(m.rows) + "x" + to_string(m.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
    return product;
}

template <typename T>
Matrix<T> Matrix<T>::operator* (T n){
    T factor = gcf(this->divisor, n);
    Matrix product = Matrix<T>(this->rows, this->cols, this->divisor / factor);
    factor = n / factor;
    for (int i = 0; i < product.rows; i++)
        for (int j = 0; j < product.cols; j++)
            product.setElem(i, j, *(*this)[i][j] * factor);
    return product;
}

template <typename T>
Matrix<T> Matrix<T>::operator/ (T n){
    Matrix quotient;
    try {
        if (n == (T)0)
            throw std::runtime_error("Division by 0.");
        else {
            if (n == (T)1)
                return *this;
            T factor = n;
            for (int i = 0; i < this->rows; i++)
                for (int j = 0; j < this->cols; j++) {
                    if (*(*this)[i][j] != (T)0)
                        factor = gcf(factor, *(*this)[i][j]);
                    if (factor == (T)1)
                        return Matrix<T>(this->rows, this->cols, this->elems, this->divisor * n);
                }
            Matrix quotient = Matrix<T>(this->rows, this->cols, this->divisor * n / factor);
            for (int i = 0; i < quotient.rows; i++)
                for (int j = 0; j < quotient.cols; j++)
                    quotient.setElem(i, j, *(*this)[i][j] / factor);
            return quotient;
        }
    } catch (const std::exception& e) {
        printf("Floating point error: %s\n", e.what());
    }
    return quotient;
}

template <typename T>
bool Matrix<T>::operator== (Matrix m){
    if (this->rows == m.rows && this->cols == m.cols) {
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                if (*(*this)[i][j] != *m[i][j])
                    return false;
            }
        }
        return true;
    }
    return false;
}

template <typename T>
T* Vector<T>::operator[] (int index){
    try {
        if (index < this->dimensions) {
            return this->elems + index;
        } else {
            throw std::runtime_error("Trying to access dimension " + to_string(index) + " in vector with " + to_string(this->dimensions) + " dimensions.");
        }
    } catch (const std::exception& e) {
        printf("Index out of range: %s\n", e.what());
    }
    return new T();
}

template <typename T>
Matrix<T> Vector<T>::operator* (Row<T> row) {
    Matrix<T> product;
    try {
        if (this->dimensions == row.cols) {
            product = Matrix<T>(this->dimensions);
            for (int i = 0; i < product.rows; i++)
                for (int j = 0; j < product.cols; j++)
                    product.setElem(i, j, *product[i][j] + *(*this)[i] * *row[j]);
            if (this->divisor != (T)1 || row.divisor != (T)1)
                product = product / (this->divisor * row.divisor);
        } else {
            throw std::runtime_error(
                "(" + to_string(this->dimensions) + " and "
                + to_string(row.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
    return product;
}

template <typename T>
Vector<T> Vector<T>::operator/ (T n) {
    Vector quotient;
    try {
        if (n == (T)0)
            throw std::runtime_error("Division by 0.");
        else {
            if (n == (T)1)
                return *this;
            T factor = n;
            for (int i = 0; i < this->dimensions; i++) {
                if (*(*this)[i] != 0)
                    factor = gcf(factor, *(*this)[i]);
                if (factor == (T)1)
                    return Vector<T>(this->dimensions, this->elems, this->divisor * n);
            }
            Vector quotient = Vector<T>(this->dimensions, this->divisor * n / factor);
            for (int i = 0; i < quotient.dimensions; i++)
                *quotient.set(i, *(*this)[i] / factor);
            return quotient;
        }
    } catch (const std::exception& e) {
        printf("Floating point error: %s\n", e.what());
    }
    return quotient;
    
}

template <typename T>
void print(Row<T> row) {
    using namespace std;
    printf("{");
    for (int i = 0; i < row.cols; i++) {
        if (i > 0)
            printf(",\t");
        printf("%s", to_string(*row[i]).c_str());
    }
    printf("\t} / %s\n", to_string(row.divisor).c_str());
}

template <typename T>
void print(Matrix<T> matrix) {
    using namespace std;
    printf("{");
    for (int i = 0; i < matrix.rows; i++) {
        if (i > 0)
            printf("\n ");
        printf("{");
        for (int j = 0; j < matrix.cols; j++) {
            if (j > 0)
                printf(",\t");
            printf("%s", to_string(*matrix[i][j]).c_str());
        }
        printf("\t}");
    }
    printf("} / %s\n", to_string(matrix.divisor).c_str());
}

template <typename T>
void print(Vector<T> vector) {
    using namespace std;
    printf("{");
    for (int i = 0; i < vector.dimensions; i++) {
        if (i > 0)
            printf("\n ");
        printf("{%s\t}", to_string(*vector[i]).c_str());
    }
    printf("} / %s\n", to_string(vector.divisor).c_str());
}

template <typename T>
Row<T> Vector<T>::transpose() {
    return Row<T>(this->dimensions, this->elems, this->divisor);
}

template <typename T>
void Row<T>::set(int col, T elem) {
    try {
        if (col >= 0 && col < this->cols)
            this->elems[col] = (T)elem;
        else
            throw std::runtime_error("Trying to set column " + to_string(col) + " in row with " + to_string(this->cols) + " columns.");
    } catch (const std::exception& e) {
        printf("Index out of range: %s\n", e.what());
    }
}

template <typename T>
void Matrix<T>::setRow(int row, Row<T> r) {
    try {
        if (this->cols == r.cols)
            for (int i = 0; i < this->cols; i++)
                this->elems[row * this->cols + i] = (T)r.elems[i];
        else {
            throw std::runtime_error(
                "(" + to_string(this->cols) + " and "
                + to_string(r.cols) + ")."
            );
        }
    } catch (const std::exception& e) {
        printf("Incompatible dimensions %s\n", e.what());
    }
}

template <typename T>
void Matrix<T>::setElem(int row, int col, T elem) {
    try {
        if (row < 0 || row >= this->rows)
            throw std::runtime_error("Trying to set row " + to_string(row) + " in matrix with " + to_string(this->rows) + " rows.");
        else if (col < 0 || col >= this->cols)
            throw std::runtime_error("Trying to set column " + to_string(col) + " in row with " + to_string(this->cols) + " columns.");
        else
            this->elems[row * this->cols + col] = (T)elem;
    } catch (const std::exception& e) {
        printf("Index out of range: %s\n", e.what());
    }
}

template <typename T>
void Vector<T>::set(int dimension, T elem) {
    try {
        if (dimension >= 0 && dimension < this->dimensions)
            this->elems[dimension] = (T)elem;
        else
            throw std::runtime_error("Trying to set dimension " + to_string(dimension) + " in vector with " + to_string(this->dimensions) + " dimensions.");
    } catch (const std::exception& e) {
        printf("Index out of range: %s\n", e.what());
    }
}

template <typename T>
Matrix<T> Matrix<T>::transpose() {
    Matrix<T> transposed = Matrix<T>(this->cols, this->rows);
    for (int i = 0; i < this->cols; i++) {
        for (int j = 0; j < this->rows; j++) {
            transposed.setElem(i, j, *(*this)[j][i]);
        }
    }
    return transposed;
}

template <typename T>
Matrix<T> Matrix<T>::invert() {
    try {
        if (this->rows == this->cols) {
            Matrix<T> initial = Matrix<T>(this->rows, this->cols, this->elems);
            Matrix<T> inverted = Identity<T>(this->rows) * this->divisor;
            for (int i = 0; i < this->rows; i++) {
                if (*initial[i][i] == 0) {
                    for (int j = i + 1; j < this->rows; j++) {
                        if (*initial[j][i] != 0) {
                            Row<T> temp = Row<T>(this->cols);
                            temp = initial[i];
                            initial.setRow(i, initial[j]);
                            initial.setRow(j, temp);
                            temp = inverted[i];
                            inverted.setRow(i, inverted[j]);
                            inverted.setRow(j, temp);
                            break;
                        }
                    }
                }
                for (int j = 0; j < this->rows; j++) {
                    if (i == j || *initial[j][i] == 0)
                        continue;
                    T factor = gcf(*initial[i][i], *initial[j][i]);
                    T left = *initial[i][i] / factor;
                    T right = *initial[j][i] / factor;
                    initial.setRow(j, initial[j] * left - initial[i] * right);
                    inverted.setRow(j, inverted[j] * left - inverted[i] * right);
                }
            }
            T scale = 1;
            for (int i = 0; i < this->rows; i++) {
                scale = lcm(scale, *initial[i][i]);
            }
            if (initial == Matrix<T>(this->rows)) {
                throw std::runtime_error("Matrix is singular.");
            } else {
                for (int i = 0; i < this->rows; i++) {
                    inverted.setRow(i, inverted[i] * (scale / *initial[i][i]));
                }
                return inverted / scale;
            }
        } else {
            throw std::runtime_error("Non-square matrix.");
        }
    } catch (const std::exception& e) {
        printf("Matrix is not invertible: %s\n", e.what());
    }
    return Matrix<T>();
}

template <typename T>
std::string to_string(T var) {
    return to_string(var);
}

template <typename T>
T gcf(T a, T b) {
    if (a == (T)0 || b == (T)0)
        return T();
    while (a > T(1) && b > T(1)) {
        if (a % b == T())
            return b;
        a = a % b;
        if (b % a == T())
            return a;
        b = b % a;
    }
    return 1;
}

template <typename T>
T lcm(T a, T b) {
    return a * b / gcf(a, b);
}

#endif