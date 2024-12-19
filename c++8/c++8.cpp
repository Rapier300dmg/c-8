#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstdlib>
#include <ctime>

template <typename T>
class Matrix {
private:
    T** data;
    size_t rows;
    size_t cols;

    void allocateMemory() {
        data = new T * [rows];
        for (size_t i = 0; i < rows; ++i) {
            data[i] = new T[cols]();
        }
    }

    void freeMemory() {
        for (size_t i = 0; i < rows; ++i) {
            delete[] data[i];
        }
        delete[] data;
    }

public:
    Matrix(size_t r, size_t c) : rows(r), cols(c) {
        allocateMemory();
    }

    ~Matrix() {
        freeMemory();
    }

    void input() {
        std::cout << "Enter matrix elements (" << rows << "x" << cols << "):\n";
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cin >> data[i][j];
            }
        }
    }

    void fillRandom(T minValue, T maxValue) {
        std::srand(std::time(nullptr));
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                data[i][j] = minValue + (std::rand() % (maxValue - minValue + 1));
            }
        }
    }

    void display() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << std::setw(5) << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    T maxElement() const {
        T maxVal = data[0][0];
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (data[i][j] > maxVal) {
                    maxVal = data[i][j];
                }
            }
        }
        return maxVal;
    }

    T minElement() const {
        T minVal = data[0][0];
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (data[i][j] < minVal) {
                    minVal = data[i][j];
                }
            }
        }
        return minVal;
    }

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrix dimensions do not allow multiplication");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix operator/(T scalar) const {
        if (scalar == 0) {
            throw std::invalid_argument("Division by zero");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] / scalar;
            }
        }
        return result;
    }
};

int main() {
    try {
        Matrix<int> mat1(3, 3);
        mat1.fillRandom(1, 10);
        std::cout << "Matrix 1:" << std::endl;
        mat1.display();

        Matrix<int> mat2(3, 3);
        mat2.fillRandom(5, 15);
        std::cout << "Matrix 2:" << std::endl;
        mat2.display();

        Matrix<int> sum = mat1 + mat2;
        std::cout << "Sum of matrices:" << std::endl;
        sum.display();

        std::cout << "Max element in Matrix 1: " << mat1.maxElement() << std::endl;
        std::cout << "Min element in Matrix 2: " << mat2.minElement() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
