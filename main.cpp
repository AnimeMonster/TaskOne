#include <iostream>
#include <iomanip>
#include <cassert>

class Matrix {
public:
    int rows, columns;
    double **value;

    Matrix(int row, int column, bool diag = false) {
        columns = column;
        rows = row;
        value = new double *[row];
        for (int s = 0; s < row; s++) {
            value[s] = new double[column];
            for (int k = 0; k < column; k++) {
                value[s][k] = 0;
                if (diag && s == k)
                    value[s][k] = 1;
            }
        }
    }

    ~Matrix() {
        for (int i = 0; i < rows; i++) {
            delete value[i];
        }
        delete value;
    }

    Matrix(const Matrix &x) {
        rows = x.rows;
        columns = x.columns;
        value = new double *[rows];
        for (int s = 0; s < rows; s++) {
            value[s] = new double[columns];
            for (int k = 0; k < columns; k++) {
                value[s][k] = x.value[s][k];
            }
        }
    }

    void rearr_rows(int a, int b) {
        for (int k = 0; k < columns; k++) {
            double v = value[a][k];
            value[a][k] = value[b][k];
            value[b][k] = v;
        }
    }

    void triangle(Matrix &affil) {
        for (int ro = 0; ro < rows; ++ro) {
            for (int dimin = ro + 1; dimin < rows; ++dimin) {
                int change = ro;
                for (int d = ro; d < rows; d++) {
                    if (value[d][ro] > value[change][ro])
                        change = d;
                }
                rearr_rows(ro, change);
                affil.rearr_rows(ro, change);
                double k = value[dimin][ro] / value[ro][ro];
                for (int i = ro; i < columns; ++i) {
                    value[dimin][i] -= k * value[ro][i];
                }
                for (int i = 0; i < affil.columns; ++i) {
                    affil.value[dimin][i] -= k * affil.value[ro][i];
                }
            }
        }
    }

    void diagonal(Matrix &affil) {
        triangle(affil);
        for (int ro = rows - 1; ro >= 0; --ro) {
            for (int dimin = 0; dimin < ro; ++dimin) {
                double k = value[dimin][ro] / value[ro][ro];
                value[dimin][ro] -= k * value[ro][ro];
                for (int i = 0; i < affil.columns; ++i) {
                    affil.value[dimin][i] -= k * affil.value[ro][i];
                }
            }
        }
    }

    Matrix multiply(Matrix &mat) {
        Matrix comp(rows, mat.columns);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < mat.columns; ++j) {
                double comp_ij = 0;
                for (int k = 0; k < columns; ++k) {
                    comp_ij += value[i][k] * mat.value[k][j];
                }
                comp.value[i][j] = comp_ij;
            }
        }
        return comp;
    }

    Matrix inverse() {
        assert(rows == columns);
        Matrix inv(rows, columns, true);
        Matrix t(*this);
        t.diagonal(inv);
        for (int i = 0; i < rows; i++) {
            if (t.value[i][i] == 0) {
                throw std::runtime_error("Virozhdenay matrica");
            }
        }
        for (int s = 0; s < rows; s++) {
            for (int k = 0; k < columns; k++) {
                inv.value[s][k] /= t.value[s][s];
            }
        }
        return inv;
    }

    void print() {
        for (int s = 0; s < rows; s++) {
            std::cout << "| ";
            for (int k = 0; k < columns; k++) {
                std::cout << value[s][k] << "\t";
            }
            std::cout << "|" << std::endl;
        }
        std::cout << std::endl;
    }

};

int main() {
    std::cout << std::setprecision(2) << std::fixed;
    int rows, columns;
    std::cout << "Vvedite razmernost' matrici a:" << std::endl;
    std::cin >> rows;
    std::cin >> columns;
    Matrix a(rows, columns);
    std::cout << "VVedite znacheniy matrici a: " << std::endl;
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.columns; j++) {
            double d;
            std::cout << i << " " << j << ":";
            std::cin >> d;
            a.value[i][j] = d;
        }
    }
    std::cout << "Vvedite znacheniy matrici b: " << std::endl;
    Matrix b(4, 1);
    for (int i = 0; i < b.rows; i++) {
        for (int j = 0; j < b.columns; j++) {
            double f;
            std::cout << i << " " << j << ":";
            std::cin >> f;
            b.value[i][j] = f;
        }
    }
    std::cout << "Obratnay: " << std::endl;
    Matrix a_inv = a.inverse();
    a_inv.print();

    Matrix x = a_inv.multiply(b);
    std::cout << "x: " << std::endl;
    x.print();
    return 0;
}
