#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>
#include <bit>
#include <random>
#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include <iomanip>
#include <Windows.h>
#include <algorithm>
#include <set>

using namespace boost::multiprecision;
using namespace std;

typedef vector<vector<int>> Matrix;
typedef vector<vector<double>> DoubleMatrix;

//diff
template<typename T>
void print(vector<vector<T>> A) {
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++) {
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
}
string create_matrix_file(int size, int amount, int lside, int rside);
string create_matrix_file(int size, int amount, double lside, double rside, const int prec);
string create_matrix_01(int size, int amount, bool sparce);
string create_once_stoch_matrix(int size, int amount, int prec);
string create_twice_stoch_matrix(int size, int amount, int prec);
string create_circ_matrix(int size, int amount, int lside, int rside);
string create_circ_matrix(int size, int amount, double lside, double rside, const int prec);
string create_perm_matrix(int size, int amount);
string create_triang_matrix(int size, int amount, int lside, int rside);

vector<Matrix> get_from_file(string filename);
vector<DoubleMatrix> get_from_file_double(string filename);
void wait_for_enter();
Matrix double_to_int(const DoubleMatrix& A, int& fix);

enum Method { RYSER = 1, LAPLAS, GLYNN, GG, AB, MC };
enum Menu { GET_DATA_FROM_FILE = 1, CALCULATE_DATA, VIEW_DATA, GENERATE };
enum Menu_Calculate {ALL = 1, PREC, APPROX, CHOOSE};
enum Menu_Generate {INT_A_TO_B = 1, FLOAT_A_TO_B, SPARCE_DENSE, ONCE_STOCH, TWICE_STOCH, INT_CIRCULANT, FLOAT_CIRCULANT, PERMATURE, TRIANGULAR};

//prec
/* 1*/ cpp_int ryser_perm(const Matrix& A);
/* 2*/ cpp_int laplas_perm(const Matrix& A);
/* 3*/ cpp_int glynn_perm(const Matrix& A);

//approx
/* 4*/ cpp_int godsil_gotman_perm(const Matrix& A, int N = 5000);
/* 5*/ cpp_int AB_perm(const Matrix& A, int N = 500);
/* 6*/ double mc_perm(const DoubleMatrix& A);

//det
cpp_int det(const Matrix& A);
double detDouble(const DoubleMatrix& A);

//the calculation
void calculation(vector<Matrix>, vector<Method> methods = {RYSER, LAPLAS, GLYNN, GG, AB});
void calculation(vector<DoubleMatrix>, vector<Method> methods = {RYSER, LAPLAS, GLYNN, MC, AB});

//Fact
const int factN = 40; 
vector<cpp_int> fact(factN + 1, 0);
void build_factorial();

//rand
static std::mt19937 gen(std::chrono::steady_clock::now().time_since_epoch().count());


int main()
{
    build_factorial();

    string matr = "None!";
    size_t n;
    vector<Matrix> user_matrix;
    vector<DoubleMatrix> user_matrix_d;


    string temp_matr;
    string v;

    bool check_for_double = false;
    int check = 0;

    int size;
    int amount;
    int lb;
    int rb;

    double lbdouble;
    double rbdouble;
    int prec;
    do {
        cout << "Permanent of Matrix" << endl;
        cout << "Current matrix datafile: " << matr << endl;
        cout << "1 - set matrix data from file" << endl;
        cout << "2 - calculation of permanent with multiple methods with time" << endl;
        cout << "3 - view current matrix" << endl;
        cout << "4 - generate a set of matrixes" << endl;

        cout << "0 - exit" << endl;
        cout << endl;

        cin >> check;

        try {
            switch (check) {
            case GET_DATA_FROM_FILE:
                cout << "Does Matrix data contains double? Y/N" << endl;
                cin >> v;
                if (v == "y" || v == "Y") {
                    check_for_double = true;
                }
                
                cout << "Input the filename or full path to .txt" << endl;
                cin >> matr;

                if (matr.length() < 4 || matr.substr(matr.length() - 4) != ".txt") {
                    matr += ".txt";
                    cout << "Added .txt extension: " << matr << endl;
                }
                if (check_for_double == true) {
                    user_matrix_d = get_from_file_double(matr);
                }
                else {
                    user_matrix = get_from_file(matr);
                }
                wait_for_enter();
                break;
            case CALCULATE_DATA:
                if (matr != "None!") {


                    cout << "Ryser" << endl;
                    cout << "Laplas" << endl;
                    cout << "Glynn" << endl;
                    cout << "Godsil Gotman (approx, only (0,1)-matrixes)" << endl;
                    cout << "Peter MC's (approx, only twice stohastic)" << endl;
                    cout << "AB (approx, only non negative)" << endl;
                    cout << endl;

                    cout << "Choose one of following options: " << endl;
                    cout << "1 - use all methods" << endl;
                    cout << "2 - use only precisious methods" << endl;
                    cout << "3 - use only approx methods" << endl;
                    cout << "4 - choose methods" << endl;
                    cin >> check;
                    switch (check) {
                    case ALL:
                        if (check_for_double == false) {
                            calculation(user_matrix);
                        }
                        else {
                            calculation(user_matrix_d);
                        }
                        break;
                    case PREC:
                        if (check_for_double == false) {
                            calculation(user_matrix, { RYSER, LAPLAS, GLYNN });
                        }
                        else {
                            calculation(user_matrix_d, {RYSER, LAPLAS, GLYNN});
                        }
                        break;
                    case APPROX:
                        if (check_for_double == false) {
                            calculation(user_matrix, { GG, AB });
                        }
                        else {
                            calculation(user_matrix_d, { AB, MC });
                        }
                        break;
                    case CHOOSE:
                        set<Method> user_methods;
                        do {
                            cout << "1 - Ryser" << endl;
                            cout << "2 - Laplas" << endl;
                            cout << "3 - Glynn" << endl;
                            cout << "4 - Godsil Gotman (approx, (0,1) only)" << endl;
                            cout << "5 - MC (approx, only twice stohastic)" << endl;
                            cout << "6 - AB (approx, only non negative)" << endl;
                            cout << "0 - done" << endl;
                            cin >> check;

                            switch (check) {
                            case RYSER:
                                user_methods.insert(RYSER);
                                cout << "Ryser method was added." << endl;
                                break;
                            case LAPLAS:
                                user_methods.insert(LAPLAS);
                                cout << "Laplas method was added." << endl;
                                break;
                            case GLYNN:
                                user_methods.insert(GLYNN);
                                cout << "Glynn method was added." << endl;
                                break;
                            case GG:
                                user_methods.insert(GG);
                                cout << "GG method was added." << endl;
                                break;
                            case MC:
                                user_methods.insert(MC);
                                cout << "MC method was added." << endl;
                                break;
                            case AB:
                                user_methods.insert(AB);
                                cout << "AB method was added." << endl;
                                break;
                            default:
                                cout << "No methods were chosen!" << endl;
                                break;
                            }
                            wait_for_enter();

                            system("cls");

                        } while (check != 0);
                        vector<Method> methods(user_methods.begin(), user_methods.end());
                        break;
                    }
                    check = 2;
                    wait_for_enter();
                }
                else {
                    cout << "Data is empty, please input a file." << endl;
                    wait_for_enter();
                }
                break;
            case VIEW_DATA:
                cout << "Current data: " << matr << endl;
                if (matr == "None!") {
                    cout << "There are no matrixes!" << endl;
                }
                else {
                    if (check_for_double) {
                        vector<DoubleMatrix> tempD = get_from_file_double(matr);
                        for (size_t i = 0; i < tempD.size(); i++) {
                            print(tempD[i]);
                            cout << endl;
                        }
                    }
                    else {
                        vector<Matrix> temp = get_from_file(matr);
                        for (size_t i = 0; i < temp.size(); i++) {
                            print(temp[i]);
                            cout << endl;
                        }
                    }
                }
                wait_for_enter();
                break;
            case GENERATE:
                do {
                    cout << "Input size of matrixes: " << endl;
                    cin >> size;
                    if (size < 1) {
                        cout << "It has to be bigger than 0!" << endl;
                    }
                } while (size < 1);

                do {
                    cout << "Input amount of matrixes: " << endl;
                    cin >> amount;
                    if (amount < 1) {
                        cout << "It has to be bigger than 0" << endl;
                    }
                } while (amount < 1);

                cout << "Generate a set of matrixes" << endl;
                cout << "1 - Matrixes with integer elements from [a;b]" << endl;
                cout << "2 - Matrixes with float elements from [a;b]" << endl;
                cout << "3 - Sparce and dense (0,1)-matrixes" << endl;
                cout << "4 - Once stohastic matrixes" << endl;
                cout << "5 - Twice stohastic matrixes" << endl;
                cout << "6 - Integer circulant matrixes from [a;b]" << endl;
                cout << "7 - Float circulant matrixes from [a;b]" << endl;
                cout << "8 - Permature matrixes" << endl;
                cout << "9 - Triangular matrixes from [a;b]" << endl;
                cout << "All distributions are uniform" << endl;
                cin >> check;
                switch (check) {
                case INT_A_TO_B:
                    cout << "Int [a;b] matrixes" << endl;
                   
                    cout << "Input left and right border of matrixes' elements" << endl;
                    cout << "Left border: " << endl;
                    cin >> lb;
                    do {
                        cout << "Right border: " << endl;
                        cin >> rb;
                        if (rb < lb) {
                            cout << "Right border has to be bigger than left!" << endl;
                        }
                    } while (rb < lb);

                    check_for_double = false;
                    temp_matr = create_matrix_file(size, amount, lb, rb);
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                case FLOAT_A_TO_B:
                    cout << "Float [a;b] matrixes" << endl;

                    cout << "Input left and right border of matrixes' elements" << endl;
                    cout << "Left border: " << endl;
                    cin >> lbdouble;
                    do {
                        cout << "Right border: " << endl;
                        cin >> rbdouble;
                        if (rbdouble < lbdouble) {
                            cout << "Right border has to be bigger than left!" << endl;
                        }
                    } while (rbdouble < lbdouble);

                    do {
                        cout << "Input amount of digits after point" << endl;
                        cin >> prec;
                        if (prec < 1) {
                            cout << "It has to be bigger than 1!" << endl;
                        }
                    } while (prec < 1);

                    check_for_double = true;
                    temp_matr = create_matrix_file(size, amount, lbdouble, rbdouble, prec);
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                case SPARCE_DENSE:
                    cout << "(0,1) Sparce or Dense matrixes" << endl;
                    cout << "Do you want for matrix to be Sparce (~75% - 0) Y/N?" << endl;
                    cout << "Otherwise it will be Dense" << endl;
                    cin >> v;
                    if (v == "y" || v == "Y") {
                        temp_matr = create_matrix_01(size, amount, true);
                    }
                    else {
                        temp_matr = create_matrix_01(size, amount, false);
                    }
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                case ONCE_STOCH:
                    cout << "Once stochastic matrixes" << endl;
                    cout << "Sum of elements might differ by +-1 * 10^(prec) from 1" << endl;
                    do {
                        cout << "Input amount of digits after point" << endl;
                        cin >> prec;
                        if (prec < 1) {
                            cout << "It has to be bigger than 1!" << endl;
                        }
                    } while (prec < 1);

                    check_for_double = true;
                    temp_matr = create_once_stoch_matrix(size, amount, prec);
                    break;
                case TWICE_STOCH:
                    cout << "Twice stochastic matrixes" << endl;
                    cout << "Sum of elements might differ by +-0.001 from 1" << endl;
                    do {
                        cout << "Input amount of digits after point" << endl;
                        cin >> prec;
                        if (prec < 1) {
                            cout << "It has to be bigger than 1!" << endl;
                        }
                    } while (prec < 1);

                    check_for_double = true;
                    temp_matr = create_twice_stoch_matrix(size, amount, prec);
                    break;
                case INT_CIRCULANT:
                    cout << "Int [a;b] circulant matrixes" << endl;

                    cout << "Input left and right border of matrixes' elements" << endl;
                    cout << "Left border: " << endl;
                    cin >> lb;
                    do {
                        cout << "Right border: " << endl;
                        cin >> rb;
                        if (rb < lb) {
                            cout << "Right border has to be bigger than left!" << endl;
                        }
                    } while (rb < lb);

                    check_for_double = false;
                    temp_matr = create_circ_matrix(size, amount, lb, rb);
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                case FLOAT_CIRCULANT:
                    cout << "Float [a;b] circulant matrixes" << endl;

                    cout << "Input left and right border of matrixes' elements" << endl;
                    cout << "Left border: " << endl;
                    cin >> lbdouble;
                    do {
                        cout << "Right border: " << endl;
                        cin >> rbdouble;
                        if (rbdouble < lbdouble) {
                            cout << "Right border has to be bigger than left!" << endl;
                        }
                    } while (rbdouble < lbdouble);

                    do {
                        cout << "Input amount of digits after point" << endl;
                        cin >> prec;
                        if (prec < 1) {
                            cout << "It has to be bigger than 1!" << endl;
                        }
                    } while (prec < 1);

                    check_for_double = true;
                    temp_matr = create_circ_matrix(size, amount, lbdouble, rbdouble, prec);
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                case PERMATURE:
                    temp_matr = create_perm_matrix(size, amount);
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                case TRIANGULAR:
                    cout << "Int [a;b] triangular matrixes" << endl;

                    cout << "Input left and right border of matrixes' elements" << endl;
                    cout << "Left border: " << endl;
                    cin >> lb;
                    do {
                        cout << "Right border: " << endl;
                        cin >> rb;
                        if (rb < lb) {
                            cout << "Right border has to be bigger than left!" << endl;
                        }
                    } while (rb < lb);

                    check_for_double = false;
                    temp_matr = create_triang_matrix(size, amount, lb, rb);
                    cout << "Created matrixes in: " << temp_matr << endl;
                    break;
                }


                cout << "Do you want to set current matrix data to matrixes from created file? Y/N" << endl;
                cin >> v;
                if (v == "y" || v == "Y") {
                    matr = temp_matr;
                    if(check_for_double){
                        user_matrix_d = get_from_file_double(matr);
                    }
                    else {
                        user_matrix = get_from_file(matr);
                    }
                }

                wait_for_enter();
                break;
            default:
                cout << "No options where chosen" << endl;
                wait_for_enter();
            }
            system("cls");
        }
        catch (exception e) {
            cout << "Exception!" << endl;
            cout << e.what() << endl;
            wait_for_enter();
            system("cls");
        }
    } while (check != 0);
}

string create_matrix_file(int n, int k, int lside, int rside) {
    string filename = "random_matrixes_size" + to_string(n) + "_amount" + to_string(k) + "l" + to_string(lside) + "r" + to_string(rside) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);
    int minVal = min(lside, rside);
    int maxVal = max(lside, rside);
    std::uniform_int_distribution<> dist(minVal, maxVal);


    srand(time(0));
    for (int v = 0; v < k; v++) {
        out << n << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out << dist(gen) << ' ';
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_matrix_file(int n, int k, double lside, double rside, const int prec) {
    string filename = "random_matrixes_size" + to_string(n) + "_amount" + to_string(k) + "l" + to_string(lside) + "r" + to_string(rside) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);
    double minVal = min(lside, rside);
    double maxVal = max(lside, rside);
    std::uniform_real_distribution<double> dist(minVal, maxVal);

    for (int v = 0; v < k; v++) {
        out << n << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out << std::fixed << std::setprecision(prec) << dist(gen) << ' ';
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_matrix_01(int n, int k, bool sparce) {
    string filename;
    if (sparce) {
        filename = "random_sparce_matrixes_size" + to_string(n) + "Amount" + to_string(k) + ".txt";
    }
    else {
        filename = "random_dense_matrixes_size" + to_string(n) + "Amount" + to_string(k) + ".txt";
    }
    ofstream out(filename, ios::out | ios::trunc);

    std::uniform_int_distribution<int> dist(0, 3);

    for (int v = 0; v < k; v++) {
        out << n << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist(gen) == 0) {
                    if (sparce) {
                        out << 1 << ' ';
                    }
                    else {
                        out << 0 << ' ';
                    }
                }
                else {
                    if (sparce) {
                        out << 1 << ' ';
                    }
                    else {
                        out << 0 << ' ';
                    }
                }
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_once_stoch_matrix(int n, int k, int prec) {
    string filename;
    filename = "random_once_stoch_matrixes_size" + to_string(n) + "_amount" + to_string(k) + ".txt";
    
    ofstream out(filename, ios::out | ios::trunc);

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int v = 0; v < k; v++) {
        out << n << endl;
        for (int i = 0; i < n; i++) {
            vector<double> row(n);
            double sum = 0.0;

            for (int j = 0; j < n; j++) {
                row[j] = dist(gen);
                sum += row[j];
            }

            for (int j = 0; j < n; j++) {
                out << fixed << setprecision(prec) << row[j] / sum << ' ';
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_twice_stoch_matrix(int n, int k, int prec) {
    string filename = "random_twice_stoch_matrixes_size" + to_string(n) + "_amount" + to_string(k) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);

    uniform_real_distribution<double> dist(0.0, 1.0);

    for (int v = 0; v < k; v++) {
        DoubleMatrix matrix(n, vector<double>(n));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = dist(gen);
            }
        }

        for (int i = 0; i < n; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < n; j++) {
                row_sum += matrix[i][j];
            }
            for (int j = 0; j < n; j++) {
                matrix[i][j] /= row_sum;
            }
        }

        for (int iter = 0; iter < 100; iter++) {
            for (int j = 0; j < n; j++) {
                double col_sum = 0.0;
                for (int i = 0; i < n; i++) {
                    col_sum += matrix[i][j];
                }
                for (int i = 0; i < n; i++) {
                    matrix[i][j] /= col_sum;
                }
            }

            for (int i = 0; i < n; i++) {
                double row_sum = 0.0;
                for (int j = 0; j < n; j++) {
                    row_sum += matrix[i][j];
                }
                for (int j = 0; j < n; j++) {
                    matrix[i][j] /= row_sum;
                }
            }
        }

        out << n << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out << fixed << setprecision(prec) << matrix[i][j] << ' ';
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_circ_matrix(int size, int amount, int lside, int rside) {
    string filename = "random_circ_matrixes_size" + to_string(size) + "_amount" + to_string(amount) + "l" + to_string(lside) + "r" + to_string(rside) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);

    int minVal = min(lside, rside);
    int maxVal = max(lside, rside);
    std::uniform_int_distribution<int> dist(minVal, maxVal);

    for (int v = 0; v < amount; v++) {
        vector<int> first_row(size);
        for (int i = 0; i < size; i++) {
            first_row[i] = dist(gen);
        }

        out << size << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                out << first_row[(j - i + size) % size] << ' ';
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_circ_matrix(int size, int amount, double lside, double rside, const int prec) {
    string filename = "random_circ_matrixes_size" + to_string(size) + "_amount" + to_string(amount) + "l" + to_string(lside) + "r" + to_string(rside) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);

    double minVal = min(lside, rside);
    double maxVal = max(lside, rside);
    std::uniform_real_distribution<double> dist(minVal, maxVal);

    for (int v = 0; v < amount; v++) {
        vector<double> first_row(size);
        for (int i = 0; i < size; i++) {
            first_row[i] = dist(gen);
        }

        out << size << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                out << fixed << setprecision(prec) << first_row[(j - i + size) % size] << ' ';
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_perm_matrix(int size, int amount) {
    string filename = "random_perm_matrixes_size" + to_string(size) + "_amount" + to_string(amount) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);

    for (int v = 0; v < amount; v++) {
        vector<int> first(size);
        for (int i = 0; i < size; i++) {
            first[i] = i;
        }
        shuffle(first.begin(), first.end(), gen);

        out << size << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (first[i] == j) {
                    out << 1 << ' ';
                }
                else {
                    out << 0 << ' ';
                }
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}
string create_triang_matrix(int size, int amount, int lside, int rside) {
    string filename = "random_triang_matrixes_size" + to_string(size) + "_amount" + to_string(amount) + "l" + to_string(lside) + "r" + to_string(rside) + ".txt";
    ofstream out(filename, ios::out | ios::trunc);

    int minVal = min(lside, rside);
    int maxVal = max(lside, rside);
    std::uniform_int_distribution<int> dist(minVal, maxVal);

    for (int v = 0; v < amount; v++) {
        vector<int> first(size);
        for (int i = 0; i < size; i++) {
            first[i] = dist(gen);
        }

        out << size << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    out << first[i] << ' ';
                }
                else {
                    out << 0 << ' ';
                }
            }
            out << endl;
        }
    }
    cout << "Done." << endl;
    return filename;
}

vector<Matrix> get_from_file(string filename) {
    ifstream in(filename);
    if (!in.is_open()) {
        throw runtime_error("Error: Cannot open file");
    }
    string l;
    vector<Matrix> get;

    size_t n;
    Matrix A;
    while (in >> l) {
        n = stoi(l);
        A = vector<vector<int>>(n, vector<int>(n, 0));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if(!(in >> A[i][j])){
                    throw runtime_error("Error: Failed to read file");
                }
            }
        }
        get.push_back(A);
    }
    return get;
}
vector<DoubleMatrix> get_from_file_double(string filename) {
    ifstream in(filename);
    if (!in.is_open()) {
        throw runtime_error("Error: Cannot open file ");
    }
    string l;
    vector<DoubleMatrix> get;

    size_t n;
    DoubleMatrix A;
    while (in >> l) {
        n = stoi(l);
        A = vector<vector<double>>(n, vector<double>(n, 0));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                in >> A[i][j];
            }
        }
        get.push_back(A);
    }
    return get;
}
void wait_for_enter() {
    cout << "\nPress Enter to continue...";
    cin.ignore(10000, '\n');
    cin.get();
}
Matrix double_to_int(const DoubleMatrix& A, int& fix) {
    Matrix B(A.size(), vector<int>(A[0].size(), 0));
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++) {
            double temp = A[i][j];
            int count = 0;
            double eps = 1e-12;

            while (abs(temp - round(temp)) > eps && count < 12) {
                temp *= 10;
                count++;
            }
            if(count > fix){
                fix = count;
            }
            B[i][j] = temp;
        }
    }
    return B;
}

cpp_int ryser_perm(const Matrix& A) {
    size_t n = A.size();
    cpp_int perm = 0;
    unsigned long long number = (1ULL << n);

    int count;
    int sign;

    for (unsigned long long S = 1; S < number; S++) {
        cpp_int product = 1;

        for (size_t i = 0; i < n; i++) {
            cpp_int row_sum = 0;
            for (size_t j = 0; j < n; j++) {
                if ((S >> j) & 1) {
                    row_sum += A[i][j];
                }
            }
            product *= row_sum;
        }

        count = popcount(S);
        sign = ((n - count) % 2 == 1) ? -1 : 1;
        perm += ((n - count) % 2 == 1) ? -product : product;

    }
    return perm;
}
cpp_int laplas_perm(const Matrix& A) {
    size_t n = A.size();

    if (n == 1) {
        return cpp_int(A[0][0]);
    }

    cpp_int perm = 0;

    Matrix B;
    for (unsigned long long j = 0; j < n; j++) {
        if (A[0][j] == 0) continue;
        B = vector<vector<int>>(n - 1, vector<int>(n - 1, 0));

        for (unsigned long long i = 0; i < n - 1; i++) {
            for (unsigned long long k = 0; k < n - 1; k++) {
                if (k < j) {
                    B[i][k] = A[i + 1][k];
                }
                else {
                    B[i][k] = A[i + 1][k + 1];
                }
            }
        }
        perm += cpp_int(A[0][j]) * laplas_perm(B);
    }

    return perm;
}
cpp_int glynn_perm(const Matrix& A) {
    size_t n = A.size();
    unsigned long long number = 1ULL << (n - 1);
    cpp_int perm = 0;

    for (unsigned long long B = 0; B < number; B++) {
        int sign = (popcount(B) & 1) ? -1 : 1;

        cpp_int product = 1;
        for (size_t i = 0; i < n; i++) {
            cpp_int sum = A[i][0];
            unsigned long long mask = B;
            for (size_t j = 1; j < n; j++) {
                if (mask & 1) sum -= A[i][j];
                else sum += A[i][j];
                mask >>= 1;
            }
            product *= sum;
        }
        perm += sign * product;
    }
    return perm / number;
}

cpp_int godsil_gotman_perm(const Matrix& A, int N) {
    size_t n = A.size();
    cpp_int perm = 0;
    std::uniform_int_distribution<> dis(0, 1);

    Matrix B;
    cpp_int detB;
    for (int k = 0; k < N; k++) {
        B = A;
        for (unsigned long long i = 0; i < n; i++) {
            for (unsigned long long j = 0; j < n; j++) {
                if (dis(gen)) {
                    B[i][j] = -B[i][j];
                }
            }
        }
        detB = det(B);
        perm += (detB * detB);
    }
    return perm/N;
}
cpp_int AB_perm(const Matrix& A, int N) {
    size_t n = A.size();
    normal_distribution<> dist(0.0, 1.0);

    cpp_int perm = 0;
    double det;
    double g;
    DoubleMatrix B;

    for (int k = 0; k < N; k++) {
        B = vector<vector<double>>(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                g = dist(gen);
                B[i][j] = sqrt((double)A[i][j]) * g;
            }
        }
        det = detDouble(B);
        perm += cpp_int(det * det);
    }

    return perm / N;
}
double mc_perm(const DoubleMatrix& A) {
    size_t n = A.size();

    DoubleMatrix AtA(n, vector<double>(n, 0.0));
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                AtA[i][j] += A[k][i] * A[k][j];
            }
        }
    } 

    //считаем I + J + A'A
    double inv_n = 1.0 / n;
    DoubleMatrix M(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M[i][j] = (i == j ? 1.0 : 0.0) + inv_n - AtA[i][j];
        }
    }

    double detM = detDouble(M);
    //для таких матриц в следствии округления определитель может быть очень маленьким
    detM = fabs(detM);
    if (detM < 1e-12) detM = 1e-12;

    double log_n_fact = 0.0;
    for (int i = 2; i <= n; ++i) log_n_fact += log(static_cast<double>(i));

    double perm_nA = exp(log_n_fact - 0.5 * log(detM));
    double perm = perm_nA / pow(n, n);

    return perm;
}

cpp_int det(const Matrix& A) {
    Matrix M = A;
    size_t n = M.size();
    cpp_int det = 1;
    int sign = 1;

    size_t pivot;
   

    for (size_t i = 0; i < n; i++) {
        long long factor;
        pivot = i;
        while (pivot < n && M[pivot][i] == 0) {
            pivot++;
        }

        if (pivot == n) {
            return 0;
        }

        if (pivot != i) {
            swap(M[i], M[pivot]);
            sign *= -1;
        }

        for (size_t j = i + 1; j < n; j++) {
            while (M[j][i] != 0) {
                if (abs(M[j][i]) < abs(M[i][i])) {
                    swap(M[j], M[i]);
                    sign *= -1;
                }
                factor = M[j][i] / M[i][i];
                for (size_t k = i; k < n; k++) {
                    M[j][k] -= factor * M[i][k];
                }
            }
        }
        det *= M[i][i];
    }
    return det * sign;
}
double detDouble(const DoubleMatrix& A) {
    size_t n = A.size();

    DoubleMatrix M = A;
    double det = 1.0;
    double sign = 1.0;
    double eps = 1e-12;

    size_t pivot;
    double max_val;

    for (size_t i = 0; i < n; i++) {
        pivot = i;
        max_val = fabs(M[i][i]);
        for (size_t k = i + 1; k < n; k++) {
            if (fabs(M[k][i]) > max_val) {
                max_val = fabs(M[k][i]);
                pivot = k;
            }
        }

        if (max_val < eps) {
            return 0.0;
        }

        if (pivot != i) {
            swap(M[i], M[pivot]);
            sign *= -1.0;
        }

        for (size_t j = i + 1; j < n; j++) {
            double factor = M[j][i] / M[i][i];
            for (size_t k = i; k < n; k++) {
                M[j][k] -= factor * M[i][k];
            }
        }

        det *= M[i][i];
    }

    return det * sign;
}

void calculation(vector<Matrix> matrixes, vector<Method> methods) {

    //vector<Matrix> matrixes = get_from_file(filename);
    bool has_laplas = false;
    bool laplas_check = false;

    for (size_t i = 0; i < methods.size() && !has_laplas; i++) {
        if (methods[i] == LAPLAS) {
            has_laplas = true;
        }
    }
    if (has_laplas) {
        for (size_t i = 0; i < matrixes.size() && !laplas_check; i++) {
            if (matrixes[i].size() > 11) {
                laplas_check = true;
            }
        }
        if (laplas_check) {
            string t;
            cout << "Methods include Laplas' method of calculating permanent." << endl;
            cout << "Data has matrix of size > 11" << endl;
            cout << "Calculating it with Laplas' method WILL take a lot of time" << endl;
            cout << "It's not recommended to test it." << endl;
            cout << "Do you want to keep the method? Y/N" << endl;
            cin >> t;
            has_laplas = (t == "y" || t == "Y") ? true : false;
        }
    }
    if (!has_laplas) {
        methods.erase(remove(methods.begin(), methods.end(), LAPLAS), methods.end());
    }

    vector<vector<cpp_int>> results(matrixes.size(), vector<cpp_int>(methods.size()));
    vector<vector<double>> times(matrixes.size(), vector<double>(methods.size()));

    const int NAME_WIDTH = 15;
    const int TIME_WIDTH = 12;
    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << "Processing matrix " << j + 1 << "/" << matrixes.size() << "...\r" << flush;
        for (size_t m = 0; m < methods.size(); m++) {
            chrono::high_resolution_clock::time_point start;

            switch (methods[m]) {

            case RYSER:
                start = chrono::high_resolution_clock::now();
                results[j][m] = ryser_perm(matrixes[j]);
                break;

            case LAPLAS:
                if (has_laplas) {
                    start = chrono::high_resolution_clock::now();
                    results[j][m] = laplas_perm(matrixes[j]);
                }
                break;

            case GLYNN:
                start = chrono::high_resolution_clock::now();
                results[j][m] = glynn_perm(matrixes[j]);
                break;

            case GG:
                for (size_t i = 0; i < matrixes[j].size(); i++) {
                    for (size_t v = 0; v < matrixes[j].size(); v++) {
                        if (matrixes[j][i][v] < 0) {
                            throw std::runtime_error("GG_perm: Matrix contains negative elements");
                        }
                    }
                }
                start = chrono::high_resolution_clock::now();
                results[j][m] = godsil_gotman_perm(matrixes[j]);
                break;

            case AB:
                for (size_t i = 0; i < matrixes[j].size(); i++) {
                    for (size_t v = 0; v < matrixes[j].size(); v++) {
                        if (matrixes[j][i][v] < 0) {
                            throw std::runtime_error("AB_perm: Matrix contains negative elements");
                        }
                    }
                }
                start = chrono::high_resolution_clock::now();
                results[j][m] = AB_perm(matrixes[j]);
                break;
            }
            auto end = chrono::high_resolution_clock::now();
            times[j][m] = chrono::duration<double>(end - start).count();
            //cout << results[j][m];
        }
    }

    cout << "\nTIME (seconds)\n";
    cout << left << setw(NAME_WIDTH) << "Matrix";
    for (auto& method : methods) {
        switch (method) {
        case RYSER: cout << setw(TIME_WIDTH) << "Ryser"; break;
        case LAPLAS: cout << setw(TIME_WIDTH) << "Laplas"; break;
        case GLYNN: cout << setw(TIME_WIDTH) << "Glynn"; break;
        case GG: cout << setw(TIME_WIDTH) << "GG"; break;
        case AB: cout << setw(TIME_WIDTH) << "AB"; break;
        }
    }
    cout << "\n";
    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    vector<double> time_sums(methods.size(), 0.0);
    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << left << setw(NAME_WIDTH) << ("Matrix " + to_string(j + 1));
        for (size_t m = 0; m < methods.size(); m++) {
            cout << setw(TIME_WIDTH) << fixed << setprecision(4) << times[j][m];
            time_sums[m] += times[j][m];
        }
        cout << endl;
    }

    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    cout << left << setw(NAME_WIDTH) << "AVERAGE";
    for (size_t m = 0; m < methods.size(); m++) {
        double avg = time_sums[m] / matrixes.size();
        cout << setw(TIME_WIDTH) << fixed << setprecision(4) << avg;
    }
    cout << endl;

    cout << "\nPERMANENT VALUES\n";
    cout << left << setw(NAME_WIDTH) << "Matrix";

    for (auto& method : methods) {
        switch (method) {
        case RYSER: cout << setw(TIME_WIDTH) << "Ryser"; break;
        case LAPLAS: cout << setw(TIME_WIDTH) << "Laplas"; break;
        case GLYNN: cout << setw(TIME_WIDTH) << "Glynn"; break;
        case GG: cout << setw(TIME_WIDTH) << "GG"; break;
        case AB: cout << setw(TIME_WIDTH) << "AB"; break;
        }
    }
    cout << "\n";

    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << left << setw(NAME_WIDTH) << ("Matrix " + to_string(j + 1));

        for (size_t m = 0; m < methods.size(); m++) {
            string val = results[j][m].str();
            if (val.length() > TIME_WIDTH) {
                val = val.substr(0, TIME_WIDTH - 3) + "...";
            }
            cout << setw(TIME_WIDTH) << val;
        }
        cout << endl;
    }
}
void calculation(vector<DoubleMatrix> matrixes, vector<Method> methods) {

    //vector<Matrix> matrixes = get_from_file(filename);
    bool has_laplas = false;
    bool laplas_check = false;

    for (size_t i = 0; i < methods.size() && !has_laplas; i++) {
        if (methods[i] == LAPLAS) {
            has_laplas = true;
        }
    }
    if (has_laplas) {
        for (size_t i = 0; i < matrixes.size() && !laplas_check; i++) {
            if (matrixes[i].size() > 11) {
                laplas_check = true;
            }
        }
        if (laplas_check) {
            string t;
            cout << "Methods include Laplas' method of calculating permanent." << endl;
            cout << "Data has matrix of size > 11" << endl;
            cout << "Calculating it with Laplas' method WILL take a lot of time" << endl;
            cout << "It's not recommended to test it." << endl;
            cout << "Do you want to keep the method? Y/N" << endl;
            cin >> t;
            has_laplas = (t == "y" || t == "Y") ? true : false;
        }
    }
    if (!has_laplas) {
        methods.erase(remove(methods.begin(), methods.end(), LAPLAS), methods.end());
    }

    vector<vector<double>> results(matrixes.size(), vector<double>(methods.size()));
    vector<vector<double>> times(matrixes.size(), vector<double>(methods.size()));

    const int NAME_WIDTH = 15;
    const int TIME_WIDTH = 12;
    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << "Processing matrix " << j + 1 << "/" << matrixes.size() << "...\r" << flush;
        int fix;
        Matrix B = double_to_int(matrixes[j], fix);
        cpp_int power = pow(cpp_int(10), fix);
        cpp_rational result = 0;
        for (size_t m = 0; m < methods.size(); m++) {
            chrono::high_resolution_clock::time_point start;

            switch (methods[m]) {

            case RYSER:
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(ryser_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;

            case LAPLAS:
                if (has_laplas) {
                    start = chrono::high_resolution_clock::now();
                    result = cpp_rational(laplas_perm(B)) / power;
                    results[j][m] = result.convert_to<double>();
                }
                break;

            case GLYNN:
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(glynn_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;

            case GG:
                for (size_t i = 0; i < matrixes[j].size(); i++) {
                    for (size_t v = 0; v < matrixes[j].size(); v++) {
                        if (matrixes[j][i][v] < 0) {
                            throw std::runtime_error("GG_perm: Matrix contains negative elements");
                        }
                        if (matrixes[j][i][v] != 0 && matrixes[j][i][v] != 1) {
                            throw std::runtime_error("GG_perm: Matrix has to contain 0 or 1s only");
                        }
                    }
                }
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(godsil_gotman_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;

            case AB:
                for (size_t i = 0; i < matrixes[j].size(); i++) {
                    for (size_t v = 0; v < matrixes[j].size(); v++) {
                        if (matrixes[j][i][v] < 0) {
                            throw std::runtime_error("AB_perm: Matrix contains negative elements");
                        }
                    }
                }
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(AB_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;
            }
            auto end = chrono::high_resolution_clock::now();
            times[j][m] = chrono::duration<double>(end - start).count();
        }
    }

    cout << "\nTIME (seconds)\n";
    cout << left << setw(NAME_WIDTH) << "Matrix";
    for (auto& method : methods) {
        switch (method) {
        case RYSER: cout << setw(TIME_WIDTH) << "Ryser"; break;
        case LAPLAS: cout << setw(TIME_WIDTH) << "Laplas"; break;
        case GLYNN: cout << setw(TIME_WIDTH) << "Glynn"; break;
        case GG: cout << setw(TIME_WIDTH) << "GG"; break;
        case AB: cout << setw(TIME_WIDTH) << "AB"; break;
        }
    }
    cout << "\n";
    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    vector<double> time_sums(methods.size(), 0.0);
    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << left << setw(NAME_WIDTH) << ("Matrix " + to_string(j + 1));
        for (size_t m = 0; m < methods.size(); m++) {
            cout << setw(TIME_WIDTH) << fixed << setprecision(4) << times[j][m];
            time_sums[m] += times[j][m];
        }
        cout << endl;
    }

    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    cout << left << setw(NAME_WIDTH) << "AVERAGE";
    for (size_t m = 0; m < methods.size(); m++) {
        double avg = time_sums[m] / matrixes.size();
        cout << setw(TIME_WIDTH) << fixed << setprecision(4) << avg;
    }
    cout << endl;

    cout << "\nPERMANENT VALUES\n";
    cout << left << setw(NAME_WIDTH) << "Matrix";

    for (auto& method : methods) {
        switch (method) {
        case RYSER: cout << setw(TIME_WIDTH) << "Ryser"; break;
        case LAPLAS: cout << setw(TIME_WIDTH) << "Laplas"; break;
        case GLYNN: cout << setw(TIME_WIDTH) << "Glynn"; break;
        case GG: cout << setw(TIME_WIDTH) << "GG"; break;
        case AB: cout << setw(TIME_WIDTH) << "AB"; break;
        }
    }
    cout << "\n";

    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << left << setw(NAME_WIDTH) << ("Matrix " + to_string(j + 1));

        for (size_t m = 0; m < methods.size(); m++) {
            string val = to_string(results[j][m]);
            if (val.length() > TIME_WIDTH) {
                val = val.substr(0, TIME_WIDTH - 3) + "...";
            }
            cout << setw(TIME_WIDTH) << val;
        }
        cout << endl;
    }
}

void build_factorial() {
    fact[0] = 1;
    for (int i = 1; i <= factN; ++i) {
        fact[i] = fact[i - 1] * i;
    }
}