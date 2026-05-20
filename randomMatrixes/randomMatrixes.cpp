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
#include <cmath>

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
bool has_negative(const Matrix& A);
bool has_negative(const DoubleMatrix& A);
bool is_01(const Matrix& A);
bool is_01(const DoubleMatrix& A, double prec = 1e-12);

vector<Matrix> get_from_file(string filename);
vector<DoubleMatrix> get_from_file_double(string filename);
void wait_for_enter();
Matrix double_to_int(const DoubleMatrix& A, int& fix);

enum Method { RYSER = 1, LAPLAS, GLYNN, EJ, GG, AB, MC };
enum Menu { GET_DATA_FROM_FILE = 1, CALCULATE_DATA, VIEW_DATA, GENERATE };
enum Menu_Calculate {ALL = 1, PREC, APPROX, CHOOSE};
enum Menu_Generate {INT_A_TO_B = 1, FLOAT_A_TO_B, SPARCE_DENSE, ONCE_STOCH, TWICE_STOCH, INT_CIRCULANT, FLOAT_CIRCULANT, PERMATURE, TRIANGULAR};

//prec
/* 1*/ cpp_int ryser_perm(const Matrix& A);
/* 2*/ cpp_int laplas_perm(const Matrix& A);
/* 3*/ cpp_int glynn_perm(const Matrix& A);
/* 4*/ cpp_int G(int j, vector<cpp_int> row_sums, vector<int> f, Matrix B);
       cpp_int ej_perm(const Matrix& A);

//approx
/* 5*/ cpp_int godsil_gotman_perm(const Matrix& A, int N = 400);
/* 6*/ cpp_int AB_perm(const Matrix& A, int N = 400);
/* 7*/ double mc_perm(const DoubleMatrix& A);

//det
cpp_int det(const Matrix& A);
double detDouble(const DoubleMatrix& A);


//the calculation
void calculation(vector<Matrix>, vector<Method> methods = {RYSER, LAPLAS, GLYNN, GG, AB, EJ});
void calculation(vector<DoubleMatrix>, vector<Method> methods = {RYSER, LAPLAS, GLYNN, MC, AB, EJ, GG});

void the_table(const vector<Method>& methods, size_t size, const vector<vector<cpp_int>>& results, const vector<vector<double>>& times);
void the_table(const vector<Method>& methods, size_t size, const vector<vector<double>>& results, const vector<vector<double>>& times);


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
                    cout << "Erix B, Joel Frankin ((0,1) only)" << endl;
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
                            calculation(user_matrix, { RYSER, LAPLAS, GLYNN, EJ });
                        }
                        else {
                            calculation(user_matrix_d, { RYSER, LAPLAS, GLYNN, EJ });
                        }
                        break;
                    case APPROX:
                        if (check_for_double == false) {
                            calculation(user_matrix, { GG, AB });
                        }
                        else {
                            calculation(user_matrix_d, { GG, AB, MC });
                        }
                        break;
                    case CHOOSE:
                        set<Method> user_methods;
                        do {
                            cout << "1 - Ryser - prec, any matrixes" << endl;
                            cout << "2 - Laplas - prec, any matrixes with n < 11" << endl;
                            cout << "3 - Glynn - prec, any matrixes" << endl;
                            cout << "4 - E. Bax, J. Frankin ((0,1) only)" << endl;
                            cout << "5 - Godsil Gutman (approx, (0,1) only)" << endl;
                            cout << "6 - Alexander Barvinok (approx, only non negative)" << endl;
                            cout << "7 - Peter McCullagh (approx, only twice stohastic)" << endl;
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
                            case EJ:
                                user_methods.insert(EJ);
                                cout << "E. Bax, J. Frankin's method was added." << endl;
                                break;
                            case GG:
                                user_methods.insert(GG);
                                cout << "Godsil Gutman's method was added." << endl;
                                break;
                            case MC:
                                if (check_for_double) {
                                    user_methods.insert(MC);
                                    cout << "Peter McCullagh's method was added." << endl;
                                }
                                else {
                                    cout << "Cannot choose that method, matrix has to be with float elements" << endl;
                                }
                                break;
                            case AB:
                                user_methods.insert(AB);
                                cout << "Alexander Barvinok's method was added." << endl;
                                break;
                            case 0:
                                break;
                            default:
                                cout << "No methods were chosen!" << endl;
                                break;
                            }
                            wait_for_enter();
                            system("cls");

                        } while (check != 0);
                        vector<Method> methods(user_methods.begin(), user_methods.end());
                        if (check_for_double && !empty(methods)) {
                            calculation(user_matrix_d, methods);
                        }
                        else if(!empty(methods)){
                            calculation(user_matrix, methods);
                        }
                        else {
                            cout << "No options were chosen." << endl;
                        }
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
                        out << 0 << ' ';
                    }
                    else {
                        out << 1 << ' ';
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
cpp_int G(int j, vector<cpp_int> row_sums, vector<int> f, Matrix B) {
    size_t n = B.size();
    for (int i = 0; i < n; ++i) {
        if (f[i] == j - 1 && row_sums[i] == 0)
            return 0;
    }
    if (j == n) {
        cpp_int prod = 1;
        for (int i = 0; i < n; ++i) {
            prod *= row_sums[i];
        }
        return prod;
    }
    int u = (j % 2 == 0) ? 1 : -1;
    for (int i = 0; i < n; ++i) {
        if (B[i][j] == 1) {
            row_sums[i] += u;
        }
    }
    cpp_int a = G(j + 1, row_sums, f, B);
    for (int i = 0; i < n; ++i) {
        if (B[i][j] == 1) {
            row_sums[i] -= u;
        }
    }
    cpp_int b = G(j + 1, row_sums, f, B);
    return (a - b) * u;
}
cpp_int ej_perm(const Matrix& A) {
    int n = A.size();
    if (n == 0) return 0;

    Matrix B = A;
    vector<int> col_order;
    vector<int> avail;
    for (int j = 0; j < n; ++j) {
        avail.push_back(j);
    }

    while (!avail.empty()) {
        int best = -1;
        int min = n + 1;
        for (int i = 0; i < n; ++i) {
            int count = 0;

            for (int c = 0; c < avail.size(); c++) {
                if (A[i][avail[c]] == 1) {
                    count++;
                }
            }
            if (count > 0 && count < min) {
                min = count;
                best = i;
            }
        }
        if (best == -1) {
            col_order.insert(col_order.end(), avail.begin(), avail.end());
            break;
        }

        vector<int> ones_cols;
        for (int c = 0; c < avail.size(); c++) {
            if (A[best][avail[c]] == 1) {
                ones_cols.push_back(avail[c]);
            }
        }
        col_order.insert(col_order.end(), ones_cols.begin(), ones_cols.end());

        vector<int> new_avail;
        for (int c = 0; c < avail.size(); c++) {
            if (find(ones_cols.begin(), ones_cols.end(), avail[c]) == ones_cols.end())
                new_avail.push_back(avail[c]);
        }
        avail.swap(new_avail);
    }
    for (int i = 0; i < n; ++i) {
        vector<int> new_row(n);
        for (int j = 0; j < n; ++j) new_row[j] = A[i][col_order[j]];
        B[i] = new_row;
    }

    vector<int> f(n, -1);
    for (int i = 0; i < n; ++i) {
        for (int j = n - 1; j >= 0; --j) {
            if (B[i][j] == 1) {
                f[i] = j;
                break;
            }
        }
    }
    vector<cpp_int> row_sums(n, 0);
    return G(0, row_sums, f, B);
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
    size_t n = A.size();

    vector<vector<cpp_int>> M(n, vector<cpp_int>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            M[i][j] = cpp_int(A[i][j]);

    cpp_int sign = 1;
    cpp_int prev = 1;

    for (size_t k = 0; k < n - 1; ++k) {
        size_t pivot_row = k;
        while (pivot_row < n && M[pivot_row][k] == 0) ++pivot_row;
        if (pivot_row == n) return cpp_int(0);

        if (pivot_row != k) {
            swap(M[k], M[pivot_row]);
            sign = -sign;
        }

        const cpp_int pivot = M[k][k];
        for (size_t i = k + 1; i < n; ++i) {
            for (size_t j = k + 1; j < n; ++j) {
                M[i][j] = (M[i][j] * pivot - M[i][k] * M[k][j]) / prev;
            }
            M[i][k] = 0;
        }
        prev = pivot;
    }
    return sign * M[n - 1][n - 1];
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
    bool has_laplas = false;
    bool has_ab = false;
    bool has_gg = false;
    bool has_ej = false;
    bool laplas_check = false;

    bool has_negative_elements = false;
    bool has_non_01 = false;

    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i] == LAPLAS) {
            has_laplas = true;
        }
        if (methods[i] == AB) {
            has_ab = true;
        }
        if (methods[i] == EJ) {
            has_ej = true;
        }
        if (methods[i] == GG) {
            has_gg = true;
        }
    }
    for (size_t i = 0; i < matrixes.size(); i++) {
        if (matrixes[i].size() > 11) {
            laplas_check = true;
        }
        if (has_negative(matrixes[i])) {
            has_negative_elements = true;
        }
        if (!is_01(matrixes[i])) {
            has_non_01 = true;
        }
    }

    if (has_negative_elements && has_ab) {
        cout << "Methods include AB' method of calculating permanent." << endl;
        cout << "Data has matrix with negative elements" << endl;
        cout << "Method will be removed." << endl;
        methods.erase(remove(methods.begin(), methods.end(), AB), methods.end());
    }
    if (has_non_01 && has_gg) {
        cout << "Methods include GG' method of calculating permanent." << endl;
        cout << "Data has matrix that isn't (0,1)" << endl;
        cout << "Method will be removed." << endl;
        methods.erase(remove(methods.begin(), methods.end(), GG), methods.end());
    }
    if (has_non_01 && has_ej) {
        cout << "Methods include EJ' method of calculating permanent." << endl;
        cout << "Data has matrix that isn't (0,1)" << endl;
        cout << "Method will be removed." << endl;
        methods.erase(remove(methods.begin(), methods.end(), EJ), methods.end());
    }

    if (has_laplas) {
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

    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << "Processing matrix " << j + 1 << "/" << matrixes.size() << "...\r" << flush;
        size_t n = matrixes[j].size();
        n = n * n;
        for (size_t m = 0; m < methods.size(); m++) {
            chrono::high_resolution_clock::time_point start;

            switch (methods[m]) {

            case RYSER:
                start = chrono::high_resolution_clock::now();
                results[j][m] = ryser_perm(matrixes[j]);
                break;

            case LAPLAS:
                start = chrono::high_resolution_clock::now();
                results[j][m] = laplas_perm(matrixes[j]);
                
                break;

            case GLYNN:
                start = chrono::high_resolution_clock::now();
                results[j][m] = glynn_perm(matrixes[j]);
                break;

            case GG:
                start = chrono::high_resolution_clock::now();
                results[j][m] = godsil_gotman_perm(matrixes[j], n);
                break;

            case AB:
                start = chrono::high_resolution_clock::now();
                results[j][m] = AB_perm(matrixes[j], n);
                break;
            case EJ:
                start = chrono::high_resolution_clock::now();
                results[j][m] = ej_perm(matrixes[j]);
                break;
            }

            auto end = chrono::high_resolution_clock::now();
            times[j][m] = chrono::duration<double>(end - start).count();
        }
    }

    the_table(methods, matrixes.size(), results, times);
}
void calculation(vector<DoubleMatrix> matrixes, vector<Method> methods) {
    bool has_laplas = false;
    bool has_ab = false;
    bool has_gg = false;
    bool has_ej = false;
    bool laplas_check = false;

    bool has_negative_elements = false;
    bool has_non_01 = false;

    for (size_t i = 0; i < methods.size(); i++) {
        if (methods[i] == LAPLAS) {
            has_laplas = true;
        }
        if (methods[i] == AB) {
            has_ab = true;
        }
        if (methods[i] == EJ) {
            has_ej = true;
        }
        if (methods[i] == GG) {
            has_gg = true;
        }
    }
    for (size_t i = 0; i < matrixes.size(); i++) {
        if (matrixes[i].size() > 11) {
            laplas_check = true;
        }
        if (has_negative(matrixes[i])) {
            has_negative_elements = true;
        }
        if (!is_01(matrixes[i])) {
            has_non_01 = true;
        }
    }

    if (has_negative_elements && has_ab) {
        cout << "Methods include AB' method of calculating permanent." << endl;
        cout << "Data has matrix with negative elements" << endl;
        cout << "Method will be removed." << endl;
        methods.erase(remove(methods.begin(), methods.end(), AB), methods.end());
    }
    if (has_non_01 && has_gg) {
        cout << "Methods include GG' method of calculating permanent." << endl;
        cout << "Data has matrix that isn't (0,1)" << endl;
        cout << "Method will be removed." << endl;
        methods.erase(remove(methods.begin(), methods.end(), GG), methods.end());
    }
    if (has_non_01 && has_ej) {
        cout << "Methods include EJ' method of calculating permanent." << endl;
        cout << "Data has matrix that isn't (0,1)" << endl;
        cout << "Method will be removed." << endl;
        methods.erase(remove(methods.begin(), methods.end(), EJ), methods.end());
    }

    if (has_laplas) {
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

    for (size_t j = 0; j < matrixes.size(); j++) {
        cout << "Processing matrix " << j + 1 << "/" << matrixes.size() << "...\r" << flush;
        int fix = 0;
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
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(laplas_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
               
                break;

            case GLYNN:
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(glynn_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;

            case GG:
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(godsil_gotman_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;

            case AB:
                start = chrono::high_resolution_clock::now();
                result = cpp_rational(AB_perm(B)) / power;
                results[j][m] = result.convert_to<double>();
                break;
            case MC:
                start = chrono::high_resolution_clock::now();
                results[j][m] = mc_perm(matrixes[j]);
                break;
            }
            auto end = chrono::high_resolution_clock::now();
            times[j][m] = chrono::duration<double>(end - start).count();
        }
    }
    the_table(methods, matrixes.size(), results, times);
}

void the_table(const vector<Method>& methods, size_t size, const vector<vector<cpp_int>>& results, const vector<vector<double>>& times) {
    const int NAME_WIDTH = 25;
    const int TIME_WIDTH = 25;

    cout << "\nTIME (seconds)\n";
    cout << left << setw(NAME_WIDTH) << "Matrix";
    for (auto& method : methods) {
        switch (method) {
        case RYSER: cout << setw(TIME_WIDTH) << "Ryser"; break;
        case LAPLAS: cout << setw(TIME_WIDTH) << "Laplas"; break;
        case GLYNN: cout << setw(TIME_WIDTH) << "Glynn"; break;
        case EJ: cout << setw(TIME_WIDTH) << "E.Bax J.Frankin"; break;
        case GG: cout << setw(TIME_WIDTH) << "Godsil-Gotman"; break;
        case AB: cout << setw(TIME_WIDTH) << "A.Barvinok"; break;
        case MC: cout << setw(TIME_WIDTH) << "P.McCullagh"; break;
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
    for (size_t j = 0; j < size; j++) {
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
        double avg = time_sums[m] / size;
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
        case EJ: cout << setw(TIME_WIDTH) << "E.Bax J.Frankin"; break;
        case GG: cout << setw(TIME_WIDTH) << "Godsil-Gotman"; break;
        case AB: cout << setw(TIME_WIDTH) << "A.Barvinok"; break;
        case MC: cout << setw(TIME_WIDTH) << "P.McCullagh"; break;
        }
    }
    cout << "\n";

    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    for (size_t j = 0; j < size; j++) {
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
void the_table(const vector<Method>& methods, size_t size, const vector<vector<double>>& results, const vector<vector<double>>& times) {
    const int NAME_WIDTH = 15;
    const int TIME_WIDTH = 12;

    cout << "\nTIME (seconds)\n";
    cout << left << setw(NAME_WIDTH) << "Matrix";
    for (auto& method : methods) {
        switch (method) {
        case RYSER: cout << setw(TIME_WIDTH) << "Ryser"; break;
        case LAPLAS: cout << setw(TIME_WIDTH) << "Laplas"; break;
        case GLYNN: cout << setw(TIME_WIDTH) << "Glynn"; break;
        case EJ: cout << setw(TIME_WIDTH) << "E.Bax J.Frankin"; break;
        case GG: cout << setw(TIME_WIDTH) << "Godsil-Gotman"; break;
        case AB: cout << setw(TIME_WIDTH) << "A.Barvinok"; break;
        case MC: cout << setw(TIME_WIDTH) << "P.McCullagh"; break;
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
    for (size_t j = 0; j < size; j++) {
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
        double avg = time_sums[m] / size;
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
        case EJ: cout << setw(TIME_WIDTH) << "E.Bax J.Frankin"; break;
        case GG: cout << setw(TIME_WIDTH) << "Godsil-Gotman"; break;
        case AB: cout << setw(TIME_WIDTH) << "A.Barvinok"; break;
        case MC: cout << setw(TIME_WIDTH) << "P.McCullagh"; break;
        }
    }
    cout << "\n";

    cout << string(NAME_WIDTH, '-') << "+";
    for (size_t i = 0; i < methods.size(); i++) {
        cout << string(TIME_WIDTH, '-');
        if (i < methods.size() - 1) cout << "+";
    }
    cout << endl;

    for (size_t j = 0; j < size; j++) {
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

bool has_negative(const Matrix& A) {
    size_t n1 = A.size();
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n1; j++) {
            if (A[i][j] < 0) {
                return true;
            }
        }
    }
    return false;
}
bool has_negative(const DoubleMatrix& A) {
    size_t n1 = A.size();
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n1; j++) {
            if (A[i][j] < 0) {
                return true;
            }
        }
    }
    return false;
}
bool is_01(const Matrix& A) {
    size_t n1 = A.size();
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n1; j++) {
            if (A[i][j] != 0 && A[i][j] != 1) {
                return false;
            }
        }
    }
    return true;
}
bool is_01(const DoubleMatrix& A, double prec) {
    size_t n1 = A.size();
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n1; j++) {
            if (fabs(A[i][j]) > prec && fabs(A[i][j] - 1) > prec) {
                return false;
            }
        }
    }
    return true;
}

void build_factorial() {
    fact[0] = 1;
    for (int i = 1; i <= factN; ++i) {
        fact[i] = fact[i - 1] * i;
    }
}