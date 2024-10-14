#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

using namespace std;

class Polynomial {
private:
    vector<double> coeffs; // Store coefficients of the polynomial
public:
    // Constructors
    Polynomial() : coeffs(1, 0.0) {}
    Polynomial(const vector<double>& coefficients) : coeffs(coefficients) {}
    Polynomial(const Polynomial& other) : coeffs(other.coeffs) {}

    // Destructor
    ~Polynomial() {}

    // Assignment operator
    Polynomial& operator=(const Polynomial& other) {
        if (this != &other) {
            coeffs = other.coeffs;
        }
        return *this;
    }

    // Arithmetic operators
    Polynomial operator+(const Polynomial& other) const {
        int maxDegree = max(degree(), other.degree());
        vector<double> result(maxDegree + 1, 0.0);
        for (int i = 0; i <= degree(); ++i) {
            result[i] += coeffs[i];
        }
        for (int i = 0; i <= other.degree(); ++i) {
            result[i] += other.coeffs[i];
        }
        return Polynomial(result);
    }

    Polynomial operator-(const Polynomial& other) const {
        int maxDegree = max(degree(), other.degree());
        vector<double> result(maxDegree + 1, 0.0);
        for (int i = 0; i <= degree(); ++i) {
            result[i] += coeffs[i];
        }
        for (int i = 0; i <= other.degree(); ++i) {
            result[i] -= other.coeffs[i];
        }
        return Polynomial(result);
    }

    Polynomial operator*(const Polynomial& other) const {
        vector<double> result(degree() + other.degree() + 1, 0.0);
        for (int i = 0; i <= degree(); ++i) {
            for (int j = 0; j <= other.degree(); ++j) {
                result[i + j] += coeffs[i] * other.coeffs[j];
            }
        }
        return Polynomial(result);
    }

    // Equality operator
    bool operator==(const Polynomial& other) const {
        return coeffs == other.coeffs;
    }

    // Output operator
    friend ostream& operator<<(ostream& out, const Polynomial& poly) {
        for (int i = poly.degree(); i >= 0; --i) {
            out << poly.coeffs[i];
            if (i > 0) out << "x^" << i << " + ";
        }
        return out;
    }

    // Utility functions
    int degree() const {
        return coeffs.size() - 1;
    }

    double evaluate(double x) const {
        double result = 0.0;
        for (int i = degree(); i >= 0; --i) {
            result = result * x + coeffs[i];
        }
        return result;
    }

    Polynomial & operator+(__gnu_cxx::__alloc_traits<allocator<double>>::value_type value) const;

    Polynomial compose(const Polynomial& q) const {
        Polynomial result;
        for (int i = degree(); i >= 0; --i) {
            result = result * q + coeffs[i];
        }
        return result;
    }

    Polynomial derivative() const {
        if (degree() == 0) return Polynomial();
        vector<double> derivCoeffs(degree());
        for (int i = 1; i <= degree(); ++i) {
            derivCoeffs[i - 1] = coeffs[i] * i;
        }
        return Polynomial(derivCoeffs);
    }

    Polynomial integral() const {
        vector<double> integralCoeffs(degree() + 2, 0.0);
        for (int i = 0; i <= degree(); ++i) {
            integralCoeffs[i + 1] = coeffs[i] / (i + 1);
        }
        return Polynomial(integralCoeffs);
    }

    double integral(double x1, double x2) const {
        Polynomial integralPoly = integral();
        return integralPoly.evaluate(x2) - integralPoly.evaluate(x1);
    }

    double getRoot(double guess = 1, double tolerance = 1e-6, int maxIter = 100) const {
        Polynomial derivPoly = derivative();
        double x0 = guess;
        for (int i = 0; i < maxIter; ++i) {
            double fx = evaluate(x0);
            double dfx = derivPoly.evaluate(x0);
            if (abs(fx) < tolerance) break;
            x0 = x0 - fx / dfx;
        }
        return x0;
    }

    void setCoefficients(const vector<double>& coefficients) {
        coeffs = coefficients;
    }

    double getCoefficient(int degree) const {
        return degree < coeffs.size() ? coeffs[degree] : 0.0;
    }
};

// Test case function provided by the user
void runTests() {
    // Test Case 1: Default constructor and degree
    Polynomial p1;
    assert(p1.degree() == 0); // p1 should be the zero polynomial

    // Test Case 2: Constructor with coefficients
    Polynomial p2({1, 2, 3}); // p2 represents 1 + 2x + 3x^2
    assert(p2.degree() == 2);

    // Test Case 3: Addition of polynomials
    Polynomial p3({1, 2, 3});
    Polynomial p4({3, 4, 5}); // p3 + p4 = (4, 6, 8)
    Polynomial pAdd = p3 + p4;
    assert(pAdd.getCoefficient(0) == 4 && pAdd.getCoefficient(1) == 6 && pAdd.getCoefficient(2) == 8);

    // Test Case 4: Subtraction of polynomials
    Polynomial pSub = p4 - p3;
    assert(pSub.getCoefficient(0) == 2 && pSub.getCoefficient(1) == 2 && pSub.getCoefficient(2) == 2);

    // Test Case 5: Multiplication of polynomials
    Polynomial pMul = p3 * p4; // (1 + 2x + 3x^2) * (3 + 4x + 5x^2)
    assert(pMul.getCoefficient(0) == 3 && pMul.getCoefficient(1) == 10 && pMul.getCoefficient(2) == 22);

    // Test Case 6: Polynomial evaluation
    double val = p3.evaluate(2); // 1 + 2*2 + 3*2^2 = 1 + 4 + 12 = 17
    assert(val == 17);

    // Test Case 7: Polynomial composition
    Polynomial pComp = p3.compose(p2); // p3(p2), composing p3 and p2
    assert(pComp.getCoefficient(0) == 1); // Check some coefficients

    // Test Case 8: Polynomial derivative
    Polynomial pDeriv = p3.derivative(); // Derivative of 1 + 2x + 3x^2 = 2 + 6x
    assert(pDeriv.getCoefficient(0) == 2 && pDeriv.getCoefficient(1) == 6);

    // Test Case 9: Polynomial integral (indefinite)
    Polynomial pInt = p3.integral(); // Integral of 1 + 2x + 3x^2 = x + x^2 + x^3
    assert(pInt.getCoefficient(0) == 0 && pInt.getCoefficient(1) == 1 && pInt.getCoefficient(2) == 1 && pInt.getCoefficient(3) == 1);

    // Test Case 10: Polynomial integral (definite)
    double intVal = p3.integral(0, 1); // Definite integral from 0 to 1
    assert(fabs(intVal - 2.5) < 1e-6); // Value is approximately 2.5

    // Test Case 11: Finding roots (Newton's method)
    Polynomial pRoot({-6, 1}); // x - 6, root at x = 6
    double root = pRoot.getRoot(0);
    assert(fabs(root - 6) < 1e-6);

    // Test Case 12: Copy constructor
    Polynomial pCopy(p3);
    assert(pCopy == p3);

    // Test Case 13: Assignment operator
    Polynomial pAssign;
    pAssign = p3;
    assert(pAssign == p3);

    // Test Case 14: Zero polynomial
    Polynomial pZero({0, 0, 0});
    assert(pZero.degree() == 0);

    // Test Case 15: Set coefficients
    pZero.setCoefficients({1, 0, -4});
    assert(pZero.getCoefficient(0) == 1 && pZero.getCoefficient(2) == -4);

    // Test Case 16: Get coefficient
    assert(p3.getCoefficient(2) == 3);

    // Test Case 17: Compare polynomials (equality)
    Polynomial pEqual1({2, 3, 4});
    Polynomial pEqual2({2, 3, 4});
    assert(pEqual1 == pEqual2);

    // Test Case 18: Multiplying by zero polynomial
       // Test Case 18: Multiplying by zero polynomial
    Polynomial pZeroMul = p3 * Polynomial({0});
    assert(pZeroMul == Polynomial({0}));

    // Test Case 19: Adding zero polynomial
    Polynomial pAddZero = p3 + Polynomial({0});
    assert(pAddZero == p3);

    // Test Case 20: Subtracting zero polynomial
    Polynomial pSubZero = p3 - Polynomial({0});
    assert(pSubZero == p3);

    // Test Case 21: Polynomial with all negative coefficients
    Polynomial pNeg({-1, -2, -3});
    Polynomial pNegSum = pNeg + Polynomial({1, 2, 3});
    assert(pNegSum == Polynomial({0, 0, 0}));

    // Test Case 22: Large degree polynomial
    Polynomial pLarge({1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}); // x^10 + 1
    assert(pLarge.degree() == 10);

    // Test Case 23: Multiplying large degree polynomials
    Polynomial pLargeMul = pLarge * pLarge;
    assert(pLargeMul.degree() == 20);

    // Test Case 24: Integral of a constant polynomial
    Polynomial pConst({3});
    double constInt = pConst.integral(0, 2); // Integral of 3 from 0 to 2 is 6
    assert(constInt == 6);

    // Test Case 25: Derivative of a constant polynomial
    Polynomial pConstDeriv = pConst.derivative();
    assert(pConstDeriv == Polynomial({0}));

    // Test Case 26: Composition with zero polynomial
    Polynomial pCompZero = p3.compose(Polynomial({0}));
    assert(pCompZero == Polynomial({1}));

    // Test Case 27: Large coefficients in polynomial
    Polynomial pLargeCoeff({1e9, -1e9});
    assert(pLargeCoeff.getCoefficient(0) == 1e9 && pLargeCoeff.getCoefficient(1) == -1e9);

    // Test Case 28: Evaluate polynomial at zero
    assert(p3.evaluate(0) == 1);

    // Test Case 29: Evaluate polynomial at one
    assert(p3.evaluate(1) == 6); // 1 + 2 + 3 = 6

    // Test Case 30: Polynomial with fractional coefficients
    Polynomial pFrac({0.5, 1.5, -2.5});
    assert(fabs(pFrac.evaluate(2) - (-5.5)) < 1e-6); // Evaluate at 2
}

int main() {
    runTests();
    cout << "All tests passed!" << endl;
    return 0;
}
