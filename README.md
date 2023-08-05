# Polynomials
A c++ class for calculating polynomials

The Polynomial class is designed to represent and perform various operations on polynomials.
A polynomial is an algebraic expression consisting of variables (usually represented as 'x') raised to non-negative integer powers, multiplied by coefficients.
The Polynomial class allows users to create, manipulate, and evaluate polynomials.

Functionality Overview

  Constructors:
    Default Constructor: Creates a polynomial with coefficients initialized to [1, 2, 1], representing the polynomial x^2 + 2x + 1.
    Empty Array Constructor: Creates a polynomial with an empty array of coefficients and a given arrSize.
    Parameterized Constructor: Creates a polynomial with given coefficients.
    Copy Constructor: Creates a copy of an existing polynomial.

  Destructors:
    The destructor is responsible for freeing the dynamically allocated memory for the polynomial coefficients.

  Operators Overloading:
    Addition (+), Subtraction (-), Multiplication (*), and Division (/) of polynomials.
    Addition and Subtraction with scalar values.
    Negation (-) of a polynomial.
    Exponentiation (^) of a polynomial with a non-negative integer exponent.

  Member Functions:
    show(): Displays the polynomial in a readable format, e.g., y = 2*x^2 - x + 3.
    declare(int position, double value): Assigns a coefficient value to the given exponent position in the polynomial.
    showRoots(std::vector<double> vect) const: Displays the roots of the polynomial.
    evaluate(double x) const: Evaluates the polynomial for a given value of x.
    differentiate() const: Calculates the derivative of the polynomial.
    integrate() const: Calculates the indefinite integral of the polynomial.
    roots() const: Finds the roots (zeros) of the polynomial (with a brute-force algorithm)

Limits and Considerations

  Exponents: The class supports integer exponents (non-negative) for polynomial terms.
  Coefficients: The coefficients are represented as double values, allowing for real numbers.
  Division by Zero: Division by a polynomial with zero coefficients (degree 0) is not supported and will throw an std::out_of_range exception.
  Numerical Precision: The class uses floating-point arithmetic for coefficients and evaluations, which may result in numerical precision errors for large or small coefficient values.
  Real Roots: The roots() function returns real roots of the polynomial only. Complex roots are not handled. The precision of finding the roots is limited.
  Exponentiation: The exponentiation operator supports only non-negative integer exponents. Negative or fractional exponents are not supported.
