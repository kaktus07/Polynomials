/*
"Problem solving with C++" by Walter Savitch, 10th edition,  chapter 11, pp.731-732:

6. Using dynamic arrays, implement a polynomial class with polynomial ad-
dition, subtraction, and multiplication.
Discussion: A variable in a polynomial does very little other than act as
a placeholder for the coefficients. Hence, the only interesting thing
about polynomials is the array of coefficients and the corresponding
exponent. Think about the polynomial
x*x*x + x + 1
One simple way to implement the polynomial class is to use an array
of doubles to store the coefficients. The index of the array is the expo-
nent of the corresponding term. Where is the term in x*x in the previ-
ous example? If a term is missing, then it simply has a zero coefficient.
There are techniques for representing polynomials of high degree
with many missing terms. These use so-called sparse polynomial tech-
niques. Unless you already know these techniques, or learn very
quickly, don’t use them.
Provide a default constructor, a copy constructor, and a parameterized
constructor that enable an arbitrary polynomial to be constructed. Also
supply an overloaded operator = and a destructor.
Provide these operations:
■■ polynomial + polynomial
■■ constant + polynomial
■■ polynomial + constant
■■ polynomial − polynomial
■■ constant − polynomial
■■ polynomial − constant
■■ polynomial * polynomial
■■ constant * polynomial
■■ polynomial * constant
Supply functions to assign and extract coefficients, indexed by exponent.
Supply a function to evaluate the polynomial at a value of type double.
You should decide whether to implement these functions as members,
friends, or stand-alone functions.
*/

//Rafał Józef Gaida
//July 9th, 2023 - August 5th, 2023

/*
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
*/

#include<iostream>
#include<cmath>
#include<vector>
#include<stdexcept>
#include<cassert>
#include<complex>

class Polynomial
{
public:
    //Constructors
    Polynomial();      //Default constructor
    Polynomial(int arrSize); //empty array constructor with 'arrSize' size
    Polynomial(double zero, double first, double second, double third);    //Parameterized constructor
    Polynomial(const Polynomial& copiedClass);    //Copy constructor
    //Destructor
    ~Polynomial();

    //Operators overloading
    friend Polynomial operator +(const Polynomial& first, const Polynomial& second);
    friend Polynomial operator +(const double first, const Polynomial& second);
    friend Polynomial operator +(const Polynomial& first, const double second);
    friend Polynomial operator -(const Polynomial& first, const Polynomial& second);
    friend Polynomial operator -(const double first, const Polynomial& second);
    friend Polynomial operator -(const Polynomial& first, const double second);
    friend Polynomial operator -(const Polynomial& first);
    friend Polynomial operator *(const Polynomial& first, const Polynomial& second);
    friend Polynomial operator *(const double first, const Polynomial& second);
    friend Polynomial operator *(const Polynomial& first, const double second);
    friend Polynomial operator /(const Polynomial& first, const double second);
    friend Polynomial operator /(const Polynomial& first, const Polynomial& second);

    //Member functions
    void show() const;
    void declare(int position, double value);
    void showRoots(std::vector<double> vect) const;
    Polynomial& operator =(const Polynomial& rightSide);
    double evaluate(double x) const;
    Polynomial differentiate() const;
    Polynomial integrate() const;
    std::vector<double> roots() const;

private:
    double* poly;
    int arraySize;
};

Polynomial operator^(const Polynomial& base, int exponent);

int main()
{
    //Example driver program

    // Create some polynomials using different constructors
    Polynomial p1;  // Default constructor
    Polynomial p2(3);  // Empty array constructor with size 4
    Polynomial p3(1.0, -2.0, 0.5, 3.0); // Parameterized constructor

    // Declare some values in the polynomials
    p2.declare(0, 1.0);
    p2.declare(2, 2.0);

    // Show the polynomials
    std::cout << "Polynomial p1: ";
    p1.show();
    std::cout << "Polynomial p2: ";
    p2.show();
    std::cout << "Polynomial p3: ";
    p3.show();

    // Addition
    Polynomial resultAddition = p1 + p2;
    std::cout << "p1 + p2: ";
    resultAddition.show();

    // Subtraction
    Polynomial resultSubtraction = p1 - p2;
    std::cout << "p1 - p2: ";
    resultSubtraction.show();

    // Negation
    Polynomial resultNegation = -p3;
    std::cout << "-p3: ";
    resultNegation.show();

    // Multiplication
    Polynomial resultMultiplication = p2 * p3;
    std::cout << "p2 * p3: ";
    resultMultiplication.show();

    // Division
    Polynomial resultFactorDivision = p3 / 2.0;
    std::cout << "p3 / 2.0: ";
    resultFactorDivision.show();

    // Division
    Polynomial resultPolynomialDivision = p3 / p2;
    std::cout << "p3 / p2: ";
    resultPolynomialDivision.show();

    // Evaluation
    double x = 1.5;
    std::cout << "p3 evaluated at x = " << x << ": " << p3.evaluate(x) << std::endl;

    // Differentiation
    Polynomial resultDifferentiation = p3.differentiate();
    std::cout << "Derivative of p3: ";
    resultDifferentiation.show();

    // Integration
    Polynomial resultIntegration = p2.integrate();
    std::cout << "Integral of p2: ";
    resultIntegration.show();

    // Roots
    std::vector<double> roots = p3.roots();
    std::cout << "Roots of p3: ";
    p3.showRoots(roots);

    // Assigning p2 to p1
    p1 = p2;
    std::cout << "After assignment, Polynomial p1: ";
    p1.show();

    // Calculate exponent of a polynomial
    int exponent = 3;
    Polynomial resultExponent = p3 ^ exponent;
    std::cout << "p3 raised to the power of " << exponent << ": ";
    resultExponent.show();

    return 0;
}

//Initializers and destructors definitions
Polynomial::Polynomial()
{
    arraySize = 3;
    poly = new double[arraySize];
    poly[0] = 1;
    poly[1] = 2;
    poly[2] = 1;
}

Polynomial::Polynomial(int arrSize)
{
    arraySize = arrSize;
    poly = new double[arraySize];
    for (int i=0; i<arraySize; i++)
        poly[i] = 0;
}

Polynomial::Polynomial(double zero, double first, double second, double third)
{
    arraySize = 4;
    poly = new double[arraySize];
    poly[0] = zero;
    poly[1] = first;
    poly[2] = second;
    poly[3] = third;
}

Polynomial::Polynomial(const Polynomial& copiedClass)
{
    arraySize = copiedClass.arraySize;
    poly = new double[arraySize];
    for (int i=0; i<arraySize; i++)
        poly[i] = copiedClass.poly[i];
}

Polynomial::~Polynomial()
{
    delete[] poly;
    //std::cout<<"Memory cleared. Exiting the class!\n";
}

//Operators overloading
Polynomial operator +(const Polynomial& first, const Polynomial& second)
{
    int arraySize = first.arraySize;
    if (second.arraySize > first.arraySize)
        arraySize = second.arraySize;

    Polynomial temp(arraySize);
    for (int i=0; i<arraySize; i++)
        temp.poly[i] = first.poly[i] + second.poly[i];
    //std::cout<<"\nArray size: "<<arraySize<<std::endl;
    return temp;
}
Polynomial operator +(const double first, const Polynomial& second)
{
    Polynomial temp(second.arraySize);
    temp.poly[0] = first + second.poly[0];
    for (int i=1; i<second.arraySize; i++)
        temp.poly[i] = second.poly[i];

    return temp;
}
Polynomial operator +(const Polynomial& first, const double second)
{
    Polynomial temp(first.arraySize);
    temp.poly[0] = first.poly[0] + second;
    for (int i=1; i<first.arraySize; i++)
        temp.poly[i] = first.poly[i];

    return temp;
}
Polynomial operator -(const Polynomial& first, const Polynomial& second)
{
    int arraySize = first.arraySize;
    if (second.arraySize > first.arraySize)
        arraySize = second.arraySize;

    Polynomial temp(arraySize);
    for (int i=0; i<arraySize; i++)
        temp.poly[i] = first.poly[i] - second.poly[i];

    return temp;
}
Polynomial operator -(const double first, const Polynomial& second)
{
    Polynomial temp(second.arraySize);
    temp.poly[0] = first - second.poly[0];
    for (int i=1; i<second.arraySize; i++)
        temp.poly[i] = second.poly[i];

    return temp;
}
Polynomial operator -(const Polynomial& first, const double second)
{
    Polynomial temp(first.arraySize);
    temp.poly[0] = first.poly[0] - second;
    for (int i=1; i<first.arraySize; i++)
        temp.poly[i] = first.poly[i];

    return temp;
}
Polynomial operator -(const Polynomial& first)
{
    Polynomial temp(first.arraySize);
    for (int i=0; i<first.arraySize; i++)
        temp.poly[i] = -first.poly[i];

    return temp;
}
Polynomial operator *(const Polynomial& first, const Polynomial& second)
{
    int arraySize = (first.arraySize-1) + (second.arraySize-1) + 1;
    Polynomial temp(arraySize);
    for (int i=0; i<first.arraySize; i++)
        for(int j=0; j<second.arraySize; j++)
            temp.poly[i+j] = temp.poly[i+j] + first.poly[i] * second.poly[j];

    return temp;
}
Polynomial operator *(const double first, const Polynomial& second)
{
    Polynomial temp(second.arraySize);
    for (int i=0; i<second.arraySize; i++)
        temp.poly[i] = first * second.poly[i];

    return temp;
}
Polynomial operator *(const Polynomial& first, const double second)
{
    Polynomial temp(first.arraySize);
    for (int i=0; i<first.arraySize; i++)
        temp.poly[i] = first.poly[i] * second;

    return temp;
}

Polynomial operator /(const Polynomial& first, const double second)
{
    if (second ==0)
        throw std::out_of_range("Do not divide by zero");
    else
    {
        Polynomial temp(first.arraySize);
        for (int i=0; i<first.arraySize; i++)
            temp.poly[i] = first.poly[i] / second;

        return temp;
    }
}

Polynomial operator /(const Polynomial& first, const Polynomial& second)
{
    if (second.arraySize > first.arraySize)
        throw std::out_of_range("The order of denominator is higher than the order of numerator. Can not perform division!\n");

    int denominatorDegree=0;
    for (int i = 0; i < second.arraySize; i++)
    {
        if (second.poly[i] != 0)
            denominatorDegree++;
    }

    if (denominatorDegree == 0)
        throw std::out_of_range("The denominator polynomial equals zero. Do not divide by zero!\n");

    int highestExponent = second.arraySize-1;
    double highestFactor = second.poly[highestExponent];

    Polynomial result(first.arraySize-highestExponent); //result stores the result
    Polynomial numerator = first;   //numerator stores the current 'first' during division&subtraction

    for (int i=numerator.arraySize-1; i>=highestExponent; i--)
    {
        result.poly[i-highestExponent] = numerator.poly[i] / highestFactor;
        //std::cout << result.poly[i];

        Polynomial temp(i-highestExponent+1);  //temp stores current "factor&exponent"
        temp.declare(i-highestExponent, result.poly[i-highestExponent]);

        numerator = numerator - (temp * second);
        numerator.arraySize--;
    }

    return result;
}

void Polynomial::show() const
{
    std::cout<<"y = ";
    for(int i=0; i<arraySize; i++)
    {
        if (poly[i] !=0)
        {
            if (poly[i] < 0)
                std::cout<<poly[i]<<"*x^"<<i<<" ";
            else if(poly[i] > 0)
                std::cout<<"+"<<poly[i]<<"*x^"<<i<<" ";
        }
    }
    std::cout<<std::endl;
}

void Polynomial::declare(int position, double value)
{
    if (position >= arraySize)
    {
        std::cout<<"Array size exceeded! Closing the program...\n";
        throw std::out_of_range("Invalid position for polynomial coefficient modification.");
    }
    else
        poly[position] = value;
}

void Polynomial::showRoots(std::vector<double> vect) const
{
    std::cout << "The approximate roots are: ";
    if (vect.empty())
        std::cout << "none" << std::endl;

    else
    {
    std::cout << vect[0];
    for (int i=1; i<vect.size(); i++)
        std::cout << "; " << vect[i];

    std::cout << std::endl;
    }
}

Polynomial& Polynomial::operator =(const Polynomial& rightSide)
{
    if (this != &rightSide)
    {
        delete[] poly;
        arraySize = rightSide.arraySize;
        poly = new double[arraySize];

        for (int i=0; i<arraySize; i++)
            poly[i] = rightSide.poly[i];
    }

    return *this;
}

//Requires cmath library
double Polynomial::evaluate(double x) const
{
    double result=0;
    for (int i=0; i<arraySize; i++)
        result = result + (poly[i] * pow(x,i));

    return result;
}

Polynomial Polynomial::differentiate() const
{
    Polynomial result(arraySize-1);

    for (int i=1; i<arraySize; i++)
        result.poly[i-1] = poly[i] * i;

    return result;
}

Polynomial Polynomial::integrate() const
{
    Polynomial integral(arraySize+1);
    integral.poly[0] = 0;  //technically should be a constant C

    for (int i=0; i<arraySize; i++)
        integral.poly[i+1] = poly[i] / (i+1);

    return integral;
}

std::vector<double> Polynomial::roots() const
{
    std::vector<double> roots;

    //Handle special cases: degree 0 or 1
    if (arraySize <= 1)
    {
        if (poly[0] != 0)
            roots.push_back(NAN);
    }

    // Handle degree 2 polynomial
    else if (arraySize == 2)
        roots.push_back(-poly[0] / poly[1]);

    //Handle degree 3 polynomial (quadratic)
    else if (arraySize == 3)
    {
        double discriminant = pow(poly[1], 2) - 4*poly[0]*poly[2];
        if (discriminant < 0)
            throw std::out_of_range("No real roots!\n");

        roots.push_back((-poly[1] - sqrt(discriminant))/(2*poly[2]));
        roots.push_back((-poly[1] + sqrt(discriminant))/(2*poly[2]));
    }

    //Handle higher degrees with a numerical approach
    else
    {
    const double min_x = -100.0; // Minimum x value to search for roots
    const double max_x = 100.0; // Maximum x value to search for roots
    const double step = 0.001; // Step size for iterating over x values
    const double tolerance = 0.002; // Tolerance for considering a value as an approximate root

        for (double x = min_x; x <= max_x; x += step)
        {
            double fx = evaluate(x);

            if (std::abs(fx) < tolerance)
            {
                // Found an approximate root within the specified tolerance
                roots.push_back(x);
            }
        }
    }

    return roots;
}

Polynomial operator^(const Polynomial& base, int exponent)
{
    if (exponent < 1)
        throw std::range_error("The exponent has to be a positive integral!\n");
    else if (exponent == 1)
        return base;
    else
    {
        Polynomial temp;
        temp = base;
        for (int i=1; i<exponent; i++)
            temp = temp * base;
        return temp;
    }
}
