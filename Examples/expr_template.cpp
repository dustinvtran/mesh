#include <iostream>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"

// USE_EXPR_TEMP == 0  --  Default to the trivial op* and op+
// USE_EXPR_TEMP == 1  --  Use lazy evaluation to prevent temporaries
#define USE_EXPR_TEMP 1


/** Curiously Recurring Template Pattern (CRTP) for Matrix Expressions
 * A Matrix Expression is anything that inherits from MatExpr<E> and implements:
 * concept E {
 *  double operator()(unsigned i, unsigned j) const;
 *  unsigned rows() const;
 *  unsigned cols() const;
 * }
 * A MatExpr presents this interface, but simply forwards to the derived class.
 * TODO: Should rename the derived/base methods to prevent override confusion.
 */
template <typename E>
struct MatExpr {
  // Cast to the derived type
  const E& derived() const {
    return static_cast<const E&>(*this);
  }
  double operator()(unsigned i, unsigned j) const {
    return derived()(i,j);
  }
  unsigned rows() const {
    return derived().rows();
  }
  unsigned cols() const {
    return derived().cols();
  }
  // Print out type information!
  friend std::ostream& operator<<(std::ostream& s, const MatExpr& m) {
    return s << m.derived();
  }
  // No member data! Only interface.
};


/** Lazy evaluation for a scaled matrix expression */
template <typename E1>
struct MatScale : public MatExpr<MatScale<E1>> {
  MatScale(double b, const MatExpr<E1>& a)
      : alpha_(b), a_(a) {
  }
  unsigned rows() const {
    return a_.rows();
  }
  unsigned cols() const {
    return a_.cols();
  }
  double operator()(unsigned i, unsigned j) const {
    return alpha_ * a_(i,j);
  }
  // Print out type information!
  friend std::ostream& operator<<(std::ostream& s, const MatScale& m) {
    return s << "MatScale<" << m.alpha_ << "," << m.a_ << ">";
  }
private:
  double alpha_;
  const MatExpr<E1>& a_;
};


/** Lazy evaluation for the sum of two matrix expressions */
template <typename E1, typename E2>
struct MatAdd : public MatExpr<MatAdd<E1, E2>> {
  MatAdd(const MatExpr<E1>& a, const MatExpr<E2>& b)
      : a_(a), b_(b) {
    assert(a_.rows() == b_.rows() && a_.cols() == b_.cols());
  }
  unsigned rows() const {
    return a_.rows();
  }
  unsigned cols() const {
    return b_.cols();
  }
  double operator()(unsigned i, unsigned j) const {
    return a_(i,j) + b_(i,j);
  }
  // Print out type information!
  friend std::ostream& operator<<(std::ostream& s, const MatAdd& m) {
    return s << "MatAdd<" << m.a_ << "," << m.b_ << ">";
  }
private:
  const MatExpr<E1>& a_;
  const MatExpr<E2>& b_;
};


/** Lazy evaluation for the product of two matrix expressions */
template <typename E1, typename E2>
struct MatMult : public MatExpr<MatMult<E1, E2>> {
  MatMult(const MatExpr<E1>& a, const MatExpr<E2>& b)
      : a_(a), b_(b) {
    assert(a_.cols() == b_.rows());
  }
  unsigned rows() const {
    return a_.rows();
  }
  unsigned cols() const {
    return b_.cols();
  }
  double operator()(unsigned i, unsigned j) const {
    double result = 0;
    for (unsigned k = 0; k < a_.cols(); ++k)
      result += a_(i,k) * b_(k,j);
    return result;
  }
  // Print out type information!
  friend std::ostream& operator<<(std::ostream& s, const MatMult& m) {
    return s << "MatMult<" << m.a_ << "," << m.b_ << ">";
  }
private:
  const MatExpr<E1>& a_;
  const MatExpr<E2>& b_;
};


/** Concrete Matrix of values -- Also a matrix expression */
struct Matrix : public MatExpr<Matrix> {
  // Construct a Matrix
  Matrix(unsigned rows, unsigned cols, double val = 0)
      : rows_(rows), cols_(cols), v_(rows_*cols_, val) {
  }
  // Construct a Matrix from a MatExpr
  template <typename E>
  Matrix(const MatExpr<E>& m)
      : rows_(m.rows()), cols_(m.cols()) {
    v_.reserve(rows() * cols());
    for (unsigned i = 0; i < rows(); ++i)
      for (unsigned j = 0; j < cols(); ++j)
        v_.push_back(m(i,j));
  }
  unsigned rows() const {
    return rows_;
  }
  unsigned cols() const {
    return cols_;
  }
  const double& operator()(unsigned i, unsigned j) const {
    return v_[i*cols() + j];
  }
  double& operator()(unsigned i, unsigned j) {
    return v_[i*cols() + j];
  }
  // Print out type information!
  friend std::ostream& operator<<(std::ostream& s, const Matrix&) {
    return s << "Matrix";
  }
 private:
  unsigned rows_, cols_;
  std::vector<double> v_;
};


//////////////////////
// MATRIX OPERATORS //
//////////////////////

#if USE_EXPR_TEMP

/** Add any two matrix expressions */
template <typename E1, typename E2>
MatAdd<E1,E2> operator+(const MatExpr<E1>& a, const MatExpr<E2>& b) {
  return MatAdd<E1,E2>(a, b);
}

/** Multiply any two matrix expressions */
template <typename E1, typename E2>
MatMult<E1,E2> operator*(const MatExpr<E1>& a, const MatExpr<E2>& b) {
  return MatMult<E1,E2>(a, b);
}

/** Scale any matrix expression */
template <typename E>
MatScale<E> operator*(double b, const MatExpr<E>& a) {
  return MatScale<E>(b, a);
}
template <typename E>
MatScale<E> operator*(const MatExpr<E>& a, double b) {
  return b * a;
}

#else

// Classic op+ for two Matrices
Matrix operator+(const Matrix& a, const Matrix& b) {
  Matrix r(a.rows(), b.cols());
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      r(i,j) = a(i,j) + b(i,j);
  return r;
}

// Classic op* for two Matrices
Matrix operator*(const Matrix& a, const Matrix& b) {
  Matrix r(a.rows(), b.cols(), 0);
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      for (unsigned k = 0; k < a.cols(); ++k)
        r(i,j) += a(i,k) * b(k,j);
  return r;
}

// Classic op* for a scalar and a Matrix
Matrix operator*(double b, const Matrix& a) {
  Matrix r(a.rows(), a.cols());
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      r(i,j) = b * a(i,j);
  return r;
}
Matrix operator*(const Matrix& a, double b) {
  return b * a;
}

#endif


int main() {
  CS207::Clock timer;
  unsigned N = 512;
  Matrix A(N, N, 1.23);
  Matrix B(N, N, 2);
  Matrix C(N, N, 5.76);
  double alpha = 3.14;

  // Print out the type information of M
  std::cout << "Type of Expr:\n  "
            << 2.73*(A*B) + alpha*(alpha*A + B + A(0,0)*C) << std::endl;

  timer.start();
  Matrix D = 2.73*(A*B) + alpha*(alpha*A + B + A(0,0)*C);
  double eval_time = timer.seconds();
  std::cout << "Compute Time: " << eval_time << " seconds" << std::endl;
}
