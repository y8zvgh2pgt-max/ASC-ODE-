#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <autodiff.hpp>

using namespace std;
using namespace ASC_ode;

class PendulumAD : public NonlinearFunction
{
  double m_length;
  double m_gravity;

public:
  PendulumAD(double length, double gravity=9.81) 
    : m_length(length), m_gravity(gravity) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    T_evaluate<double>(x, f);
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    Vector<AutoDiff<2>> x_ad(2);
    Vector<AutoDiff<2>> f_ad(2);

    x_ad(0) = Variable<0>(x(0)); 
    x_ad(1) = Variable<1>(x(1));

    T_evaluate<AutoDiff<2>>(x_ad, f_ad);

    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++)
         df(i,j) = f_ad(i).deriv()[j];
  }

  template <typename T>
  void T_evaluate (VectorView<T> x, VectorView<T> f) const
  {
    f(0) = x(1);

    f(1) = - T(m_gravity / m_length) * sin(x(0));
  }
};

int main()
{
  double l = 1.0;
  double tend = 10.0;
  int steps = 1000;
  double tau = tend/steps;

  Vector<> y = { 3.0, 0.0 }; 
  
  auto rhs = make_shared<PendulumAD>(l);
  
  CrankNicolson stepper(rhs);

  ofstream out("output_pendulum.txt");
  
  for (int i = 0; i < steps; i++)
  {
     stepper.DoStep(tau, y);
     out << (i+1)*tau << " " << y(0) << " " << y(1) << endl;
  }
}