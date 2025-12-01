#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;
using namespace std;

class RCCircuit : public NonlinearFunction
{
  double R;
  double C;
  
public:
  RCCircuit(double resistance, double capacity) 
    : R(resistance), C(capacity) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }

  double U0(double t) const 
  {
    return std::cos(100* M_PI *t); 
  }

  double U0_prime(double t) const
  {
    return -100.0 * M_PI * std::sin(100 * M_PI * t);
  }

  void evaluate (VectorView<double> y, VectorView<double> f) const override
  {
    double Uc = y(0);
    double t  = y(1);
    f(0) = (U0(t) - Uc) / (R * C);
    f(1) = 1.0;
  }

  void evaluateDeriv (VectorView<double> y, MatrixView<double> df) const override
  {
    double t = y(1);
  
    df(0, 0) = -1.0 / (R * C);
    df(0, 1) = U0_prime(t) / (R * C);

    df(1, 0) = 0.0;
    df(1, 1) = 0.0;
  }
};

int main(int argc, char* argv[])
{
  double R = 1.0;
  double C = 1.0;

  // default values
  double tend = 0.1;
  int steps = 100;
  std::string method = "improved";

  if (argc > 1) steps = std::stoi(argv[1]); 
  if (argc > 2) method = argv[2];
  
  double tau = tend/steps;

  Vector<> y = { 0, 0 }; 
  auto rhs = std::make_shared<RCCircuit>(R, C);
  
  std::unique_ptr<ASC_ode::TimeStepper> stepper;

  if (method == "explicit_euler") {
    stepper = std::make_unique<ASC_ode::ExplicitEuler>(rhs);
  } 
  else if (method == "implicit_euler") {
    stepper = std::make_unique<ASC_ode::ImplicitEuler>(rhs);
  } 
  else if (method == "improved_euler") {
    stepper = std::make_unique<ASC_ode::ImprovedEuler>(rhs);
  } 
  else if (method == "crank_nicolson") {
    stepper = std::make_unique<ASC_ode::CrankNicolson>(rhs);
  }

  ofstream outfile("output_rc.txt");
  
  cout << "Simuliere RC-Kreis..." << endl;

  outfile << y(1) << " " << y(0) << " " << rhs->U0(y(1)) << endl;

  for (int i = 0; i < steps; i++)
  {
     stepper->DoStep(tau, y);
     outfile << y(1) << " " << y(0) << " " << rhs->U0(y(1)) << endl;
  }
  
  cout << "Fertig. Daten in output_rc.txt"<<endl;
}