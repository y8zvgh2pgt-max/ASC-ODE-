#include <iostream>
#include <fstream> 
#include <string>
#include <memory>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>

using namespace ASC_ode;

class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};


int main(int argc, char* argv[])
{
  // default values
  double tend = 4*M_PI;
  int steps = 100;
  std::string method = "improved";
  int stages = 3;

  if (argc > 1) steps = std::stoi(argv[1]); 
  if (argc > 2) tend = std::stod(argv[2]);
  if (argc > 3) method = argv[3];
  if (argc > 4) stages = std::stoi(argv[4]);   // only used for gauss_legendre / radau

  double tau = tend/steps;

  Vector<> y = { 1, 0 }; 
  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

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
  // fixed Gauss-Legendre 2-stage (order 4)
  else if (method == "implicit_rk_gauss2") {
    stepper = std::make_unique<ASC_ode::ImplicitRungeKutta>(rhs, Gauss2a, Gauss2b, Gauss2c);
  }
  // fixed Gauss-Legendre 3-stage (order 6)
  else if (method == "implicit_rk_gauss3") {
    auto [Gauss3a, Gauss3b] = computeABfromC(Gauss3c);
    stepper = std::make_unique<ASC_ode::ImplicitRungeKutta>(rhs, Gauss3a, Gauss3b, Gauss3c);
  }
  // arbitrary-order Gauss-Legendre (implicit RK, order 2*stages)
  else if (method == "implicit_rk_gauss_legendre") {
    Vector<> c(stages), w(stages);
    GaussLegendre(c, w);
    auto [a, b] = computeABfromC(c);
    stepper = std::make_unique<ASC_ode::ImplicitRungeKutta>(rhs, a, b, c);
  }
  // arbitrary-order Radau IIA (implicit RK, order 2*stages-1)
  else if (method == "implicit_rk_radau") {
    Vector<> c(stages), w(stages);
    GaussRadau(c, w);
    auto [a, b] = computeABfromC(c);
    stepper = std::make_unique<ASC_ode::ImplicitRungeKutta>(rhs, a, b, c);
  }
  // explicit RK methods via ExplicitRungeKutta
  else if (method == "explicit_rk2") {          // explicit midpoint (2-stage RK2)
    stepper = std::make_unique<ExplicitRungeKutta>(rhs, ERK2A, ERK2b, ERK2c);
  }
  else if (method == "explicit_rk4") {          // classical RK4
    stepper = std::make_unique<ExplicitRungeKutta>(rhs, ERK4A, ERK4b, ERK4c);
  }

  std::ofstream outfile ("output_test_ode.txt");
  outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     stepper->doStep(tau, y);
     outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
}