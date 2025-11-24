#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

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

  if (argc > 1) {
      steps = std::stoi(argv[1]); 
  }

  if (argc > 2) {
    tend = std::stod(argv[2]);
  }
  
  if (argc > 3) {
    method = argv[3];
  }



  
  double tau = tend/steps;

  Vector<> y = { 1, 0 };  // initializer list
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



  std::ofstream outfile ("output_test_ode.txt");
  std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
  outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     stepper->DoStep(tau, y);

     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
     outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }

}
