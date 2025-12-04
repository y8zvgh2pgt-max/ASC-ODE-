#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void doStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void doStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void doStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };


  class ImprovedEuler : public TimeStepper
  {
    Vector<> m_k1;
    Vector<> m_k2;
    Vector<> m_ytmp;

  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), 
      m_k1(rhs->dimF()), 
      m_k2(rhs->dimF()), 
      m_ytmp(rhs->dimX()) 
    {}

    void doStep(double tau, VectorView<double> y) override
    {
      m_rhs->evaluate(y, m_k1);

      m_ytmp = y;
      m_ytmp += 0.5 * tau * m_k1;

      m_rhs->evaluate(m_ytmp, m_k2);
      y += tau * m_k2;
    }

  };


  class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<ConstantFunction> m_f_old;
    Vector<> m_tmp_f;

  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs),
      m_tau(std::make_shared<Parameter>(0.0)),
      m_yold(std::make_shared<ConstantFunction>(rhs->dimX())),
      m_f_old(std::make_shared<ConstantFunction>(rhs->dimF())),
      m_tmp_f(m_rhs->dimF())
    {      
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      auto half_tau_f_new = 0.5 * (m_tau * m_rhs);
      auto half_tau_f_old = 0.5 * (m_tau * m_f_old);
      
      m_equ = ynew - m_yold - half_tau_f_old - half_tau_f_new;
    }

    void doStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      
      m_rhs->evaluate(y, m_tmp_f);
      m_f_old->set(m_tmp_f);

      NewtonSolver(m_equ,y);
    }
  };
}

#endif
