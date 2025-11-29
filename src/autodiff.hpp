#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <cstddef> 
#include <ostream> 
#include <cmath>   
#include <array>  


namespace ASC_ode
{

  template <size_t N, typename T = double>
  class Variable 
  {
    private:
      T m_val;
    public:
      Variable (T v) : m_val(v) {}
      T value() const { return m_val; }
  };

  template <typename T = double>
  auto derivative (T v, size_t /*index*/) { return T(0); } 


  template <size_t N, typename T = double>
  class AutoDiff
  {
  private:
    T m_val;
    std::array<T, N> m_deriv;
  public: 
    AutoDiff () : m_val(0), m_deriv{} {}
    AutoDiff (T v) : m_val(v), m_deriv{} 
    {
      for (size_t i = 0; i < N; i++)
        m_deriv[i] = derivative(v, i);
    }
    
    template <size_t I>
    AutoDiff (Variable<I, T> var) : m_val(var.value()), m_deriv{} 
    {
      m_deriv[I] = 1.0;
    }

    T value() const { return m_val; }
    std::array<T, N>& deriv() { return m_deriv; }
    const std::array<T, N>& deriv() const { return m_deriv; }
  };


  template <size_t N, typename T = double>
  auto derivative (AutoDiff<N, T> v, size_t index) 
  {
    return v.deriv()[index];
  }



  template <size_t N, typename T>
  std::ostream & operator<< (std::ostream& os, const AutoDiff<N, T>& ad)
  {
    os << "Value: " << ad.value() << ", Deriv: [";
    for (size_t i = 0; i < N; i++)
    {
      os << ad.deriv()[i];
      if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
  }

  template <size_t N, typename T = double>
  AutoDiff<N, T> operator+ (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
  {
     AutoDiff<N, T> result(a.value() + b.value());
     for (size_t i = 0; i < N; i++)
        result.deriv()[i] = a.deriv()[i] + b.deriv()[i];
       return result;
   }

   template <size_t N, typename T = double>
   auto operator+ (T a, const AutoDiff<N, T>& b) { return AutoDiff<N, T>(a) + b; }


   template <size_t N, typename T = double>
   AutoDiff<N, T> operator* (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
   {
       AutoDiff<N, T> result(a.value() * b.value());
       for (size_t i = 0; i < N; i++)
          result.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
       return result;
   }

   using std::sin;
   using std::cos;

   template <size_t N, typename T = double>
   AutoDiff<N, T> sin(const AutoDiff<N, T> &a)
   {
       AutoDiff<N, T> result(sin(a.value()));
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = cos(a.value()) * a.deriv()[i];
       return result;
   }

   // --- Ergänzungen für autodiff.hpp ---

   // Subtraktion (a - b)
   template <size_t N, typename T = double>
   AutoDiff<N, T> operator- (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
   {
       AutoDiff<N, T> result(a.value() - b.value());
       for (size_t i = 0; i < N; i++)
          result.deriv()[i] = a.deriv()[i] - b.deriv()[i];
       return result;
   }

   // Unäres Minus (-a)
   template <size_t N, typename T = double>
   AutoDiff<N, T> operator- (const AutoDiff<N, T>& a)
   {
       AutoDiff<N, T> result(-a.value());
       for (size_t i = 0; i < N; i++)
          result.deriv()[i] = -a.deriv()[i];
       return result;
   }

   // Division (a / b) - Quotientenregel!
   template <size_t N, typename T = double>
   AutoDiff<N, T> operator/ (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
   {
       T val = a.value() / b.value();
       AutoDiff<N, T> result(val);
       T b_sq = b.value() * b.value();
       for (size_t i = 0; i < N; i++)
          // (u'v - uv') / v^2
          result.deriv()[i] = (a.deriv()[i] * b.value() - a.value() * b.deriv()[i]) / b_sq;
       return result;
   }

   // Cosinus (Kettenregel: cos(u)' = -sin(u) * u')
   template <size_t N, typename T = double>
   AutoDiff<N, T> cos(const AutoDiff<N, T> &a)
   {
       AutoDiff<N, T> result(cos(a.value()));
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = -sin(a.value()) * a.deriv()[i];
       return result;
   }

   // Exponentialfunktion (exp(u)' = exp(u) * u')
   template <size_t N, typename T = double>
   AutoDiff<N, T> exp(const AutoDiff<N, T> &a)
   {
       T exp_val = exp(a.value());
       AutoDiff<N, T> result(exp_val);
       for (size_t i = 0; i < N; i++)
           result.deriv()[i] = exp_val * a.deriv()[i];
       return result;
   }
   
   // Misch-Operatoren (AutoDiff +/-/* double) fehlen oft noch, 
   // aber für das Pendel reicht meistens das oben.


} // namespace ASC_ode

#endif
