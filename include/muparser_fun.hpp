#include <muParser.h>
#include <memory>
#include <string>

#ifndef MUPARSER_FUN_HPP
#define MUPARSER_FUN_HPP
class MuparserFun
{
public:
  MuparserFun(const MuparserFun &m)
    : m_parser(m.m_parser)
  {
    m_parser.DefineVar("x", &m_varx);
    m_parser.DefineVar("y", &m_vary);
  };

  MuparserFun(const std::string &s)
  {
    try
      {
        m_parser.DefineVar("x", &m_varx);
        m_parser.DefineVar("y", &m_vary);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x, const double &y)
  {
    m_varx = x;
    m_vary = y;
    double f = m_parser.Eval();
    return f;
  };

private:
  double     m_varx;
  double     m_vary;
  mu::Parser m_parser;
};

#endif