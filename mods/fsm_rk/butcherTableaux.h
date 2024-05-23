struct butcherTableau{
    long unsigned int s; //number of stages
    std::string name;
    std::vector<double> b;
    std::vector<double> A;
};

#include <functional>
#include <map>

#define id(eps,sig) (eps-1)*(eps-2)/2+(sig-1)

void fillMethod(long unsigned int ns, butcherTableau &method, std::vector<double> b, std::vector<double> A, std::string name)
{
    method.s = ns;
    method.b.resize(ns);
    method.A.resize(ns*(ns-1)/2);

    method.name = name;

    method.b = b;
    method.A = A;
}

double c(uint64_t i, butcherTableau method)
{
  double c=0.0;
  for(uint64_t j=1;j<=i-1;++j) c+= method.A.at(id(i,j));
  return c;
}

butcherTableau paramRK2(double theta)
{
    long unsigned int ns=2;
    butcherTableau method;
    std::vector<double> b = {1-theta,theta};
    std::vector<double> A = {1.0/(2*theta)};
    fillMethod(ns,method,b,A,"default");
    return method;
}

butcherTableau heunRK2(void)
{
    butcherTableau method = paramRK2(0.5);
    method.name = "Heun RK2";
    return method;
}

butcherTableau stdRK2(void)
{
    butcherTableau method = paramRK2(1.0);
    method.name = "Standard RK2";
    return method;
}

butcherTableau paramEuler(void)
{
    long unsigned int ns=1;
    butcherTableau method;
    std::vector<double> b = {1.0};
    std::vector<double> A = {0.0};
    fillMethod(ns,method,b,A,"paramEuler");
    return method;
}

butcherTableau stdRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> A = {0.5,-1.0,2.0};
    std::vector<double> b = {1./6.0,2./3.0,1./6.0};
    fillMethod(ns,method,b,A,"Standard RK3");
    return method;
}

butcherTableau wrayRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> b = {0.25,0.0,0.75};
    std::vector<double> A = {8./15.0,0.25,5./12.0};
    fillMethod(ns,method,b,A, "Wray RK3");
    return method;
}

butcherTableau nystromRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> b = {0.25,3.0/8.0,3.0/8.0};
    std::vector<double> A = {2.0/3.0,0.0,2.0/3.0};
    fillMethod(ns,method,b,A,"Nystrom RK3");
    return method;
}

butcherTableau SSPRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> b = {1./6.0,1./6.0,2.0/3.0};
    std::vector<double> A = {1.0,0.25,0.25};
    fillMethod(ns,method,b,A,"Ralston RK3");
    return method;
}

butcherTableau ralstonRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> b = {2./9.0,1./3.0,4.0/9.0};
    std::vector<double> A = {0.5,0.0,2.0/3.0};
    fillMethod(ns,method,b,A,"Ralston RK3");
    return method;
}

butcherTableau kuntzmannRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> b = {0.2071768,0.3585646,0.4342585};
    std::vector<double> A = {0.5658162,-0.0581020,0.8256939};
    fillMethod(ns,method,b,A,"Kuntzmann RK3");
    return method;
}

butcherTableau heunRK3(void)
{
    long unsigned int ns=3;
    butcherTableau method;
    std::vector<double> b = {0.25,0.0,0.75};
    std::vector<double> A = {1./3.0,0,2./3.0};
    fillMethod(ns,method,b,A,"Heun RK3");
    return method;
}

butcherTableau genRK3(double x, double y)
{
    long unsigned int ns=3;
    butcherTableau method;
    double a21 = x;
    double a32 = y;
    double a31 = 0.5*(x-2*y-sqrt(pow(2*y-x,2)-4*(3*y*x*(x-1)+pow(y,2))));
    std::vector<double> A = {a21,a31,a32};
    double b3 = 1/(6*x*y);
    double b2 = (1/x)*(0.5-(a31+y)*b3);
    double b1 = 1-b2-b3;
    std::vector<double> b = {b1,b2,b3};
    fillMethod(ns,method,b,A, "Generic RK3");
    return method;
}

butcherTableau ps3p5q(double c3)
{
    long unsigned int ns=4;
    butcherTableau method;

    double a21 = (c3-1)/(4*c3-3);
    double a32 = (2*c3-1)/(2*a21);
    double a31 = c3-a32;
    double a42 = 1./3.0*((4*c3-3)/(a21*(2*c3-1))+0.5*pow(c3,2.0)/((c3-1)*(2*c3-1)));
    double a43 = a21/(2*c3-1);
    double a41 = 1-a42-a43;
    std::vector<double> A = {a21,a31,a32,a41,a42,a43};
    double b1 = 1/(12*(c3-1));
    double b2 = b1*pow(4*c3-3,2)/(2*c3-1);
    double b3 = -b1/(2*c3-1);
    double b4 = b1*(4*c3-3);
    std::vector<double> b = {b1,b2,b3,b4};
    fillMethod(ns,method,b,A,"ps3p5q("+std::to_string(c3)+")");
    return method;
}

butcherTableau stdRK4(void)
{
    long unsigned int ns=4;
    butcherTableau method;
    std::vector<double> A = {0.5,0.0,0.5,0.0,0.0,1.0};
    std::vector<double> b = {1./6.0,1./3.0,1./3.0,1./6.0};
    fillMethod(ns,method,b,A,"Standard RK4");
    return method;
}

butcherTableau varRK4(void)
{
    long unsigned int ns=4;
    butcherTableau method;
    std::vector<double> A = {1/3.0,-1/3.0,1.0,1.0,-1.0,1.0};
    std::vector<double> b = {1./8.,3./8.,3./8.,1./8.};
    fillMethod(ns,method,b,A,"Kutta RK4");
    return method;
}

butcherTableau gillRK4(void)
{
    long unsigned int ns=4;
    butcherTableau method;
    std::vector<double> A = {0.5,(sqrt(2)-1)/2.0,(2-sqrt(2))/2.0,0.0,-sqrt(2)/2.0,(1+sqrt(2))/2};
    std::vector<double> b = {1./6.,(2-sqrt(2))/6.,(2+sqrt(2))/6.,1./6.};
    fillMethod(ns,method,b,A,"Gill RK4");
    return method;
}

butcherTableau kuntzmannRK4(void)
{
    long unsigned int ns=4;
    butcherTableau method;
    std::vector<double> A = {88/220.0,-33/220.0,165/220.0,95/220.0,-75/220.0,200/220.0};
    std::vector<double> b = {55./360.,125./360.,125./360.,55./360.};
    fillMethod(ns,method,b,A,"Kuntzmann RK4");
    return method;
}

butcherTableau ps4p7q(void)
{
    long unsigned int ns=6;
    butcherTableau method;
    double a21 = 0.23593376536651968050;
    double a31 = 0.347507356584235168;
    double a32 = -0.135619353983464433;
    double a41 = -0.20592852403227;
    double a42 = 1.891790766221084;
    double a43 = -0.89775024478958;
    double a51 = -0.094354932814554;
    double a52 = 1.756171412237619;
    double a53 = -0.967078504769475;
    double a54 = 0.069328259979890148;
    double a61 = 0.14157883255197;
    double a62 = -1.17039696277833;
    double a63 = 1.30579112376331;
    double a64 = -2.203541368552894;
    double a65 = 2.9265683750159476;
    std::vector<double> A = {a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65};
    double b1 = 0.07078941627598264;
    double b2 = 0.87808570611880957;
    double b3 = -0.448875122394792210;
    double b4=b3;
    double b5=b2;
    double b6=b1;
    std::vector<double> b = {b1,b2,b3,b4,b5,b6};
    fillMethod(ns,method,b,A,"ps4p7q");
    return method;
}

using FunctionTypeVoid = std::function<butcherTableau()>;
using FunctionTypeDouble = std::function<butcherTableau(double)>;
using FunctionTypeDoubleDouble = std::function<butcherTableau(double,double)>;

std::map<std::string, FunctionTypeVoid> functionMapVoid = {
  {"heunRK2",heunRK2},
  {"stdRK2",stdRK2},
  {"paramEuler",paramEuler},
  {"stdRK3",stdRK3},
  {"wrayRK3",wrayRK3},
  {"nystromRK3",nystromRK3},
  {"ralstonRK3",ralstonRK3},
  {"SSPRK3",SSPRK3},
  {"kuntzmannRK3",kuntzmannRK3},
  {"heunRK3",heunRK3},
  {"stdRK4",stdRK4},
  {"varRK4",varRK4},
  {"gillRK4",gillRK4},
  {"kuntzmannRK4",kuntzmannRK4},
  {"ps4p7q",ps4p7q}
};

std::map<std::string, FunctionTypeDouble> functionMapDouble = {
  {"paramRK2",paramRK2},
  {"ps3p5q",ps3p5q}
};

std::map<std::string, FunctionTypeDoubleDouble> functionMapDoubleDouble = {
  {"genRK3",genRK3}
};

butcherTableau intScheme(std::string functionName, double param1 = 0.0, double param2 = 0.0)
{
  auto itVoid = functionMapVoid.find(functionName);
  if (itVoid != functionMapVoid.end()){
    FunctionTypeVoid function = itVoid->second;
    return function();
  }

  auto itDouble = functionMapDouble.find(functionName);
  if (itDouble != functionMapDouble.end()){
    FunctionTypeDouble function = itDouble->second;
    return function(param1);
  }

  auto itDoubleDouble = functionMapDoubleDouble.find(functionName);
  if (itDoubleDouble != functionMapDoubleDouble.end()){
    FunctionTypeDoubleDouble function = itDoubleDouble->second;
    return function(param1,param2);
  }

  FunctionTypeDoubleDouble function = itDoubleDouble->second;
  return function(param1,param2);
}

void displayMethod(butcherTableau method)
{
  std::cout << "Scheme used: " << method.name << std::endl;
  std::cout << "\tNumber of stages: " << method.s << std::endl;
  std::cout << "\tButcher's tableau:" << std::endl;
  std::cout << "\t\t b: ( ";
  for (const auto &element: method.b)
    std::cout << element << " ";
  std::cout << ")" << std::endl;
  std:: cout << "\t\t A: ( ";
  for (const auto &element :method.A)
    std::cout << element << " ";
  std::cout << ")" << std::endl;
}
