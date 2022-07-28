// $Id$
//==============================================================================
//!
//! \file PoroMaterial.C
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poro-elastic material models.
//!
//==============================================================================


#include "PoroMaterial.h"

#include "Functions.h"
#include "IFEM.h"
#include "tinyxml.h"
#include "Utilities.h"
#include "Vec3.h"


template<>
RealFunc* PoroMaterial::FuncConstPair<RealFunc>::parse(const char* val,
                                                       const std::string& type)
{
  return utl::parseRealFunc(val, type);
}


template<>
ScalarFunc* PoroMaterial::FuncConstPair<ScalarFunc>::parse(const char* val,
                                                           const std::string& type)
{
  return utl::parseTimeFunc(val, type);
}


template<>
VecFunc* PoroMaterial::FuncConstPair<VecFunc>::parse(const char* val,
                                                     const std::string& type)
{
  return utl::parseVecFunc(val, type);
}


  template<class T>
static bool propertyParse(PoroMaterial::FuncConstPair<T>& data,
                          const TiXmlElement* elem,
                          const std::string& attr,
                          const std::string& tag)
{
  std::string constant;
  if (utl::getAttribute(elem,attr.c_str(),constant)) {
    std::stringstream str;
    str << constant;
    str >> data.constant;
    return true;
  }

  const TiXmlElement* child = elem->FirstChildElement(tag);
  if (child) {
    IFEM::cout <<" ";
    std::string type;
    utl::getAttribute(child,"type",type,true);
    const TiXmlNode* aval;
    if ((aval = child->FirstChild()))
      data.function = data.parse(aval->Value(),type);

    return data.function != nullptr;
  }

 return false;
}


void PoroMaterial::parse(const TiXmlElement* elem)
{
  propertyParse(Emod, elem, "E", "stiffness");
  propertyParse(nu, elem, "nu", "poisson");
  propertyParse(rhof, elem, "rhof", "fluiddensity");
  propertyParse(rhos, elem, "rhos", "soliddensity");

  propertyParse(fheatcapacity, elem, "cpf", "fluidheatcapacity");
  propertyParse(sheatcapacity, elem, "cps", "solidheatcapacity");
  propertyParse(fconductivity, elem, "kappaf", "fluidconductivity");
  propertyParse(sconductivity, elem, "kappas", "solidconductivity");
  propertyParse(sexpansion, elem, "alphas", "solidexpansion");

  propertyParse(porosity, elem, "poro", "porosity");
  propertyParse(permeability, elem, "perm", "permeability");
  propertyParse(bulkf, elem, "Kf", "fluidbulk");
  propertyParse(bulks, elem, "Ks", "solidbulk");
  propertyParse(bulkm, elem, "Ko", "mediumbulk");
}


void PoroMaterial::printLog() const
{
  IFEM::cout << "\tConstitutive Properties: "
             << "\n\t\tYoung's Modulus, E = " << Emod.constant
             << "\n\t\tPoisson's Ratio, nu = " << nu.constant << std::endl;
  IFEM::cout << "\tDensities: "
             << "\n\t\tDensity of Fluid, rhof = " << rhof.constant
             << "\n\t\tDensity of Solid, rhos = " << rhos.constant << std::endl;
  IFEM::cout << "\tBulk Moduli: "
             << "\n\t\tBulk Modulus of Fluid, Kf = " << bulkf.constant
             << "\n\t\tBulk Modulus of Solid, Ks = " << bulks.constant
             << "\n\t\tBulk Modulus of Medium, Ko = " << bulkm.constant << std::endl;
  IFEM::cout <<"\tPorosity, n = " << porosity.constant << std::endl;
  IFEM::cout << "\tPermeability, k = " << permeability.constant << std::endl;
  IFEM::cout << "\tHeat capacities: "
             << "\n\t\tHeat capacity of Fluid, cpf= " << fheatcapacity.constant
             << "\n\t\tHeat capacity of Solid, cps= " << sheatcapacity.constant << std::endl;
  IFEM::cout << "\tThermal Conductivities: "
             << "\n\t\tThermal Conductivity of Fluid, kappaf= " << fconductivity.constant
             << "\n\t\tThermal Conductivity of Solid, kappas= " << sconductivity.constant
             << std::endl;
}


double PoroMaterial::getMassDensity(const Vec3& X) const
{
  return getSolidDensity(X)*(1.0-getPorosity(X)) +
         getFluidDensity(X)*getPorosity(X);
}


double PoroMaterial::getHeatCapacity(double T) const
{
  Vec3 X;
  return getPorosity(X)*getFluidDensity(X)*getFluidHeatCapacity(T) +
         (1.0-getPorosity(X))*getSolidDensity(X)*getSolidHeatCapacity(T);
}


double PoroMaterial::getFluidHeatCapacity (double T) const
{
  return fheatcapacity.evaluate(T);
}


double PoroMaterial::getSolidHeatCapacity (double T) const
{
  return sheatcapacity.evaluate(T);
}


double PoroMaterial::getFluidThermalConductivity(double T) const
{
  return fconductivity.evaluate(T);
}


double PoroMaterial::getSolidThermalConductivity(double T) const
{
  return sconductivity.evaluate(T);
}


double PoroMaterial::getThermalConductivity(double T) const
{
  Vec3 X;
  return pow(getFluidThermalConductivity(T),getPorosity(X))*
         pow(getSolidThermalConductivity(T),1.0-getPorosity(X));
}


double PoroMaterial::getSolidThermalExpansion(double T) const
{
  return sexpansion.evaluate(T);
}


double PoroMaterial::getPorosity(const Vec3& X) const
{
  return porosity.evaluate(X);
}


Vec3 PoroMaterial::getPermeability(const Vec3& X) const
{
  return permeability.evaluate(X);
}


double PoroMaterial::getFluidDensity(const Vec3& X) const
{
  return rhof.evaluate(X);
}


double PoroMaterial::getSolidDensity(const Vec3& X) const
{
  return rhos.evaluate(X);
}


double PoroMaterial::getBulkFluid(const Vec3& X) const
{
  return bulkf.evaluate(X);
}


double PoroMaterial::getBulkSolid(const Vec3& X) const
{
  return bulks.evaluate(X);
}


double PoroMaterial::getBulkMedium(const Vec3& X) const
{
  return bulkm.evaluate(X);
}


double PoroMaterial::getStiffness(const Vec3& X) const
{
  return Emod.evaluate(X);
}


double PoroMaterial::getPoisson(const Vec3& X) const
{
  return nu.evaluate(X);
}


bool PoroMaterial::formBmatrix (Matrix& Bmat, const Matrix& dNdX, size_t nsd) const
{
  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);
  if (dNdX.cols() < nsd)
  {
    std::cerr << " *** PoroMaterial::formBmatrix: Invalid dimension on dN1dX, "
              << dNdX.rows() << "x" << dNdX.cols() << "." << std::endl;
    return false;
  }

#define INDEX(i,j) i+nstrc*(j-1)

  // Strain-displacement matrix for 2D elements
  for (size_t i = 1; i <= nenod; i++)
  {
    Bmat(INDEX(1,1),i) = dNdX(i,1);
    Bmat(INDEX(2,2),i) = dNdX(i,2);
    Bmat(INDEX(3,1),i) = dNdX(i,2);
    Bmat(INDEX(3,2),i) = dNdX(i,1);
  }

#undef INDEX

  Bmat.resize(nstrc,nsd*nenod);

  return true;
}


bool PoroMaterial::formElasticMatrix(Matrix& Cmat, const Vec3& X, size_t nsd) const
{
  double E = getStiffness(X);
  double nu =getPoisson(X);

  double C33 = E/(2.0+2.0*nu);
  double C12 = (E*nu)/((1.0+nu)*(1.0-2.0*nu));
  double C11 = C12 + 2.0*C33;

  const size_t nstrc = nsd*(nsd+1)/2;

  Cmat.resize(nstrc,nstrc,true);

  Cmat(1,1) = C11;
  Cmat(1,2) = C12;
  Cmat(2,1) = C12;
  Cmat(2,2) = C11;
  Cmat(3,3) = C33;

  Cmat.resize(nstrc,nstrc);

  return true;
}
