// $Id$
//==============================================================================
//!
//! \file PoroMaterial.h
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poro-elastic material models.
//!
//==============================================================================

#ifndef _PORO_MATERIAL_H
#define _PORO_MATERIAL_H

#include "Function.h"
#include "MatVec.h"
#include "Vec3.h"
#include "Vec3Oper.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a poroelastic problem.
*/

class PoroMaterial
{
public:
  //! \brief Empty constructor.
  PoroMaterial() {}

  //! \brief Empty destructor.
  ~PoroMaterial() {}

  //! \brief Parses material parementers from an XML element.
  void parse(const TiXmlElement*);

  //! \brief Prints out material parameters to the log stream.
  void printLog() const;

  //! \brief Returns mean mass density
  double getMassDensity(const Vec3&) const;
  //! \brief Evaluates the mass density of the fluid at current point.
  double getFluidDensity(const Vec3&) const;
  //! \brief Evaluates the mass density of the solid at current point.
  double getSolidDensity(const Vec3&) const;
  //! \brief Evaluates the effective heat capacity at the current point
  double getHeatCapacity(double T) const;
  //! \brief Evaluates the heat capacity at the current point
  double getFluidHeatCapacity(double T) const;
  //! \brief Evaluates the heat capacity at the current point
  double getSolidHeatCapacity(double T) const;
  //! \brief Evaluates the effective thermal conductivity at the current point
  double getThermalConductivity(double T) const;
  //! \brief Evaluates the thermal conductivity of the fluid at the current point
  double getFluidThermalConductivity(double T) const;
  //! \brief Evaluates the thermal conductivity of the solid at the current point
  double getSolidThermalConductivity(double T) const;
  //! \brief Evaluates the thermal expansion of the solid at the current point
  double getSolidThermalExpansion(double T) const;
  //! \brief Returns porosity at the current point.
  double getPorosity(const Vec3& X) const;
  //! \brief Returns permeability at the current point.
  Vec3 getPermeability(const Vec3& X) const;
  //! \brief Returns bulk modulus of the fluid at the current point.
  double getBulkFluid(const Vec3& X) const;
  //! \brief Returns bulk modulus of the solid at the current point.
  double getBulkSolid(const Vec3& X) const;
  //! \brief Returns bulk modulus of the medium at the current point.
  double getBulkMedium(const Vec3& X) const;
  //! \brief Returns stiffness at the current point.
  double getStiffness(const Vec3& X) const;
  //! \brief Returns Poisson's ratio at the current point.
  double getPoisson(const Vec3& X) const;

  //! \brief Calculates the strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] dNdX First basis function gradients at current point
  bool formBmatrix(Matrix& Bmat, const Matrix& dNdX, size_t nsd) const;

  //! \brief Evalutates the constitutive matrix at an integration point
  //! \param[out] Cmat Constitutive matrix at current point
  //! \param[in] E Young's modulus
  //! \param[in] nu Poisson's ratio
  bool formElasticMatrix(Matrix& Cmat, const Vec3& X, size_t nsd) const;

protected:
  /*! \brief Helper template for wrapping a constant/function pair */
    template<class Function>
  struct FuncConstPair
  {
    Function* function;                 //!< Function definition
    typename Function::Output constant; //!< Constant

    //! \brief Constructor.
    FuncConstPair() { function = nullptr; constant = 0.0; }

    //! \brief Parse an XML element. Specialized per type
    Function* parse(const char* val, const std::string& type) { return nullptr; }

    //! \brief Evaluate function.
    //! \param[in] X the value to evaluate at.
    typename Function::Output evaluate(const typename Function::Input& X) const
    {
      return function ? (*function)(X) : constant;
    }
  };

  FuncConstPair<RealFunc> Emod; //!< Young's modulus
  FuncConstPair<RealFunc> nu;   //!< Poisson's ratio
  FuncConstPair<RealFunc> rhof; //!< Fluid density
  FuncConstPair<RealFunc> rhos; //!< Solid density

  FuncConstPair<ScalarFunc> fheatcapacity; //!< Specific heat capacity for fluid
  FuncConstPair<ScalarFunc> sheatcapacity; //!< Specific heat capacity for solid
  FuncConstPair<ScalarFunc> fconductivity; //!< Thermal conductivity
  FuncConstPair<ScalarFunc> sconductivity; //!< Thermal conductivity
  FuncConstPair<ScalarFunc> sexpansion;    //!< Thermal expansion

  FuncConstPair<RealFunc> porosity;     //!< Porosity
  FuncConstPair<VecFunc>  permeability; //!< Permeability
  FuncConstPair<RealFunc> bulkf;        //!< Bulk modulus of fluid
  FuncConstPair<RealFunc> bulks;        //!< Bulk modulus of solid
  FuncConstPair<RealFunc> bulkm;        //!< Bulk modulus of medium
};

#endif
