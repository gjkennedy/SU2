/*!
 * \file CIncEulerVariable.hpp
 * \brief Class for defining the variables of the incompressible Euler solver.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CVariable.hpp"

/*!
 * \class CIncEulerVariable
 * \brief Class for defining the variables of the incompressible Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncEulerVariable : public CVariable {
public:
  static constexpr size_t MAXNVAR = 12;

protected:
  VectorType Velocity2;                     /*!< \brief Square of the velocity vector. */
  MatrixType Primitive;                     /*!< \brief Primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) in incompressible flows. */
  CVectorOfMatrix Gradient_Primitive;       /*!< \brief Gradient of the primitive variables (P, vx, vy, vz, T, rho, beta). */
  CVectorOfMatrix& Gradient_Reconstruction; /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  CVectorOfMatrix Gradient_Aux;             /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */
  MatrixType Limiter_Primitive;             /*!< \brief Limiter of the primitive variables (P, vx, vy, vz, T, rho, beta). */

  /*--- NS Variables declared here to make it easier to re-use code between compressible and incompressible solvers. ---*/
  MatrixType Vorticity;       /*!< \brief Vorticity of the fluid. */
  VectorType StrainMag;       /*!< \brief Magnitude of rate of strain tensor. */

  VectorType Streamwise_Periodic_RecoveredPressure,    /*!< \brief Recovered/Physical pressure [Pa] for streamwise periodic flow. */
             Streamwise_Periodic_RecoveredTemperature; /*!< \brief Recovered/Physical temperature [K] for streamwise periodic flow. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_pressure - value of the pressure.
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] temperature - Value of the temperature (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncEulerVariable(su2double pressure, const su2double *velocity, su2double temperature,
                    unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CIncEulerVariable() override = default;

  /*!
   * \brief Get the primitive variable gradients for all points.
   * \return Reference to primitive variable gradient.
   */
  inline CVectorOfMatrix& GetGradient_Primitive() final { return Gradient_Primitive; }
  inline const CVectorOfMatrix& GetGradient_Primitive() const final { return Gradient_Primitive; }

  /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction() final { return Gradient_Reconstruction; }
  inline const CVectorOfMatrix& GetGradient_Reconstruction() const final { return Gradient_Reconstruction; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
    return Gradient_Primitive(iPoint,iVar,iDim);
  }

  /*!
   * \brief Get the primitive variables limiter.
   * \return Primitive variables limiter for the entire domain.
   */
  inline MatrixType& GetLimiter_Primitive() final {return Limiter_Primitive; }
  inline const MatrixType& GetLimiter_Primitive() const final {return Limiter_Primitive; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetLimiter_Primitive(unsigned long iPoint, unsigned long iVar) const final {
    return Limiter_Primitive(iPoint,iVar);
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the primitive variables gradient.
   */
  inline CMatrixView<su2double> GetGradient_Primitive(unsigned long iPoint, unsigned long iVar=0) final {
    return Gradient_Primitive(iPoint,iVar);
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the primitive variables gradient.
   */
  inline su2double *GetLimiter_Primitive(unsigned long iPoint) final { return Limiter_Primitive[iPoint]; }

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) final { return Gradient_Reconstruction[iPoint]; }

  /*!
   * \brief Set the value of the pressure.
   * \param[in] iPoint - Point index.
   */
  inline void SetPressure(unsigned long iPoint) final { Primitive(iPoint,0) = Solution(iPoint,0); }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive() const final { return Primitive; }

  /*!
   * \brief Get the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
  inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final { return Primitive(iPoint,iVar); }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) final { Primitive(iPoint,iVar) = val_prim; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, const su2double *val_prim) final {
    for (unsigned long iVar = 0; iVar < nPrimVar; iVar++) Primitive(iPoint,iVar) = val_prim[iVar];
  }

  /*!
   * \brief Get the primitive variables of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitive(unsigned long iPoint) final { return Primitive[iPoint]; }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline bool SetDensity(unsigned long iPoint, su2double val_density) final {
    Primitive(iPoint,nDim+2) = val_density;
    if (Primitive(iPoint,nDim+2) > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline void SetVelocity(unsigned long iPoint) final {
    Velocity2(iPoint) = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      Primitive(iPoint,iDim+1) = Solution(iPoint,iDim+1);
      Velocity2(iPoint) += pow(Primitive(iPoint,iDim+1),2);
    }
  }

  /*!
   * \brief Set the value of the temperature for incompressible flows with energy equation.
   * \param[in] iPoint - Point index.
   */
  inline bool SetTemperature(unsigned long iPoint, su2double val_temperature) final {
    Primitive(iPoint,nDim+1) = val_temperature;
    if (Primitive(iPoint,nDim+1) > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the beta coeffient for incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline void SetBetaInc2(unsigned long iPoint, su2double val_betainc2) final { Primitive(iPoint,nDim+3) = val_betainc2; }

  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  inline su2double GetVelocity2(unsigned long iPoint) const final { return Velocity2(iPoint); }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(unsigned long iPoint) const final { return Primitive(iPoint,0); }

  /*!
   * \brief Get the value of beta squared for the incompressible flow
   * \return Value of beta squared.
   */
  inline su2double GetBetaInc2(unsigned long iPoint) const final { return Primitive(iPoint,nDim+3); }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(unsigned long iPoint) const final { return Primitive(iPoint,nDim+2); }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(unsigned long iPoint) const final { return Primitive(iPoint,nDim+1); }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final { return Primitive(iPoint,iDim+1); }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const final {
    su2double ProjVel = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      ProjVel += Primitive(iPoint,iDim+1)*val_vector[iDim];
    return ProjVel;
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Old(iPoint,iDim+1) = val_velocity[iDim];
  }

  /*!
   * \brief Set the momentum part of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVel_ResTruncError_Zero(unsigned long iPoint) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Res_TruncError(iPoint,iDim+1) = 0.0;
  }

  /*!
   * \brief Set all the primitive variables for incompressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) final;

  /*!
   * \brief Set the specific heat Cp.
   */
  inline void SetSpecificHeatCp(unsigned long iPoint, su2double val_Cp) final { Primitive(iPoint, nDim+7) = val_Cp; }

  /*!
   * \brief Set the specific heat Cv.
   */
  inline void SetSpecificHeatCv(unsigned long iPoint, su2double val_Cv) final { Primitive(iPoint, nDim+8) = val_Cv; }

  /*!
   * \brief Get the specific heat at constant P of the flow.
   * \return Value of the specific heat at constant P of the flow.
   */
  inline su2double GetSpecificHeatCp(unsigned long iPoint) const final { return Primitive(iPoint, nDim+7); }

  /*!
   * \brief Get the specific heat at constant V of the flow.
   * \return Value of the specific heat at constant V of the flow.
   */
  inline su2double GetSpecificHeatCv(unsigned long iPoint) const final { return Primitive(iPoint, nDim+8); }

  /*!
   * \brief Get the value of the vorticity.
   * \return Value of the vorticity.
   */
  inline su2double *GetVorticity(unsigned long iPoint) final { return Vorticity[iPoint]; }

  /*!
   * \brief Get the value of the magnitude of rate of strain.
   * \return Value of the rate of strain magnitude.
   */
  inline su2double GetStrainMag(unsigned long iPoint) const final { return StrainMag(iPoint); }
  inline su2activevector& GetStrainMag() { return StrainMag; }

  /*!
   * \brief Set the recovered pressure for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \param[in] val_pressure - pressure value.
   */
  inline void SetStreamwise_Periodic_RecoveredPressure(unsigned long iPoint, su2double val_pressure) final {
    Streamwise_Periodic_RecoveredPressure(iPoint) = val_pressure;
  }

  /*!
   * \brief Get the recovered pressure for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \return Recovered/Physical pressure for streamwise periodic flow.
   */
  inline su2double GetStreamwise_Periodic_RecoveredPressure(unsigned long iPoint) const final {
    return Streamwise_Periodic_RecoveredPressure(iPoint);
  }

  /*!
   * \brief Set the recovered temperature for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \param[in] val_temperature - temperature value.
   */
  inline void SetStreamwise_Periodic_RecoveredTemperature(unsigned long iPoint, su2double val_temperature) final {
    Streamwise_Periodic_RecoveredTemperature(iPoint) = val_temperature;
  }

  /*!
   * \brief Get the recovered temperature for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \return Recovered/Physical temperature for streamwise periodic flow.
   */
  inline su2double GetStreamwise_Periodic_RecoveredTemperature(unsigned long iPoint) const final {
    return Streamwise_Periodic_RecoveredTemperature(iPoint);
  }

  /*!
   * \brief Specify a vector to set the velocity components of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionVector(unsigned long iPoint, const su2double *val_vector) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Solution(iPoint, iDim+1) = val_vector[iDim];
  }

};
