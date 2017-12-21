/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::heTabularPsiThermo

Description
    Enthalpy for a mixture based on compressibility

SourceFiles
    heTabularPsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heTabularPsiThermo_H
#define heTabularPsiThermo_H

#include "basicPsiThermo.H"
#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heTabularPsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class heTabularPsiThermo
:
    public basicPsiThermo,
    public MixtureType
{
    // Private data

        //- Enthalpy field
        volScalarField h_;

        //- Internal Energy field
        volScalarField e_;


    // Private member functions

        //- Calculate the thermo variables
        void calculate();
        //- Calculate the thermo variables with rho and e: Janry
        void calculate2();

        //- Construct as copy (not implemented)
        heTabularPsiThermo(const heTabularPsiThermo<MixtureType>&);


public:

    //- Runtime type information
    TypeName("heTabularPsiThermo");


    // Constructors

        //- Construct from mesh
        heTabularPsiThermo(const fvMesh&, const objectRegistry& obj);


    //- Destructor
    virtual ~heTabularPsiThermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual basicMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const basicMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }



            //- Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& e()
            {
                return e_;
            }

            //- Internal Energy [J/kg]
            virtual const volScalarField& e() const
            {
                return e_;
            }




        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;


            //- Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Internal energy for patch [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;



            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
             const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for cell-set [J/kg/K]
            virtual tmp<scalarField> Cp
            (
             const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
             const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "heTabularPsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
