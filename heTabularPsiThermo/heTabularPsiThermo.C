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

\*---------------------------------------------------------------------------*/

#include "heTabularPsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::heTabularPsiThermo<MixtureType>::calculate()
{
    const scalarField& hCells = h_.internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.TH(hCells[celli], pCells[celli], TCells[celli]);
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alpha(pCells[celli], TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = h_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.H(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alpha(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.TH(ph[facei], pp[facei], pT[facei]);
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alpha(pp[facei], pT[facei]);
            }
        }
    }
}

template<class MixtureType>
void Foam::heTabularPsiThermo<MixtureType>::calculate2()
{
    const scalarField& eCells = e_.internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    


    //scalarField  rhoCells = h_.mesh().lookupObject<volScalarField>("rho").internalField();//add by Janry
    //scalarField  eCells   = h_.mesh().lookupObject<volScalarField>("e").internalField();  //add by Janry

    //Info << "rhoCells[4]=" << rhoCells[4] << " eCells[4] = " << eCells[4] << endl;
    //Info << "rhoCells[7]=" << rhoCells[7] << " eCells[7] = " << eCells[7] << endl;
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        //TCells[celli] = mixture_.TH(hCells[celli], pCells[celli], TCells[celli]);
        TCells[celli] = mixture_.TEP(eCells[celli], pCells[celli],TCells[celli]);// add by Janry
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alpha(pCells[celli], TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = h_.boundaryField()[patchi];
        fvPatchScalarField& pe = e_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        //fvPatchScalarField  prho = h_.mesh().lookupObject<volScalarField>("rho").boundaryField()[patchi];//add by Janry
        //fvPatchScalarField  pe   = h_.mesh().lookupObject<volScalarField>("e").boundaryField()[patchi];  //add by Janry

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.H(pp[facei], pT[facei]);
                pe[facei] = mixture_.E(pp[facei], pT[facei]);// should I add this one?
                //Info << "ph[" << facei << "]=" << ph[facei] <<  endl;
                //Info << "pe[" << facei << "]=" << pe[facei] <<  endl;
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alpha(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.TEP(pe[facei], pp[facei], pT[facei]);// Add by Janry
                //Info << "T["<<facei<<"] =" << pT[facei] << " e = " << pe[facei] << " p = " << pp[facei] << endl;
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alpha(pp[facei], pT[facei]);
            }
        }
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::heTabularPsiThermo<MixtureType>::heTabularPsiThermo(const fvMesh& mesh, const objectRegistry& obj)
:
    basicPsiThermo(mesh, obj),
    MixtureType(*this, mesh, obj),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->hBoundaryTypes()
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->eBoundaryTypes()
    )



{
    scalarField& hCells = h_.internalField();
    scalarField& eCells = e_.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(hCells, celli)
    {
      hCells[celli] = this->cellMixture(celli).H(pCells[celli], TCells[celli]);
    }

    forAll(eCells, celli)
    {
      eCells[celli] = this->cellMixture(celli).E(pCells[celli], TCells[celli]);
    }


    forAll(h_.boundaryField(), patchi)
    {
        h_.boundaryField()[patchi] ==
            h(
              this->p_.boundaryField()[patchi],
              this->T_.boundaryField()[patchi],
              patchi
              );
    }

    forAll(e_.boundaryField(), patchi)
    {
        e_.boundaryField()[patchi] ==
            e(
              this->p_.boundaryField()[patchi],
              this->T_.boundaryField()[patchi],
              patchi
              );
    }



    hBoundaryCorrection(h_);
    eBoundaryCorrection(e_);

    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::heTabularPsiThermo<MixtureType>::~heTabularPsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::heTabularPsiThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heTabularPsiThermo<MixtureType>::correct()" << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate2(); // starting use rho and e

    if (debug)
    {
        Info<< "exiting heTabularPsiThermo<MixtureType>::correct()" << endl;
    }
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::h
(
 const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, celli)
    {
      h[celli] = this->cellMixture(cells[celli]).H(p[celli], T[celli]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::h
(
 const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, facei)
    {
      h[facei] = this->patchFaceMixture(patchi, facei).H(p[facei], T[facei]);
    }

    return th;
}

template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::e                 //Add internal energy here
(
 const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> te(new scalarField(T.size()));
    scalarField& e = te();

    forAll(T, celli)
    {
      e[celli] = this->cellMixture(cells[celli]).E(p[celli], T[celli]);
    }

    return te;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::e
(
 const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> te(new scalarField(T.size()));
    scalarField& e = te();

    forAll(T, facei)
    {
      e[facei] = this->patchFaceMixture(patchi, facei).E(p[facei], T[facei]);
    }

    return te;
}








template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::Cp
(
 const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, facei)
    {
      cp[facei] = this->patchFaceMixture(patchi, facei).Cp(p[facei], T[facei]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::Cp
(
 const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, celli)
    {
      cp[celli] = this->cellMixture(cells[celli]).Cp(p[celli], T[celli]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heTabularPsiThermo<MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();

    forAll(this->T_, celli)
    {
      cp[celli] = this->cellMixture(celli).Cp(this->p_[celli], this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];

        forAll(pT, facei)
        {
          pCp[facei] = this->patchFaceMixture(patchi, facei).Cp(pp[facei], pT[facei]);
        }
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heTabularPsiThermo<MixtureType>::Cv
(
 const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv();

    forAll(T, facei)
    {
      cv[facei] = this->patchFaceMixture(patchi, facei).Cv(p[facei], T[facei]);
    }

    return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heTabularPsiThermo<MixtureType>::Cv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cv = tCv();

    forAll(this->T_, celli)
    {
      cv[celli] = this->cellMixture(celli).Cv(this->p_[celli], this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] =
            Cv(
               this->p_.boundaryField()[patchi],
               this->T_.boundaryField()[patchi],
               patchi
               );
    }

    return tCv;
}


template<class MixtureType>
bool Foam::heTabularPsiThermo<MixtureType>::read()
{
    if (basicPsiThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
