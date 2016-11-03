/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "symGaussSeidelSmootherKNL.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symGaussSeidelSmootherKNL, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<symGaussSeidelSmootherKNL>
        addsymGaussSeidelSmootherKNLSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<symGaussSeidelSmootherKNL>
        addsymGaussSeidelSmootherKNLAsymMatrixConstructorToTable_;
}


namespace {
    
inline void gs_impl(
    const Foam::label len,
    const Foam::label fStart, 
    const Foam::label celli,
    const Foam::label* const uPtr,
    const Foam::scalar* const upperPtr,
    const Foam::scalar* const lowerPtr,
    const Foam::scalar* const diagRcpPtr,
    Foam::scalar* psiPtr,
    Foam::scalar* bPrimePtr
) {
    
    // Get the accumulated neighbour side
    Foam::scalar psii = bPrimePtr[celli];
    
    // Accumulate the owner product side
    for (Foam::label i=0; i<len; i++)
    {
        psii -= upperPtr[fStart+i]*psiPtr[uPtr[fStart+i]];
    }

    // Finish psi for this cell
    psii *= diagRcpPtr[celli];

    // Distribute the neighbour side using psi for this cell
    #pragma ivdep
    for (Foam::label i=0; i<len; i++)
    {
        bPrimePtr[uPtr[fStart+i]] -= lowerPtr[fStart+i]*psii;
    }

    psiPtr[celli] = psii;
}

template<int N>
inline void gs_impl_tpl(
    const Foam::label fStart, 
    const Foam::label celli,
    const Foam::label* const uPtr,
    const Foam::scalar* const upperPtr,
    const Foam::scalar* const lowerPtr,
    const Foam::scalar* const diagRcpPtr,
    Foam::scalar* psiPtr,
    Foam::scalar* bPrimePtr
) {
    gs_impl(N, fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr);
}

} // end of private namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::symGaussSeidelSmootherKNL::symGaussSeidelSmootherKNL
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    rD_(matrix.diag().size())
{
    const label nCells = rD_.size();
    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    scalar* __restrict__ rDPtr = rD_.begin();
    for (label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/diagPtr[cell];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::symGaussSeidelSmootherKNL::smooth
(
    const word& fieldName_,
    scalarField& psi,
    const lduMatrix& matrix_,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs_,
    const lduInterfaceFieldPtrsList& interfaces_,
    const direction cmpt,
    const label nSweeps
) const
{
    scalar* __restrict__ psiPtr = psi.begin();

    const label nCells = psi.size();

    scalarField bPrime(nCells);
    scalar* __restrict__ bPrimePtr = bPrime.begin();

    const scalar* const __restrict__ diagRcpPtr = rD_.begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    FieldField<Field, scalar>& mBouCoeffs =
        const_cast<FieldField<Field, scalar>&>
        (
            interfaceBouCoeffs_
        );

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }


    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        bPrime = source;

        matrix_.initMatrixInterfaces
        (
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        matrix_.updateMatrixInterfaces
        (
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        label fStart;
        label fEnd = ownStartPtr[0];

        for (label celli=0; celli<nCells; celli++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = ownStartPtr[celli + 1];

            const label fLen = fEnd - fStart;
            switch(fLen)
            {
            case 0: gs_impl_tpl<0>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 1: gs_impl_tpl<1>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 2: gs_impl_tpl<2>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 3: gs_impl_tpl<3>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 4: gs_impl_tpl<4>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 5: gs_impl_tpl<5>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 6: gs_impl_tpl<6>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 7: gs_impl_tpl<7>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 8: gs_impl_tpl<8>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 9: gs_impl_tpl<9>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 10: gs_impl_tpl<10>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 11: gs_impl_tpl<11>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 12: gs_impl_tpl<12>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 13: gs_impl_tpl<13>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 14: gs_impl_tpl<14>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 15: gs_impl_tpl<15>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 16: gs_impl_tpl<16>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            default: gs_impl(fLen, fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            }
        }

        fStart = ownStartPtr[nCells];

        for (label celli=nCells-1; celli>=0; celli--)
        {
            // Start and end of this row
            fEnd = fStart;
            fStart = ownStartPtr[celli];

            const label fLen = fEnd - fStart;
            switch(fLen)
            {
            case 0: gs_impl_tpl<0>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 1: gs_impl_tpl<1>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 2: gs_impl_tpl<2>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 3: gs_impl_tpl<3>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 4: gs_impl_tpl<4>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 5: gs_impl_tpl<5>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 6: gs_impl_tpl<6>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 7: gs_impl_tpl<7>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 8: gs_impl_tpl<8>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 9: gs_impl_tpl<9>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 10: gs_impl_tpl<10>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 11: gs_impl_tpl<11>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 12: gs_impl_tpl<12>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 13: gs_impl_tpl<13>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 14: gs_impl_tpl<14>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 15: gs_impl_tpl<15>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            case 16: gs_impl_tpl<16>(fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            default: gs_impl(fLen, fStart, celli, uPtr, upperPtr, lowerPtr, diagRcpPtr, psiPtr, bPrimePtr); break;
            }
        }
    }

    // Restore interfaceBouCoeffs_
    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }
}


void Foam::symGaussSeidelSmootherKNL::smooth
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        fieldName_,
        psi,
        matrix_,
        source,
        interfaceBouCoeffs_,
        interfaces_,
        cmpt,
        nSweeps
    );
}


// ************************************************************************* //
