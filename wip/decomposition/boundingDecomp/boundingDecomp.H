/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::boundingDecomp

Description

SourceFiles
    boundingDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef boundingDecomp_H
#define boundingDecomp_H

#include "geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boundingDecomp Declaration
\*---------------------------------------------------------------------------*/

class boundingDecomp
:
    public geomDecomp
{
    // Private Member Functions

        labelList decomposeOneProc(const pointField& points) const;

        labelList decomposeOneProc
        (
            const pointField& points,
            const scalarField& weights
        ) const;

        //- Disallow default bitwise copy construct and assignment
        void operator=(const boundingDecomp&);
        boundingDecomp(const boundingDecomp&);


public:

    //- Runtime type information
    TypeName("bounding");


    // Constructors

        //- Construct given the decomposition dictionary
        boundingDecomp(const dictionary& decompositionDict);


    //- Destructor
    virtual ~boundingDecomp()
    {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // simpleDecomp sends all points to the master which does
            // the decomposition.
            return false;
        }

        virtual labelList decompose(const pointField&)
        {
            Info << "not implemented" << endl;
            return labelList();
        }

        virtual labelList decompose(const pointField&, const scalarField&)
        {
            Info << "not implemented" << endl;
            return labelList();
        }

        virtual labelList decompose
        (
            const polyMesh& polyMesh,
            const pointField& points
        );

        virtual labelList decompose
        (
            const polyMesh& polyMesh,
            const pointField& points,
            const scalarField& pointWeights
        )
        {
            Info << "virtual labelList decompose(const polyMesh&, const pointField& points, const scalarField& pointWeights)" << endl;
            return decompose(polyMesh, points);
        }

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights,
            labelList& bCells
        );

        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        )
        {
            Info << "boudningDecomp: should not get called." << endl;
            return labelList::null();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
